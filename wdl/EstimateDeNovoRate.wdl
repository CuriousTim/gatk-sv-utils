version 1.0

import "GatherGenotypes.wdl" as gg

workflow EstimateDeNovoRate {
  input {
    Array[File] vcfs
    Array[File] vcf_indicies
    File pedigree
    String? vcf_include_filter
    Int? genotypes_per_shard
    String runtime_docker
  }

  parameter_meta {
    vcfs: "The VCFs are usually split by contig, but not required."
  }

  call GatherTrios {
    input:
      pedigree = pedigree,
      runtime_docker = runtime_docker
  }

  scatter (i in range(length(vcfs))) {
    call gg.GatherGenotypes {
      input:
        vcf = vcfs[i],
        vcf_index = vcf_indicies[i],
        vcf_include_filter = vcf_include_filter,
        genotypes_per_shard = genotypes_per_shard,
        runtime_docker = runtime_docker
    }

    call CountDeNovos {
      input:
        genotypes = GatherGenotypes.genotypes,
        trios = GatherTrios.trios_output,
        output_prefix = i,
        runtime_docker = runtime_docker
    }
  }

  call MergeDeNovoCounts {
    input:
      counts_dbs = CountDeNovos.counts_db,
      runtime_docker = runtime_docker
  }

  output {
    File denovo_counts = MergeDeNovoCounts.merged_denovo_counts
    File merged_db = MergeDeNovoCounts.merged_db
  }
}

task GatherTrios {
  input {
    File pedigree
    String runtime_docker

    Float? memory_gib
    Int? boot_disk_gb
    Int? cpus
    Int? disk_gb
    Int? max_retries
    Int? preemptible_tries
  }

  Float disk_size = size(pedigree, "GB") * 1.5 + 16
  String trios = "trios.tsv"

  runtime {
    bootDiskSizeGb: select_first([boot_disk_gb, 8])
    cpus: select_first([cpus, 1])
    disks: "local-disk ${select_first([disk_gb, ceil(disk_size)])} HDD"
    docker: runtime_docker
    maxRetries: select_first([max_retries, 1])
    memory: "${select_first([memory_gib, 1])} GiB"
    preemptible: select_first([preemptible_tries, 3])
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    # Missing parents are usually represented as a 0-length string or the literal
    # 0, both of which are treated as false by awk.
    gawk -F'\t' '!/^#/ && $2 && $3 && $4 {print $2"\t"$3"\t"$4}' '~{pedigree}' > temp

    read -r trio_n _ < <(wc -l temp)
    if (( trio_n == 0 )); then
      printf 'no trios found\n' >&2
      exit 1
    else
      printf 'found %d trios\n' "${trio_n}" >&2
    fi

    sort -u -k 1,1 temp > '~{trios}'
    read -r uniq_trio_n _ < <(wc -l '~{trios}')
    if (( uniq_trio_n != trio_n )); then
      printf 'warning: trios (%d) does not equal unique trios (%d)\n' "${trio_n}" "${uniq_trio_n}" >&2
    fi
  >>>

  output {
    File trios_output = "${trios}"
  }
}

task CountDeNovos {
  input {
    Array[File] genotypes
    File trios
    String output_prefix
    String runtime_docker

    Float? memory_gib
    Int? boot_disk_gb
    Int? cpus
    Int? disk_gb
    Int? max_retries
    Int? preemptible_tries
  }

  # Need space to make genotypes database and counts database.
  Float disk_size = size(genotypes, "GB") * 2.5 + size(trios, "GB") + 16
  Int min_cpus = 2
  Int max_cpus = 8
  # The choice of 1 CPU per GB of genotypes assumes that the genotypes files have a fixed
  # number of genotypes to file size ratio.
  Int genotypes_size = floor(size(genotypes, "GB"))
  Int predicted_cpus = if (genotypes_size * 2) < min_cpus then min_cpus else (genotypes_size * 2)
  Int default_cpus = if predicted_cpus > max_cpus then max_cpus else predicted_cpus
  # DuckDB recommends 3-4 GB of memory per thread for join-heavy workloads.
  Int default_mem = default_cpus * 4

  String output_db = "${output_prefix}-dn_counts.duckdb"

  runtime {
    bootDiskSizeGb: select_first([boot_disk_gb, 8])
    cpus: select_first([cpus, default_cpus])
    disks: "local-disk ${select_first([disk_gb, ceil(disk_size)])} SSD"
    docker: runtime_docker
    maxRetries: select_first([max_retries, 1])
    memory: "${select_first([memory_gib, default_mem])} GiB"
    preemptible: select_first([preemptible_tries, 3])
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    genotypes='~{write_lines(genotypes)}'
    trios='~{trios}'
    output_db='~{output_db}'

    # Need to create a temp directory instead of a file because DuckDB will fail if
    # given an existing file that is not a DuckDB database, even if the file is
    # empty.
    tmpdir="$(mktemp -d -p "${PWD}" tmp_XXXXXXXXX)"
    genotypes_db="${tmpdir}/genotypes.duckdb"
    trap 'rm -rf "${tmpdir}"' EXIT

    genotypes_dir="${tmpdir}/links"
    mkdir "${genotypes_dir}"
    cat "${genotypes}" | xargs -L 1 ln -s -t "${genotypes_dir}"

    duckdb -bail -echo "${genotypes_db}" <<EOF
      PRAGMA disable_progress_bar;
      CREATE TYPE svclass AS ENUM ('DEL', 'DUP', 'INS', 'INV', 'CPX', 'CTX', 'CNV', 'BND');
      CREATE TABLE genotypes (
        sid VARCHAR,
        vid VARCHAR,
        svtype svclass,
        svlen INTEGER,
        gt UTINYINT,
        gq UTINYINT,
      );
      COPY genotypes FROM '${genotypes_dir}/*.tsv.zst' (DELIM '\t', NULLSTR '.');
    EOF

    duckdb -bail -echo "${output_db}" <<EOF
      PRAGMA disable_progress_bar;
      CREATE TYPE svclass AS ENUM ('DEL', 'DUP', 'INS', 'INV', 'CPX', 'CTX', 'CNV', 'BND');
      CREATE TABLE trios (
        offspring VARCHAR,
        father VARCHAR,
        mother VARCHAR
      );
      COPY trios FROM '${trios}' (DELIM '\t', HEADER false);

      CREATE TABLE inheritance (
        offspring VARCHAR,
        father VARCHAR,
        mother VARCHAR,
        vid VARCHAR,
        svtype svclass,
        svlen INTEGER,
        of_gt UTINYINT,
        of_gq UTINYINT,
        fa_gt UTINYINT,
        fa_gq UTINYINT,
        mo_gt UTINYINT,
        mo_gq UTINYINT
      );
      ATTACH '${genotypes_db}' AS gts_db;
      INSERT INTO inheritance BY NAME
      SELECT l.offspring, l.father, l.mother, r.vid, r.svtype, r.svlen, r.gt AS of_gt, r.gq AS of_gq
      FROM trios l
      JOIN (SELECT sid, vid, svtype, svlen, gt, gq FROM gts_db.genotypes WHERE gt > 0) r
      ON (l.offspring = r.sid);

      UPDATE inheritance l
      SET fa_gt = r.gt, fa_gq = r.gq
      FROM gts_db.genotypes r
      WHERE l.father = r.sid AND l.vid = r.vid;

      UPDATE inheritance l
      SET mo_gt = r.gt, mo_gq = r.gq
      FROM gts_db.genotypes r
      WHERE l.mother = r.sid AND l.vid = r.vid;
    EOF
  >>>

  output {
    File counts_db = output_db
  }
}

task MergeDeNovoCounts {
  input {
    Array[File] counts_dbs
    String runtime_docker

    Float? memory_gib
    Int? boot_disk_gb
    Int? cpus
    Int? disk_gb
    Int? max_retries
    Int? preemptible_tries
  }

  # All the databases are merged into one so we need double the space.
  Float disk_size = size(counts_dbs, "GB") * 2.5 + 16

  runtime {
    bootDiskSizeGb: select_first([boot_disk_gb, 8])
    cpus: select_first([cpus, 2])
    disks: "local-disk ${select_first([disk_gb, ceil(disk_size)])} HDD"
    docker: runtime_docker
    maxRetries: select_first([max_retries, 1])
    memory: "${select_first([memory_gib, 4])} GiB"
    preemptible: select_first([preemptible_tries, 3])
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    dbs_list='~{write_lines(counts_dbs)}'
    merged_counts='merged_denovo_counts.tsv.gz'
    merged_db='inheritance.duckdb'

    sql_cmds="$(mktemp -p "${PWD}" tmp_XXXXXXXXX)"
    tmp_merge="$(mktemp -p "${PWD}" tmp_XXXXXXXXX)"
    tmp_db_dir="$(mktemp -d -p "${PWD}" tmp_XXXXXXXXX)"
    trap 'rm -rf "${sql_cmds}" "${tmp_merge}" "${tmp_db_dir}"' EXIT

    declare -i i=0
    while read -r f; do
      dest="${tmp_db_dir}/inheritance-$(printf '%02d' "${i}").duckdb"
      ln -s "${f}" "${dest}"
      printf "ATTACH '%s' AS db_%d;\n" "${dest}" "${i}" >> "${sql_cmds}"
      (( ++i ))
    done < "${dbs_list}"

    declare -i j=0
    printf 'CREATE TABLE merged_inh AS SELECT * FROM (\n' >> "${sql_cmds}"
    for (( j=0; j<i - 1; ++j )); do
      printf 'SELECT * FROM db_%d.inheritance UNION ALL\n' "${j}" >> "${sql_cmds}"
    done
    printf 'SELECT * FROM db_%d.inheritance);\n' "${j}" >> "${sql_cmds}"

    cat >> "${sql_cmds}" <<EOF
    COPY (
      SELECT
        offspring,
        svtype,
        count() FILTER (of_gt > 0) AS sv_count,
        count() FILTER (of_gt = 1 AND fa_gt = 0 AND mo_gt = 0) AS dn_count
      FROM merged_inh
      GROUP BY offspring, svtype
    ) TO '${tmp_merge}' (DELIM '\t', HEADER true, FORMAT 'CSV', COMPRESSION 'gzip');
    EOF

    duckdb -bail -echo "${merged_db}" < "${sql_cmds}"
    mv "${tmp_merge}" "${merged_counts}"
  >>>

  output {
    File merged_denovo_counts = "merged_denovo_counts.tsv.gz"
    File merged_db = "inheritance.duckdb"
  }
}
