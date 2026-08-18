version 1.0

# Inspect PE/SR/RD evidence to verify SV de novo status
workflow CheckDeNovoEvidence {
  input {
    # The variants file can be a TSV with the following format:
    # The file has a header line as the first line and there must be columns named chr, start,
    # end, svlen, vid, svtype, and sample. These columns can be in any order relative to each
    # other and there can be other columns.
    File variants
    Int variants_per_batch = 100

    # TSV with sample ID, sample set ID
    File sample_table
    File pedigree
    Array[File] merged_pe
    Array[File] merged_pe_index
    Array[File] merged_sr
    Array[File] merged_sr_index
    Array[File] merged_rd
    Array[File] merged_rd_index
    Array[File] median_cov
    Array[String] sample_set_id

    String output_prefix

    String base_docker
    String r_docker
  }

  output {
    File checked_variants = MergeCheckedVariants.merged_checked_variants
  }

  call BatchVariants {
    input:
      variants = variants,
      variants_per_batch = variants_per_batch,
      base_docker = base_docker
  }

  call MakeEvidenceManifest {
    input:
      sample_set_id = sample_set_id,
      merged_pe = merged_pe,
      merged_sr = merged_sr,
      merged_rd = merged_rd,
      median_cov = median_cov,
      base_docker = base_docker
  }

  scatter (i in range(length(BatchVariants.batched_variants))) {
    if (size(BatchVariants.batched_variants[i]) > 0) {
      call CheckVariants {
        input:
          variants = BatchVariants.batched_variants[i],
          evidence_manifest = MakeEvidenceManifest.evidence_manifest,
          output_prefix = "batch_${i}",
          sample_table = sample_table,
          pedigree = pedigree,
          r_docker = r_docker,
          max_svlen_file = BatchVariants.max_svlens[i]
      }
    }
  }

  call MergeCheckedVariants {
    input:
      variants = select_all(CheckVariants.checked_variants),
      output_prefix = output_prefix,
      base_docker = base_docker
  }
}

task BatchVariants {
  input {
    File variants
    Int variants_per_batch
    String base_docker
  }

  output {
    Array[File] batched_variants = glob("batches/*.tsv")
    Array[File] max_svlens = glob("mems/*")
  }

  Float disk_size = size(variants, "GB") * 2 + 16

  runtime {
    bootDiskSizeGb: 8
    cpu: 1
    disks: "local-disk ${ceil(disk_size)} HDD"
    docker: base_docker
    maxRetries: 1
    memory: "3.75 GiB"
    preemptible: 3
    noAddress: true
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    variants='~{variants}'
    variants_per_batch='~{variants_per_batch}'

    mv "${variants}" variants.tsv.gz
    # sort by SV length so that SVs of similar size will end up in the same VM
    duckdb ':memory:' \
      "COPY (SELECT * FROM read_csv('variants.tsv.gz', delim = '\t') ORDER BY svlen) TO 'temp.tsv' (FORMAT CSV, DELIM '\t');"

    mkdir batches
    mkdir mems
    gawk -F'\t' -v n="${variants_per_batch}" '
      BEGIN { OFS = "\t" }
      FNR == 1 {
        header = $0
        for (i = 1; i <= NF; ++i) {
          if ($i == "svlen") {
            svlen_field = i
          }
        }
        if (!svlen_field) {
          print "variants file does not have \"svlen\" columns" > "/dev/stderr"
          exit 1
        }
        next
      }
      j % n == 0 {
        if (outpath) {
          close(outpath)
          mem_outpath = sprintf("mems/%06d", k - 1)
          print max_svlen > mem_outpath
          close(mem_outpath)
          max_svlen = 0
        }
        outpath = sprintf("batches/%06d.tsv", k++)
        print header > outpath
      }
      {
        print > outpath
        max_svlen = max_svlen >= $svlen_field ? max_svlen : $svlen_field
        ++j
      }
      END {
        if (k > 0) {
          mem_outpath = sprintf("mems/%06d", k - 1)
          print max_svlen > mem_outpath
          close(mem_outpath)
        }
      }
    ' temp.tsv
  >>>
}

task MakeEvidenceManifest {
  input {
    Array[String] sample_set_id
    Array[String] merged_pe
    Array[String] merged_sr
    Array[String] merged_rd
    Array[String] median_cov
    String base_docker
  }

  output {
    File evidence_manifest = "evidence_manifest.tsv"
  }

  runtime {
    bootDiskSizeGb: 8
    cpu: 1
    disks: "local-disk 64 HDD"
    docker: base_docker
    maxRetries: 1
    memory: "3.75 GiB"
    preemptible: 3
    noAddress: true
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    sample_set_id='~{write_lines(sample_set_id)}'
    merged_pe='~{write_lines(merged_pe)}'
    merged_sr='~{write_lines(merged_sr)}'
    merged_rd='~{write_lines(merged_rd)}'
    median_cov='~{write_lines(median_cov)}'

    paste "${sample_set_id}" "${merged_pe}" "${merged_sr}" "${merged_rd}" \
      "${median_cov}" > 'evidence_manifest.tsv'
  >>>
}

task CheckVariants {
  input {
    File variants
    File evidence_manifest
    File sample_table
    String output_prefix
    File pedigree
    String r_docker
    File max_svlen_file

    Float? mem_gib
    Int? disk_gb
    Int? cpus
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
  }

  output {
    File checked_variants = "${output_prefix}.tsv.gz"
  }

  Int max_svlen = read_int(max_svlen_file)
  Int default_cpus = if max_svlen >= 10000000 then if max_svlen >= 50000000 then 16 else 4 else 2
  Float default_mem_gib = if max_svlen >= 10000000 then if max_svlen >= 50000000 then 128 else 32 else 16

  runtime {
    bootDiskSizeGb: select_first([boot_disk_gb, 8])
    cpu: select_first([cpus, default_cpus])
    disks: "local-disk " + select_first([disk_gb, 128]) + " HDD"
    docker: r_docker
    maxRetries: select_first([max_retries, 1])
    memory: select_first([mem_gib, default_mem_gib]) + " GiB"
    preemptible: select_first([preemptible_tries, 3])
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    variants="~{variants}"
    evidence_manifest="~{evidence_manifest}"
    sample_table="~{sample_table}"
    output_prefix="~{output_prefix}"
    pedigree="~{pedigree}"

    Rscript /opt/gatk-sv-utils/scripts/check_sv_evidence.R \
      "${variants}" "${pedigree}" "${evidence_manifest}" "${sample_table}" \
      "${output_prefix}.tsv.gz"
  >>>
}

task MergeCheckedVariants {
  input {
    Array[File] variants
    String output_prefix
    String base_docker
  }

  output {
    File merged_checked_variants = "${output_prefix}.tsv.gz"
  }

  Float disk_size = size(variants, "GB") * 2 + 32

  runtime {
    bootDiskSizeGb: 8
    cpu: 1
    disks: "local-disk ${ceil(disk_size)} HDD"
    docker: base_docker
    maxRetries: 1
    memory: "3.75 GiB"
    preemptible: 3
    noAddress: true
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    variants="~{write_lines(variants)}"
    output_prefix="~{output_prefix}"

    x=1
    while read -r f; do
      if (( x == 1 )); then
        gzip -cd "${f}"
        x=0
      else
        gzip -cd "${f}" | gawk 'NR>1'
      fi
    done < "${variants}" \
      | gzip -c > "${output_prefix}.tsv.gz"
  >>>
}
