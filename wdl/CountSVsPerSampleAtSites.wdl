version 1.0

workflow CountSVsPerSampleAtSites {
  input {
    Array[File] vcfs
    File variant_ids
    String base_docker
  }

  scatter (vcf in vcfs) {
    call CountSVsPerSample {
      input:
        vcf = vcf,
        variant_ids = variant_ids,
        base_docker = base_docker
    }
  }

  call MergeCounts {
    input:
      counts = CountSVsPerSample.counts,
      base_docker = base_docker
  }

  output {
    File counts = MergeCounts.merged_counts
  }
}

task CountSVsPerSample {
  input {
    File vcf
    File variant_ids
    String base_docker

    Float? memory_gib
    Int? disk_gb
    Int? cpus
    Int? preemptible_tries
    Int? max_retries
    Int? boot_disk_gb
  }

  Float disk_size = size(vcf, "GB") * 2 + 50

  runtime {
    bootDiskSizeGb: select_first([boot_disk_gb, 8])
    cpus: select_first([cpus, 1])
    disks: "local-disk ${select_first([disk_gb, ceil(disk_size)])} HDD"
    docker: base_docker
    maxRetries: select_first([max_retries, 1])
    memory: "${select_first([memory_gib, 4])} GiB"
    preemptible: select_first([preemptible_tries, 3])
  }

  String counts_name = basename(vcf, ".vcf.gz") + "-counts.tsv.gz"

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    vcf='~{vcf}'
    variant_ids='~{variant_ids}'
    counts_name='~{counts_name}'

    mv "${variant_ids}" variants

    bcftools query --include 'ID=@variants' --format '[%SAMPLE\t%GT\n]' "${vcf}" \
      | gawk '$2 ~ /\./{next} {split($2, a, "/"); b[$1]+=a[1] + 0 + a[2]} END{for(i in b){print i "\t" b[i]}}' \
      | gzip -c > "${counts_name}"
  >>>

  output {
    File counts = counts_name
  }
}

task MergeCounts {
  input {
    Array[File] counts
    String base_docker

    Float? memory_gib
    Int? disk_gb
    Int? cpus
    Int? preemptible_tries
    Int? max_retries
    Int? boot_disk_gb
  }

  Float disk_size = size(counts, "GB") * 2 + 16

  runtime {
    bootDiskSizeGb: select_first([boot_disk_gb, 8])
    cpus: select_first([cpus, 1])
    disks: "local-disk ${select_first([disk_gb, ceil(disk_size)])} SSD"
    docker: base_docker
    maxRetries: select_first([max_retries, 1])
    memory: "${select_first([memory_gib, 4])} GiB"
    preemptible: select_first([preemptible_tries, 3])
  }

  String merged_counts_name = "svs_per_sample.tsv.gz"

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    counts_files='~{write_lines(counts)}'
    merged_counts_name='~{merged_counts_name}'

    mkdir counts

    i=0
    while read -r p; do
      mv "${p}" "counts/${i}.tsv.gz"
      i=$((i + 1))
    done < "${counts_files}"

    duckdb -bail ':memory:' "COPY (SELECT sid, sum(count) AS svs_per_sample FROM read_csv('counts/*.tsv.gz', names = ['sid', 'count'], ignore_errors = true, union_by_name = true) GROUP BY sid) TO 'svs_per_sample.tsv.gz' (DELIMITER '\\t');"
  >>>

  output {
    File merged_counts = merged_counts_name
  }
}
