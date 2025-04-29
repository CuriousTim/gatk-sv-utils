version 1.0

workflow EstimateDeNovoRate {
  input {
    Array[File] vcfs
    File pedigree
    String? vcf_include_filter
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

  scatter (vcf in vcfs) {
    call GatherGenotypes {
      input:
        vcf = vcf,
        vcf_include_filter = vcf_include_filter,
        runtime_docker = runtime_docker
    }

    call CountDeNovos {
      input:
        genotypes = GatherGenotypes.genotypes_tsv,
        trios = GatherTrios.trios_output,
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
    File counts_dbs_tar = MergeDeNovoCounts.counts_dbs_tar
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
    /opt/task_scripts/EstimateDeNovoRate/GatherTrios '~{pedigree}' '~{trios}'
  >>>

  output {
    File trios_output = "${trios}"
  }
}

task GatherGenotypes {
  input {
    File vcf
    String vcf_include_filter = 'FILTER == "PASS" && SVTYPE != "CNV" && SVTYPE != "BND"'
    String runtime_docker

    Float? memory_gib
    Int? boot_disk_gb
    Int? cpus
    Int? disk_gb
    Int? max_retries
    Int? preemptible_tries
  }

  # The genotypes file get pretty large even when compressed.
  Float disk_size = size(vcf, "GB") * 4 + 16
  String genotypes = basename(vcf, ".vcf.gz") + "-genotypes.tsv.zst"

  runtime {
    bootDiskSizeGb: select_first([boot_disk_gb, 8])
    cpus: select_first([cpus, 1])
    disks: "local-disk ${select_first([disk_gb, ceil(disk_size)])} HDD"
    docker: runtime_docker
    maxRetries: select_first([max_retries, 1])
    memory: "${select_first([memory_gib, 2])} GiB"
    preemptible: select_first([preemptible_tries, 3])
  }

  command <<<
    /opt/task_scripts/EstimateDeNovoRate/GatherGenotypes '~{vcf}' '~{vcf_include_filter}' '~{genotypes}'
  >>>

  output {
    File genotypes_tsv = genotypes
  }
}

task CountDeNovos {
  input {
    File genotypes
    File trios
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
  String output_db = basename(genotypes, "-genotypes.tsv.zst") + "-dn_counts.duckdb"

  runtime {
    bootDiskSizeGb: select_first([boot_disk_gb, 8])
    cpus: select_first([cpus, 2])
    disks: "local-disk ${select_first([disk_gb, ceil(disk_size)])} SSD"
    docker: runtime_docker
    maxRetries: select_first([max_retries, 1])
    memory: "${select_first([memory_gib, 8])} GiB"
    preemptible: select_first([preemptible_tries, 3])
  }

  command <<<
    /opt/task_scripts/EstimateDeNovoRate/CountDeNovos '~{genotypes}' '~{trios}' '~{output_db}'
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

  # Each database is copied into the tar directory so need double the space.
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
    /opt/task_scripts/EstimateDeNovoRate/MergeDeNovoCounts \
      '~{write_lines(counts_dbs)}' \
      'merged_denovo_counts.tsv.gz' \
      'counts_dbs.tar'
  >>>

  output {
    File merged_denovo_counts = "merged_denovo_counts.tsv.gz"
    File counts_dbs_tar = "counts_dbs.tar"
  }
}
