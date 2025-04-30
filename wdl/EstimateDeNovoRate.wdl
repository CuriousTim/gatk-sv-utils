version 1.0

import "GatherGenotypes.wdl"  as gg

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
  Int predicted_cpus = if genotypes_size < min_cpus then min_cpus else genotypes_size
  Int default_cpus = if predicted_cpus > max_cpus then max_cpus else predicted_cpus
  # DuckDB recommends 3-4 GB of memory per thread for join-heavy workloads.
  Int default_mem = default_cpus * 3

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
    /opt/task_scripts/EstimateDeNovoRate/CountDeNovos '~{write_lines(genotypes)}' \
      '~{trios}' '~{output_db}'
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
