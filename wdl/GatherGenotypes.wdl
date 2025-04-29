version 1.0

workflow GatherGenotypes {
  input {
    File vcf
    File vcf_index
    String vcf_include_filter = 'FILTER == "PASS" && SVTYPE != "CNV" && SVTYPE != "BND"'
    Int genotypes_per_shard = 100000000
    String runtime_docker
  }

  parameter_meta {
    vcf: "BCF file works too."
    vcf_index: "The index must have per contig stats, i.e. bcftools index --stats must work."
    vcf_include_filter: "Only used when gathering genotypes, not splitting VCF samples."
    genotype_per_shard: "Number of genotypes per shard. This approximate because divisibility and site filters."
  }

  call SplitVcfSamples {
    input:
      vcf = vcf,
      vcf_index = vcf_index,
      genotypes_per_shard = genotypes_per_shard,
      runtime_docker = runtime_docker
  }

  scatter (i in range(length(SplitVcfSamples.shards))) {
    call GatherGenotypesForSamples {
      input:
        vcf = vcf,
        samples = SplitVcfSamples.shards[i],
        output_prefix = i,
        vcf_include_filter = vcf_include_filter,
        runtime_docker = runtime_docker
    }
  }

  call TarGenotypes {
    input:
      genotypes = GatherGenotypesForSamples.genotypes,
      runtime_docker = runtime_docker
  }

  output {
    File genotypes_tar = TarGenotypes.genotypes_tar
  }
}

task SplitVcfSamples {
  input {
    File vcf
    File vcf_index
    Int genotypes_per_shard
    String runtime_docker

    Float? memory_gib
    Int? boot_disk_gb
    Int? cpus
    Int? disk_gb
    Int? max_retries
    Int? preemptible_tries
  }

  Float disk_size = size(vcf, "GB") + 16

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
    /opt/task_scripts/GatherGenotypes/SplitVcfSamples '~{vcf}' \
      '~{genotypes_per_shard}' 'shards'
  >>>

  output {
    Array[File] shards = glob("shards/*.list")
  }
}

task GatherGenotypesForSamples {
  input {
    File vcf
    File samples
    Int output_prefix
    String vcf_include_filter
    String runtime_docker

    Float? memory_gib
    Int? boot_disk_gb
    Int? cpus
    Int? disk_gb
    Int? max_retries
    Int? preemptible_tries
  }

  Float disk_size = size(vcf, "GB") * 2 + 16
  String output_file = "${output_prefix}-genotypes.tsv.zst"

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
    /opt/task_scripts/GatherGenotypes/GatherGenotypesForSamples '~{vcf}' \
      '~{vcf_include_filter}' '~{samples}' '~{output_file}'
  >>>

  output {
    File genotypes = output_file
  }
}

task TarGenotypes {
  input {
    Array[File] genotypes
    String runtime_docker

    Float? memory_gib
    Int? boot_disk_gb
    Int? cpus
    Int? disk_gb
    Int? max_retries
    Int? preemptible_tries
  }

  # The genotype files are copied to the tar directory.
  Float disk_size = size(genotypes, "GB") * 2 + 16

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
    /opt/task_scripts/GatherGenotypes/TarGenotypes '~{write_lines(genotypes)}' \
      'genotypes'
  >>>

  output {
    File genotypes_tar = 'genotypes.tar'
  }
}
