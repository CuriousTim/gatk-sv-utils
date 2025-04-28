version 1.0

workflow EstimateDeNovoRate {
  input {
    Array[File] vcfs
    File pedigree
    String? vcf_filter
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
        vcf_filter = vcf_filter,
        runtime_docker = runtime_docker
    }
  }

  output {
    File trios = GatherTrios.trios_output
    Array[File] genotypes = GatherGenotypes.genotypes_tsv
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
    String vcf_filter = 'FILTER == "PASS" && SVTYPE != "CNV"'
    String runtime_docker

    Float? memory_gib
    Int? boot_disk_gb
    Int? cpus
    Int? disk_gb
    Int? max_retries
    Int? preemptible_tries
  }

  # Need space for the intermediate compressed genotypes file and the DuckDB file.
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
    /opt/task_scripts/EstimateDeNovoRate/GatherGenotypes '~{vcf}' '~{vcf_filter}' '~{genotypes}'
  >>>

  output {
    File genotypes_tsv = genotypes
  }
}
