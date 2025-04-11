version 1.0

# Index the clustered VCFs output by 05-ClusterBatch in GATK-SV.
workflow IndexClusteredVcfs {
  input {
    File? clustered_depth_vcf
    File? clustered_manta_vcf
    File? clustered_melt_vcf
    File? clustered_scramble_vcf
    File? clustered_wham_vcf
    String runtime_docker
  }

  if (defined(clustered_depth_vcf)) {
    call IndexVcf as index_depth_vcf {
      input:
        vcf = select_first([clustered_depth_vcf]),
        runtime_docker = runtime_docker
    }
  }

  if (defined(clustered_manta_vcf)) {
    call IndexVcf as index_manta_vcf {
      input:
        vcf = select_first([clustered_manta_vcf]),
        runtime_docker = runtime_docker
    }
  }

  if (defined(clustered_melt_vcf)) {
    call IndexVcf as index_melt_vcf {
      input:
        vcf = select_first([clustered_melt_vcf]),
        runtime_docker = runtime_docker
    }
  }

  if (defined(clustered_wham_vcf)) {
    call IndexVcf as index_wham_vcf {
      input:
        vcf = select_first([clustered_wham_vcf]),
        runtime_docker = runtime_docker
    }
  }

  if (defined(clustered_scramble_vcf)) {
    call IndexVcf as index_scramble_vcf {
      input:
        vcf = select_first([clustered_scramble_vcf]),
        runtime_docker = runtime_docker
    }
  }

  output {
    File? clustered_depth_vcf_index = index_depth_vcf.vcf_index
    File? clustered_manta_vcf_index = index_manta_vcf.vcf_index
    File? clustered_melt_vcf_index = index_melt_vcf.vcf_index
    File? clustered_wham_vcf_index = index_wham_vcf.vcf_index
    File? clustered_scramble_vcf_index = index_scramble_vcf.vcf_index
  }
}

task IndexVcf {
  input {
    File vcf
    String runtime_docker

    Float? memory_gib
    Int? disk_gb
    Int? cpus
    Int? preemptible_tries
    Int? max_retries
    Int? boot_disk_gb
  }

  Float disk_size = size(vcf, "GB") + 16

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
    /opt/task_scripts/IndexClusteredVcfs/IndexVcf '~{vcf}'
  >>>

  output {
    File vcf_index = basename(vcf) + ".tbi"
  }

  meta {
    author: "Justin Lim"
    email: "jlim21@mgh.harvard.edu"
    description: "Index a VCF using tabix."
  }
}
