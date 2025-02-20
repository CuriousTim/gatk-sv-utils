version 1.0

# Index the clustered VCFs output by 05-ClusterBatch in GATK-SV.
workflow IndexClusteredVcfs {
  input {
    File? clustered_depth_vcf
    File? clustered_manta_vcf
    File? clustered_scramble_vcf
    File? clustered_wham_vcf
    File indexvcf_cmd_script
    String runtime_docker
  }

  if (defined(clustered_depth_vcf)) {
    call IndexVcf as index_depth_vcf {
      input:
        vcf = select_first([clustered_depth_vcf, ""]),
        cmd_script = indexvcf_cmd_script,
        runtime_docker = runtime_docker
    }
  }

  if (defined(clustered_manta_vcf)) {
    call IndexVcf as index_manta_vcf {
      input:
        vcf = select_first([clustered_manta_vcf, ""]),
        cmd_script = indexvcf_cmd_script,
        runtime_docker = runtime_docker
    }
  }

  if (defined(clustered_wham_vcf)) {
    call IndexVcf as index_wham_vcf {
      input:
        vcf = select_first([clustered_wham_vcf, ""]),
        cmd_script = indexvcf_cmd_script,
        runtime_docker = runtime_docker
    }
  }

  if (defined(clustered_scramble_vcf)) {
    call IndexVcf as index_scramble_vcf {
      input:
        vcf = select_first([clustered_scramble_vcf, ""]),
        cmd_script = indexvcf_cmd_script,
        runtime_docker = runtime_docker
    }
  }

  output {
    File? clustered_depth_vcf_index = index_depth_vcf.vcf_index
    File? clustered_manta_vcf_index = index_manta_vcf.vcf_index
    File? clustered_wham_vcf_index = index_wham_vcf.vcf_index
    File? clustered_scramble_vcf_index = index_scramble_vcf.vcf_index
  }
}

task IndexVcf {
  input {
    File vcf
    File cmd_script
    String runtime_docker
  }

  runtime {
    memory: "1 GiB"
    disks: "local-disk " + ceil(size(vcf, "GB")) + " HDD"
    cpus: 1
    preemptible: 3
    maxRetries: 1
    docker: runtime_docker
    bootDiskSizeGb: 16
  }

  command <<<
    bash '~{cmd_script}' '~{vcf}'
  >>>

  output {
    File vcf_index = vcf + ".tbi"
  }
}
