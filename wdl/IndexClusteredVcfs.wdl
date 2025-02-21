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
    # Usage: IndexVcf.bash <VCF>
    # Index a VCF using tabix

    set -o errexit
    set -o nounset
    set -o pipefail

    vcf='~{vcf}'
    vcf_bn="$(basename "${vcf}")"

    # The index file will be placed in the current directory so that it will be
    # delocalized to the top level task output directory.
    printf 'indexing %s\n' "${vcf}" >&2
    tabix "${vcf}"
    mv "${vcf}.tbi" "${vcf_bn}.tbi"
    printf 'done\n' >&2
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
