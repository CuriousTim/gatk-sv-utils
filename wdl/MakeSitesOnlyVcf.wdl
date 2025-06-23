version 1.0

# Drop all sample-level information from a VCF
workflow MakeSitesOnlyVcf {
  input {
    Array[File] vcfs
    String runtime_docker
  }

  scatter (vcf in vcfs) {
    call DropSamplesFromVcf {
      input:
        vcf = vcf,
        runtime_docker = runtime_docker
    }
  }

  output {
    Array[File] sites_only_vcf = DropSamplesFromVcf.sites_only_vcf
    Array[File] sites_only_vcf_index = DropSamplesFromVcf.sites_only_vcf_index
  }
}

task DropSamplesFromVcf {
  input {
    File vcf
    String runtime_docker

    Float? memory_gib
    Int? boot_disk_gb
    Int? cpus
    Int? disk_gb
    Int? max_retries
    Int? preemptible_tries
  }

  Float disk_size = size(vcf, "GB") * 2 + 16
  String output_vcf = "sites_only-${basename(vcf)}"
  String output_vcf_index = "${output_vcf}.tbi"

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

    in_vcf='~{vcf}'
    out_vcf='~{output_vcf}'

    bcftools view --drop-genotypes --output-type z --output "${out_vcf}" \
      --write-index=tbi --no-update "${in_vcf}"
  >>>

  output {
    File sites_only_vcf = output_vcf
    File sites_only_vcf_index = output_vcf_index
  }
}
