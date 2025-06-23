version 1.0

# Add END2 to BND and CTX variants in VCF
workflow AddEnd2ToVcf {
  input {
    Array[File] vcfs
    String runtime_docker
  }

  scatter (vcf in vcfs) {
    call AddEnd2 {
      input:
        vcf = vcf,
        runtime_docker = runtime_docker
    }
  }

  output {
    Array[File] modified_vcfs = AddEnd2.modified_vcf
    Array[File] modified_vcf_indicies = AddEnd2.modified_vcf_index
  }
}

task AddEnd2 {
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

  Float disk_size = size(vcf, "GB") * 2.2 + 16
  String output_vcf = "end2-${basename(vcf)}"
  String output_vcf_index = "${output_vcf}.tbi"

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
    set -o errexit
    set -o nounset
    set -o pipefail

    # faster to skip bcftools parsing
    bgzip -cd '~{vcf}' \
      | gawk -f /opt/gatk-sv-utils/scripts/add_end2_to_vcf.awk - \
      | bgzip -c > '~{output_vcf}'
    bcftools index --tbi '~{output_vcf}'
  >>>

  output {
    File modified_vcf = "${output_vcf}"
    File modified_vcf_index = "${output_vcf_index}"
  }
}
