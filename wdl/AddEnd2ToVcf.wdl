version 1.0

# Add END2 field to BND and CTX variants in VCF
workflow AddEnd2ToVcf {
  input {
    Array[File] vcfs
    String runtime_docker
  }

  scatter (f in vcfs) {
    call AddEnd2 {
      input:
        vcf = f,
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
    Int? disk_gb
    Int? cpus
    Int? preemptible_tries
    Int? max_retries
    Int? boot_disk_gb
  }

  Float disk_size = size(vcf, "GB") * 2.2 + 16
  String output_vcf = "end2-${basename(vcf)}"
  String output_vcf_index = "${output_vcf}.tbi"

  runtime {
    memory: "${select_first([memory_gib, 2])} GiB"
    disks: "local-disk ${select_first([disk_gb, ceil(disk_size)])} HDD"
    cpus: select_first([cpus, 1])
    preemptible: select_first([preemptible_tries, 3])
    docker: runtime_docker
    bootDiskSizeGb: select_first([boot_disk_gb, 16])
  }

  command <<<
    bash /opt/task-scripts/AddEnd2ToVcf/AddEnd2.bash '~{vcf}' '~{output_vcf}'
  >>>

  output {
    File modified_vcf = "${output_vcf}"
    File modified_vcf_index = "${output_vcf_index}"
  }
}
