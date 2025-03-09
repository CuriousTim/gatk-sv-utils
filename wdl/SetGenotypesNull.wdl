version 1.0

# Set genotypes in VCF to null for given samples on given contigs.
# If no contigs are given, the genotypes for the given samples on all the
# contigs will be set to null.
workflow SetGenotypesNull {
  input {
    Array[File] vcfs
    File samples_list
    Array[String] contigs = ["chrX", "chrY"]

    String runtime_docker
  }

  scatter (vcf in vcfs) {
    call SetGenotypesNullTask {
      input:
        vcf = vcf,
        samples_list = samples_list,
        contigs = contigs,
        runtime_docker = runtime_docker
    }
  }

  output {
    Array[File] nulled_vcfs = SetGenotypesNullTask.nulled_vcf
    Array[File] nulled_vcf_indicies = SetGenotypesNullTask.nulled_vcf_index
  }
}

task SetGenotypesNullTask {
  input {
    File vcf
    File samples_list
    Array[String]? contigs
    String runtime_docker

    Float? memory_gib
    Int? disk_gb
    Int? cpus
    Int? preemptible_tries
    Int? max_retries
    Int? boot_disk_gb
  }

  Float disk_size = size(vcf, "GB") * 2.2 + size(samples_list, "GB") + 16
  String output_vcf = "nulled-${basename(vcf)}"
  String output_vcf_index = "${output_vcf}.tbi"
  Array[String] contigs_set = select_first([contigs, []])

  runtime {
    memory: "${select_first([memory_gib, 2])} GiB"
    disks: "local-disk ${select_first([disk_gb, ceil(disk_size)])} HDD"
    cpus: select_first([cpus, 1])
    preemptible: select_first([preemptible_tries, 3])
    docker: runtime_docker
    bootDiskSizeGb: select_first([boot_disk_gb, 16])
  }

  command <<<
    bash /opt/task-scripts/SetGenotypesNull/SetGenotypesNullTask.bash '~{vcf}' \
      '~{output_vcf}' '~{samples_list}' '~{write_lines(contigs_set)}'
  >>>

  output {
    File nulled_vcf = output_vcf
    File nulled_vcf_index = output_vcf_index
  }
}
