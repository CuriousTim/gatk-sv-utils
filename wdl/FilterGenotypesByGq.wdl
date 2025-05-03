version 1.0

# Set genotypes to missing where the GQ is lower than the threshold. The
# thresholds are specific to an SV type.
workflow FilterGenotypesByGq {
  input {
    Array[File] vcfs
    Array[String] svtypes
    Array[Int] min_gqs
    Array[String]? output_prefix_list
    File? output_prefix_file

    String runtime_docker
  }

  Array[String] output_prefix = select_first([output_prefix_list, read_lines(select_first([output_prefix_file]))])

  scatter (i in range(length(vcfs))) {
    call UpdateGenotypes {
      input:
        vcf = vcfs[i],
        svtypes = svtypes,
        min_gqs = min_gqs,
        output_prefix = output_prefix[i],
        runtime_docker = runtime_docker
    }
  }

  output {
    Array[File] filtered_vcfs = UpdateGenotypes.updated_vcf
    Array[File] filtered_vcf_indicies = UpdateGenotypes.updated_vcf_index
  }
}

task UpdateGenotypes {
  input {
    File vcf
    Array[String] svtypes
    Array[Int] min_gqs
    String output_prefix
    String runtime_docker

    Float? memory_gib
    Int? disk_gb
    Int? cpus
    Int? preemptible_tries
    Int? max_retries
    Int? boot_disk_gb
  }

  Float disk_size = size(vcf, "GB") * 2 + 16
  String output_vcf = "${output_prefix}-gq_filtered.vcf.gz"
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
    /opt/task_scripts/FilterGenotypesByGq/UpdateGenotypes '~{vcf}' \
      '~{output_vcf}' '~{write_lines(svtypes)}' '~{write_lines(min_gqs)}'
  >>>

  output {
    File updated_vcf = output_vcf
    File updated_vcf_index = output_vcf_index
  }
}
