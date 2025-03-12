version 1.0

# Move BOTHSIDES_SUPPORT, HIGH_SR_BACKGROUND, and PESR_GT_OVERDISPERSION values
# from the FILTER field to the INFO field.
workflow MoveVcfFiltersToInfo {
  input {
    Array[File] vcfs
    # Either the same length as `vcfs` or length 1.
    # This is mostly to accomodate the case where the VCFs are split / named by
    # contig and the output should be named by contig as well.
    # One of `output_prefix_list` or `output_prefix_file` must be given.
    # `output_prefix_list` overrides `output_prefix_file`.
    Array[String]? output_prefix_list
    File? output_prefix_file
    String runtime_docker
  }

  Array[String] output_prefix = select_first([output_prefix_list, read_lines(select_first([output_prefix_file]))])

  scatter (i in range(length(vcfs))) {
    call MoveVcfFilters {
      input:
        vcf = vcfs[i],
        output_prefix = if (length(output_prefix) == 1) then output_prefix[0] else output_prefix[i],
        runtime_docker = runtime_docker
    }
  }

  output {
    Array[File] filters_moved_vcfs = MoveVcfFilters.filters_moved_vcf
    Array[File] filters_moved_vcf_indicies = MoveVcfFilters.filters_moved_vcf_index
  }
}

task MoveVcfFilters {
  input {
    File vcf
    String output_prefix
    String runtime_docker

    Float? memory_gib
    Int? disk_gb
    Int? cpus
    Int? preemptible_tries
    Int? max_retries
    Int? boot_disk_gb
  }

  Float disk_size = size(vcf, "GB") * 2.2 + 16
  String output_vcf = "${output_prefix}-filters_moved.vcf.gz"
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
    bash /opt/task-scripts/MoveVcfFiltersToInfo/MoveVcfFilters.bash '~{vcf}' '~{output_vcf}'
  >>>

  output {
    File filters_moved_vcf = output_vcf
    File filters_moved_vcf_index = output_vcf_index
  }
}
