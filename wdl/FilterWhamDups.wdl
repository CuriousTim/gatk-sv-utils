version 1.0

# Filter Wham-only DUPs
workflow FilterWhamDups {
  input {
    Array[File] vcfs
    File blacklist_regions
    # Wham-only DUPs matching `extra_filters` will be checked for overlap
    # against `blacklist_regions`. The ones overlapping the regions will
    # be removed from the VCF.
    String? extra_filters
    String runtime_docker
  }

  scatter (vcf in vcfs) {
    call FilterSites {
      input:
        vcf = vcf,
        blacklist = blacklist_regions,
        extra_filters = extra_filters,
        runtime_docker = runtime_docker
    }
  }

  output {
    Array[File] filtered_vcfs = FilterSites.filtered_vcf
    Array[File] filtered_vcf_indicies = FilterSites.filtered_vcf_index
  }
}

task FilterSites {
  input {
    File vcf
    File blacklist
    String extra_filters = '(EVIDENCE == "SR" || EVIDENCE == "RD,SR")'
    String runtime_docker

    Float? memory_gib
    Int? boot_disk_gb
    Int? cpus
    Int? disk_gb
    Int? max_retries
    Int? preemptible_tries
  }

  Float disk_size = size(vcf, "GB") * 2.2  + size(blacklist, "GB") * 2 + 16
  String output_vcf = "filtered-${basename(vcf)}"
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
    /opt/task_scripts/FilterWhamDups/FilterSites '~{vcf}' \
      '~{blacklist}' '~{output_vcf}' \
      '~{if defined(extra_filters) then extra_filters else ""}'
  >>>

  output {
    File filtered_vcf = "${output_vcf}"
    File filtered_vcf_index = "${output_vcf_index}"
  }
}
