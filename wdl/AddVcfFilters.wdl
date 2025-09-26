version 1.0

# Add FILTER(s) to a VCF
workflow AddVcfFilters {
  input {
    Array[File] vcfs
    File new_filters
    File? new_headers
    String base_docker
  }

  scatter (vcf in vcfs) {
    call AddFilters {
      input:
        vcf = vcf,
        new_filters = new_filters,
        new_headers = new_headers,
        base_docker = base_docker
    }
  }

  output {
    Array[File] output_vcfs = AddFilters.modified_vcf
    Array[File] output_vcf_indices = AddFilters.modified_vcf_index
  }
}

task AddFilters {
  input {
    File vcf
    File new_filters
    File? new_headers
    String base_docker
  }

  Float disk_size = size(vcf, "GB") * 2 + 16

  runtime {
    bootDiskSizeGb: 8
    cpus: 1
    disks: "local-disk ${ceil(disk_size)}  HDD"
    docker: base_docker
    maxRetries: 1
    memory: "4 GiB"
    preemptible: 3
  }

  String output_vcf = basename(vcf)
  String output_vcf_index = basename(vcf) + ".tbi"

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    vcf='~{vcf}'
    new_filters='~{new_filters}'
    new_headers='~{if defined(new_headers) then new_headers else ""}'
    output_vcf='~{modified_vcf}'

    if [[ -n "${new_headers}" ]]; then
      add_vcf_filters "${vcf}" "${output_vcf}" "${new_filters}" "${new_headers}"
    else
      add_vcf_filters "${vcf}" "${output_vcf}" "${new_filters}"
    fi

    bcftools index --tbi "${output_vcf}"
  >>>

  output {
    File modified_vcf = "${output_vcf}"
    File modified_vcf_index = "${output_vcf_index}"
  }
}
