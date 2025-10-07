version 1.0

# Add INFO(s) to a VCF
workflow AddVcfInfos {
  input {
    Array[File] vcfs
    File new_infos
    File? new_headers
    String python_docker
  }

  scatter (vcf in vcfs) {
    call AddInfos {
      input:
        vcf = vcf,
        new_infos = new_infos,
        new_headers = new_headers,
        python_docker = python_docker
    }
  }

  output {
    Array[File] output_vcfs = AddInfos.modified_vcf
    Array[File] output_vcf_indices = AddInfos.modified_vcf_index
  }
}

task AddInfos {
  input {
    File vcf
    File new_infos
    File? new_headers
    String python_docker
  }

  Float disk_size = size(vcf, "GB") * 2 + 50

  runtime {
    bootDiskSizeGb: 8
    cpus: 1
    disks: "local-disk ${ceil(disk_size)} SSD"
    docker: python_docker
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
    new_infos='~{new_infos}'
    new_headers='~{if defined(new_headers) then new_headers else ""}'
    output_vcf='~{output_vcf}'

    if [[ -n "${new_headers}" ]]; then
      python /opt/gatk-sv-utils/scripts/add_vcf_infos.py \
        "${vcf}" "${output_vcf}" "${new_infos}" "${new_headers}"
    else
      python /opt/gatk-sv-utils/scripts/add_vcf_infos.py \
        "${vcf}" "${output_vcf}" "${new_infos}"
    fi
  >>>

  output {
    File modified_vcf = "${output_vcf}"
    File modified_vcf_index = "${output_vcf_index}"
  }
}
