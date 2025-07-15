version 1.0

# Split a VCF by contig
workflow SplitVcfByContig {
  input {
    File vcf
    File vcf_index
    # One of `contig_list` or `contig_list_file` must be given. `contig_list` overrides
    # `contig_list_file`. There will be one VCF generated for each contig in the list.
    Array[String]? contig_list
    File? contig_list_file
    String runtime_docker
  }

  Array[String] contigs = select_first([contig_list, read_lines(select_first([contig_list_file]))])
  scatter (i in range(length(contigs))) {
    call GetContigFromVcf {
      input:
        vcf = vcf,
        vcf_index = vcf_index,
        contig = contigs[i],
        runtime_docker = runtime_docker
    }
  }

  output {
    Array[File] contig_vcfs = GetContigFromVcf.contig_vcf
    Array[File] contig_vcf_indicies = GetContigFromVcf.contig_vcf_index
  }
}

task GetContigFromVcf {
  input {
    File vcf
    File vcf_index
    String contig
    String runtime_docker

    Float? memory_gib
    Int? disk_gb
    Int? cpus
    Int? preemptible_tries
    Int? max_retries
    Int? boot_disk_gb
  }

  Float disk_size = size([vcf, vcf_index], "GB") * 1.2 + 16
  String output_vcf = "${contig}-${basename(vcf)}"
  String output_vcf_index = "${output_vcf}.tbi"

  runtime {
    bootDiskSizeGb: select_first([boot_disk_gb, 16])
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

    in_vcf='~{vcf}'
    contig='~{contig}'
    out_vcf='~{output_vcf}'

    { bcftools index --stats "${in_vcf}" | cut -f 1 | grep -F -w "${contig}" -; } \
      || { printf '%s not in vcf\n' "${contig}" >&2 && exit 1; }
    trap 'rm -f "${tmpfile}"' EXIT
    tmpfile="$(mktemp -p "${PWD}")"
    bcftools head "${in_vcf}" \
      | gawk '/^##contig=/{match($0, /<ID=([^,>]+)/, a); if (RSTART && t == a[1]){print} next} 1' t="${contig}" - > "${tmpfile}"
    bcftools view --output-type u --regions "${contig}" "${in_vcf}" \
      | bcftools reheader --header "${tmpfile}" - \
      | bcftools view --output-type z --output "${out_vcf}" -
    bcftools index --tbi "${out_vcf}"
  >>>

  output {
    File contig_vcf = "${output_vcf}"
    File contig_vcf_index = "${output_vcf_index}"
  }
}
