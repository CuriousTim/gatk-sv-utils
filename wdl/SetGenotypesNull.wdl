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

    in_vcf='~{vcf}'
    out_vcf='~{output_vcf}'
    samples='~{samples_list}'
    contigs='~{write_lines(contigs_set)}'

    wc -l "${samples}" | read -r nsamples other
    if (( nsamples == 0 )); then
      printf 'expected at least one sample. none given.\n' > /dev/stderr
      exit 1
    fi

    bgzip -cd "${in_vcf}" \
    | gawk -f /opt/gatk-sv-utils/scripts/set_gt_null.awk "${contigs}" "${samples}" - \
    | bgzip -c > "${out_vcf}"
    bcftools index --tbi "${out_vcf}"
  >>>

  output {
    File nulled_vcf = output_vcf
    File nulled_vcf_index = output_vcf_index
  }
}
