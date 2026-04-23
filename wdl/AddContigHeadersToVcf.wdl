version 1.0

# Add chr1-22,X,Y contig headers to a VCF
workflow AddContigHeadersToVcf {
  input {
    Array[File] vcfs
    # index file (.fai) for an hg38 FASTA file
    File reference_index
    String base_docker
  }

  scatter (vcf in vcfs) {
    call AddHeaders {
      input:
        vcf = vcf,
        reference_index = reference_index,
        base_docker = base_docker
    }
  }

  output {
    Array[File] headered_vcfs = AddHeaders.headered_vcf
    Array[File] headered_vcf_idxs = AddHeaders.headered_vcf_index
  }
}

task AddHeaders {
  input {
    File vcf
    File reference_index
    String base_docker
  }

  Int disk_size = ceil(size(vcf, "GB") * 2) + 16

  runtime {
    bootDiskSizeGb: 8
    cpus: 1
    disks: "local-disk ${disk_size} HDD"
    docker: base_docker
    maxRetries: 1
    memory: "4 GiB"
    preemptible: 3
  }

  String headered_vcf_name = basename(vcf)

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    vcf='~{vcf}'
    reference_index='~{reference_index}'
    headered_vcf_name='~{headered_vcf_name}'

    bcftools reheader --fai "${reference_index}" --output "${headered_vcf_name}" "${vcf}"
    bcftools index --tbi "${headered_vcf_name}"
  >>>

  output {
    File headered_vcf = "${headered_vcf_name}"
    File headered_vcf_index = "${headered_vcf_name}.tbi"
  }
}
