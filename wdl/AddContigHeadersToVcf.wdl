version 1.0

# Add chr1-22,X,Y contig headers to a VCF
workflow AddContigHeadersToVcf {
  input {
    Array[File] vcfs
    File contig_fai
    String base_docker
  }

  scatter (vcf in vcfs) {
    call AddHeaders {
      input:
        vcf = vcf,
        contig_fai = contig_fai,
        base_docker = base_docker
    }
  }

  output {
    Array[File] output_vcfs = AddHeaders.headered_vcf
    Array[File] output_vcf_indices = AddHeaders.headered_vcf_index
  }
}

task AddHeaders {
  input {
    File vcf
    File contig_fai
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
    contig_fai='~{contig_fai}'
    output_vcf='~{output_vcf}'

    bcftools reheader --fai "${contig_fai}" --output "${output_vcf}" \
      "${vcf}"
    bcftools index --tbi "${output_vcf}"
  >>>

  output {
    File headered_vcf = "${output_vcf}"
    File headered_vcf_index = "${output_vcf_index}"
  }
}
