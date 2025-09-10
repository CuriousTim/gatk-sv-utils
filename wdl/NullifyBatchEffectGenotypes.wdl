version 1.0

# Nullify genotypes driven by batch effects
workflow NullifyBatchEffectGenotypes {
  input {
    Array[File] vcfs
    File batch_table
    File variant_table
    String vcf_prefix
    String python_docker
  }

  scatter (vcf in vcfs) {
    call NullifyGenotypes {
      input:
        vcf = vcf,
        batch_table = batch_table,
        variant_table = variant_table,
        vcf_prefix = vcf_prefix,
        python_docker = python_docker
    }
  }

  output {
    Array[File] nulled_vcfs = NullifyGenotypes.nulled_vcf
    Array[File] nulled_vcf_indicies = NullifyGenotypes.nulled_vcf_index
  }
}

task NullifyGenotypes {
  input {
    File vcf
    File batch_table
    File variant_table
    String vcf_prefix
    String python_docker
  }

  Float disk_size = size(vcf, "GB") * 2 + 16

  runtime {
    bootDiskSizeGb: 8
    cpus: 1
    disks: "local-disk ${ceil(disk_size)} HDD"
    docker: python_docker
    maxRetries: 1
    memory: "4 GiB"
    preemptible: 3
  }

  String nulled_vcf_name = vcf_prefix + basename(vcf)
  String nulled_vcf_index_name = nulled_vcf_name + ".tbi"

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    vcf='~{vcf}'
    batch_table='~{batch_table}'
    variant_table='~{variant_table}'
    nulled_vcf_name='~{nulled_vcf_name}'

    python3 /opt/gatk-sv-utils/scripts/nullify_batch_effect_gts.py \
      "${vcf}" "${nulled_vcf_name}" "${batch_table}" "${variant_table}"
  >>>

  output {
    File nulled_vcf = nulled_vcf_name
    File nulled_vcf_index = nulled_vcf_index_name
  }
}
