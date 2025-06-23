version 1.0

workflow GatherGenotypes {
  input {
    File vcf
    File vcf_index
    String vcf_include_filter = 'FILTER == "PASS" && SVTYPE != "CNV" && SVTYPE != "BND" && SVTYPE != "CTX"'
    Int genotypes_per_shard = 100000000
    String runtime_docker
  }

  parameter_meta {
    vcf: "BCF file works too."
    vcf_index: "The index must have per contig stats, i.e. bcftools index --stats must work."
    vcf_include_filter: "Only used when gathering genotypes, not splitting VCF samples."
    genotype_per_shard: "Number of genotypes per shard. This approximate because divisibility and site filters."
  }

  call SplitVcfSamples {
    input:
      vcf = vcf,
      vcf_index = vcf_index,
      genotypes_per_shard = genotypes_per_shard,
      runtime_docker = runtime_docker
  }

  scatter (i in range(length(SplitVcfSamples.shards))) {
    call GatherGenotypesForSamples {
      input:
        vcf = vcf,
        samples = SplitVcfSamples.shards[i],
        output_prefix = i,
        vcf_include_filter = vcf_include_filter,
        runtime_docker = runtime_docker
    }
  }

  output {
    Array[File] genotypes = GatherGenotypesForSamples.genotypes
  }
}

task SplitVcfSamples {
  input {
    File vcf
    File vcf_index
    Int genotypes_per_shard
    String runtime_docker

    Float? memory_gib
    Int? boot_disk_gb
    Int? cpus
    Int? disk_gb
    Int? max_retries
    Int? preemptible_tries
  }

  Float disk_size = size(vcf, "GB") + 16

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
    set -o errexit
    set -o nounset
    set -o pipefail

    gawk -f /opt/gatk-sv-utils/scripts/split_vcf_samples.awk \
        '~{vcf}' '~{genotypes_per_shard}' shards
  >>>

  output {
    Array[File] shards = glob("shards/*.list")
  }
}

task GatherGenotypesForSamples {
  input {
    File vcf
    File samples
    Int output_prefix
    String vcf_include_filter
    String runtime_docker

    Float? memory_gib
    Int? boot_disk_gb
    Int? cpus
    Int? disk_gb
    Int? max_retries
    Int? preemptible_tries
  }

  Float disk_size = size(vcf, "GB") * 2 + 16
  String output_file = "${output_prefix}-genotypes.tsv.zst"

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

    vcf='~{vcf}'
    filter='~{vcf_include_filter}'
    samples='~{samples}'
    output='~{output_file}'

    bcftools query --include "${filter}" --samples-file "${samples}" \
      --format '[%SAMPLE\t%ID\t%INFO/SVTYPE\t%INFO/SVLEN\t%GT\t%GQ\n]' "${vcf}" \
      | gawk -i logging \
          'BEGIN{FS="\t"; OFS="\t"}
           $5 == "0/0" {$5 = 0; print; next}
           $5 == "1/0" || $5 == "0/1" {$5 = 1; print; next}
           $5 == "1/1" {$5 = 2; print; next}
           $5 == "./." {$5 = "."; print; next}
           {logging::log_err("unknown genotype in record: " $0); exit 1}' \
      | zstd -c > "${output}"
  >>>

  output {
    File genotypes = output_file
  }
}
