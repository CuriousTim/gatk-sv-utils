version 1.0

# Pull all variants from cohort and batch VCFs that overlap a set of SVs
workflow ExtractSVsEvidenceFromVcf {
  input {
    File variants
    Array[String] sample_set_ids
    Array[File] cohort_vcfs
    Array[File] cohort_vcf_indexes
    Array[File] clustered_wham_vcfs
    Array[File] clustered_wham_vcf_indexes
    Array[File] clustered_manta_vcfs
    Array[File] clustered_manta_vcf_indexes
    Array[File] clustered_melt_vcfs
    Array[File] clustered_melt_vcf_indexes
    Array[File] clustered_depth_vcfs
    Array[File] clustered_depth_vcf_indexes
    String base_docker
  }

  scatter (i in range(length(cohort_vcfs))) {
    call SubsetVcf as subset_cohort {
      input:
        variants = variants,
        vcf = cohort_vcfs[i],
        vcf_index = cohort_vcf_indexes[i],
        vcf_label = "cohort",
        base_docker = base_docker
    }
  }

  scatter (i in range(length(sample_set_ids))) {
    call SubsetVcf as subset_manta {
      input:
        variants = variants,
        vcf = clustered_manta_vcfs[i],
        vcf_index = clustered_manta_vcf_indexes[i],
        vcf_label = "Manta",
        base_docker = base_docker
    }
  }

  scatter (i in range(length(sample_set_ids))) {
    call SubsetVcf as subset_melt {
      input:
        variants = variants,
        vcf = clustered_melt_vcfs[i],
        vcf_index = clustered_melt_vcf_indexes[i],
        vcf_label = "MELT",
        base_docker = base_docker
    }
  }

  scatter (i in range(length(sample_set_ids))) {
    call SubsetVcf as subset_wham {
      input:
        variants = variants,
        vcf = clustered_wham_vcfs[i],
        vcf_index = clustered_wham_vcf_indexes[i],
        vcf_label = "Wham",
        base_docker = base_docker
    }
  }

  scatter (i in range(length(sample_set_ids))) {
    call SubsetVcf as subset_depth {
      input:
        variants = variants,
        vcf = clustered_depth_vcfs[i],
        vcf_index = clustered_depth_vcf_indexes[i],
        vcf_label = "Depth",
        base_docker = base_docker
    }
  }

  call MergeVariants {
    input:
      cohort = subset_cohort.variants,
      wham = subset_wham.variants,
      manta = subset_manta.variants,
      melt = subset_melt.variants,
      depth = subset_depth.variants,
      base_docker = base_docker
  }

  output {
    File merged_variants = MergeVariants.merged_variants
  }
}

task SubsetVcf {
  input {
    # chr,start,end,vid,svtype,sample
    File variants
    File vcf
    File vcf_index
    String vcf_label
    String base_docker
  }

  Float disk_size = size(vcf, "GB") * 4 + 16

  runtime {
    bootDiskSizeGb: 8
    cpus: 1
    disks: "local-disk ${ceil(disk_size)}  HDD"
    docker: base_docker
    maxRetries: 1
    memory: "8 GiB"
    preemptible: 3
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    variants='~{variants}'
    vcf='~{vcf}'
    vcf_label='~{vcf_label}'

    cut -f 6 "${variants}" | LC_ALL=C sort -u > samples
    bcftools query --list-samples "${vcf}" | LC_ALL=C sort -u > vcf_samples
    LC_ALL=C comm -12 samples vcf_samples > common_samples
    if [[ ! -s common_samples ]]; then
      exit 0
    fi
    cut -f 1,2,3 "${variants}" | LC_ALL=C sort -k1,1 -k2,2n > coordinates

    bcftools query --samples-file common_samples --include 'GT="alt"' \
      --regions-file coordinates --format '[%CHROM\t%POS0\t%INFO/END\t%ID\t%INFO/SVTYPE\tSAMPLE]\n' \
      | LC_ALL=C sort -k1,1 -k2,2n \
      | gzip -c > database.bed.gz
    bedtools intersect -a "${variants}" -b "${database}" -wo -sorted \
      | gawk -F'\t' '
          BEGIN {OFS = "\t"; OFMT = "%.0f"}
          $5 != $11 && $6 != $12 { next }
          {
            print $1,$2 + 1,$3,$4,$5,$6,$7,$8 + 1,$9,$10,$11,vcf_label
          }' vcf_label="${vcf_label}" - \
      | gzip -c > variants.tsv.gz
  >>>

  output {
    File? variants = "variants.tsv.gz"
  }
}

task MergeVariants {
  input {
    Array[File?] cohort
    Array[File?] wham
    Array[File?] manta
    Array[File?] melt
    Array[File?] depth
    String base_docker
  }

  Array[File] cohort_variants = select_all(cohort)
  Array[File] wham_variants = select_all(wham)
  Array[File] manta_variants = select_all(manta)
  Array[File] melt_variants = select_all(melt)
  Array[File] depth_variants = select_all(depth)

  Float variants_size = size(cohort_variants, "GB") +
    size(wham_variants, "GB") +
    size(manta_variants, "GB") +
    size(melt_variants, "GB") +
    size(depth_variants, "GB")
  Float disk_size = variants_size * 2 + 16

  runtime {
    bootDiskSizeGb: 8
    cpus: 1
    disks: "local-disk ${ceil(disk_size)}  HDD"
    docker: base_docker
    maxRetries: 1
    memory: "4 GiB"
    preemptible: 3
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    cohort_variants='~{write_lines(cohort_variants)}'
    wham_variants='~{write_lines(wham_variants)}'
    manta_variants='~{write_lines(manta_variants)}'
    melt_variants='~{write_lines(melt_variants)}'
    depth_variants='~{write_lines(depth_variants)}'

    cat "${cohort}" "${wham_variants}" "${manta_variants}" "${melt_variants}" "${depth_variants}" \
      | xargs cat > merged_variants.tsv.gz
  >>>

  output {
    File merged_variants = "merged_variants"
  }
}
