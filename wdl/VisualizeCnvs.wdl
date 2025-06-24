version 1.0

# Plots CNV depth profiles across batches
workflow VisualizeCnvs {
  input{
    Array[File] vcfs
    String plot_prefix
    Int min_size
    Int variants_per_shard = 40

    File sample_table # TSV with sample_id, sample_set_id
    Array[String] sample_set_ids
    Array[File] median_files
    Array[File] bincov_files
    Array[File] bincov_index_files

    String base_docker
    String r_docker
  }

  scatter (vcf in vcfs) {
    call ExtractVariants {
      input:
        vcf = vcf,
        min_size = min_size,
        variants_per_shard = variants_per_shard,
        base_docker = base_docker
    }
  }

  scatter (i in range(length(bincov_files))) {
    call SubsetBincovMatrix {
      input:
        bincov = bincov_files[i],
        bincov_index = bincov_index_files[i],
        intervals = ExtractVariants.intervals,
        base_docker = base_docker
    }
  }

  Array[File] flat_shards = flatten(ExtractVariants.variants)
  scatter (shard in flat_shards) {
    call MakePlots {
      input:
        variants = shard,
        sample_set_ids = sample_set_ids,
        bincov_files = SubsetBincovMatrix.bincov_subset,
        bincov_index_files= SubsetBincovMatrix.bincov_subset_index,
        median_files = median_files,
        sample_table = sample_table,
        r_docker = r_docker
    }
  }

  call MergePlotTars {
    input:
      plot_tars = MakePlots.plots,
      plot_prefix = plot_prefix,
      base_docker = base_docker
  }

  output {
    File cnv_plots = MergePlotTars.plots
  }
}

task ExtractVariants {
  input {
    File vcf
    Int min_size
    Int variants_per_shard
    String base_docker
  }

  Float disk_size = size(vcf, "GB") * 2 + 16

  runtime {
    bootDiskSizeGb: 8
    cpu: 1
    disks: "local-disk ${ceil(disk_size)} HDD"
    docker: base_docker
    maxRetries: 1
    memory: "4 GiB"
    preemptible: 3
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    vcf='~{vcf}'
    min_size='~{min_size}'
    variants_per_shard='~{variants_per_shard}'

    bcftools query --include "(SVTYPE == \"DEL\" || SVTYPE == \"DUP\") & FILTER ~ \"PASS\" & SVLEN >= ${min_size} & GT ~ \"1\"" \
      --format '%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVTYPE\t[%SAMPLE,]\n' "${vcf}" \
      | awk -F'\t' '{sub(/,$/, "", $6); print}' OFS='\t' > cnvs.tsv

    # merge intervals that are close to prevent duplicates when retrieving with
    # tabix
    awk -F'\t' '{print $1"\t"($2 - 1)"\t"$3}' cnvs.tsv \
      | LC_ALL=C sort -k1,1 -k2,2n \
      | bedtools merge -i stdin -d 101 > merged_intervals.bed

    split -l ${variants_per_shard} cnvs.tsv cnvs_
  >>>

  output {
    Array[File] variants = glob("cnvs_*")
    File intervals = "merged_intervals.bed"
  }
}

task SubsetBincovMatrix {
  input {
    File bincov
    File bincov_index
    Array[File] intervals
    String base_docker
  }

  Float disk_size = size(bincov, "GB") * 2 + 16

  runtime {
    bootDiskSizeGb: 8
    cpu: 1
    disks: "local-disk ${ceil(disk_size)} HDD"
    docker: base_docker
    maxRetries: 1
    memory: "4 GiB"
    preemptible: 3
  }

  String bincov_bn = basename(bincov)

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    intervals='~{write_lines(intervals)}'
    bincov='~{bincov}'
    bincov_bn='~{bincov_bn}'

    # need to expand the intervals so that padding can be added when plotting
    cat "${intervals}" \
      | xargs cat \
      | awk -F'\t' '{s=$3 - $2; pad=int(s / 2) + 1; min=$2 - pad; $2=min < 0 ? 0 : min; $3+=pad; print}' OFS='\t' \
      | sort -k1,1 -k2,2n \
      | bedtools merge -i stdin -d 101 \
      | awk -F'\t' '{$2+=1} 1' OFS="\t" > merged_intervals.tsv
    tabix --print-header --regions merged_intervals.tsv "${bincov}" \
      | bgzip > "${bincov_bn}"
    tabix --zero-based --begin 2 --comment '#' --end 3 --sequence 1 "${bincov_bn}"
  >>>

  output {
    File bincov_subset = bincov_bn
    File bincov_subset_index = bincov_bn + ".tbi"
  }
}

task MakePlots {
  input {
    File variants
    Array[String] sample_set_ids
    Array[File] bincov_files
    Array[File] bincov_index_files
    Array[File] median_files
    File sample_table
    String r_docker
  }

  Int variant_count = length(read_lines(variants))
  Float input_size = size(bincov_files, "GB")
    + size(bincov_index_files, "GB")
    + size(median_files, "GB")
    + size(sample_table, "GB")
    + size(variants, "GB")
  Float disk_size = input_size + variant_count * 0.01

  runtime {
    bootDiskSizeGb: 8
    cpu: 2
    disks: "local-disk ${ceil(disk_size)} HDD"
    docker: r_docker
    maxRetries: 1
    memory: "8 GiB"
    preemptible: 3
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    batch_ids='~{write_lines(sample_set_ids)}'
    bincov_files='~{write_lines(bincov_files)}'
    median_files='~{write_lines(median_files)}'
    variants='~{variants}'
    sample_table='~{sample_table}'

    paste "${batch_ids}" "${bincov_files}" > bincov_map.tsv
    paste "${batch_ids}" "${median_files}" > median_map.tsv

    Rscript /opt/gatk-sv-utils/scripts/visualize_cnvs.R \
      --cnvs "${variants}" \
      --sample-batches "${sample_table}" \
      --coverage-paths bincov_map.tsv \
      --medians-paths median_map.tsv \
      --output plots

    tar -czf plots.tar.gz plots
  >>>

  output {
    File plots = "plots.tar.gz"
  }
}

task MergePlotTars {
  input {
    Array[File] plot_tars
    String plot_prefix
    String base_docker
  }

  Float disk_size = size(plot_tars, "GB") * 3 + 16

  runtime {
    bootDiskSizeGb: 8
    cpu: 1
    disks: "local-disk ${ceil(disk_size)} HDD"
    docker: base_docker
    maxRetries: 1
    memory: "1 GiB"
    preemptible: 3
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    plot_tars='~{write_lines(plot_tars)}'
    plot_prefix='~{plot_prefix}'

    mkdir temp
    while read -r f; do
      tar -xzf "${f}"
      mv -t temp plots/*.jpg
      rm -r plots
    done < "${plot_tars}"

    mv temp "${plot_prefix}"
    tar -czf "${plot_prefix}.tar.gz" "${plot_prefix}"
  >>>

  output {
    File plots = "${plot_prefix}.tar.gz"
  }
}
