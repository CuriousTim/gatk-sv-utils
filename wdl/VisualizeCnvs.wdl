version 1.0

# Plots CNV depth profiles across batches
workflow VisualizeCnvs {
  input {
    Array[File] vcfs
    String plot_prefix
    Int min_size
    Int variants_per_shard = 40

    File sample_table # TSV with sample_id, sample_set_id
    Array[String] sample_set_ids
    Array[File] medians_files
    Array[File] bincov_files
    Array[File] bincov_index_files

    String base_docker
    String r_docker
  }

  scatter (vcf in vcfs) {
    call ExtractVariants {
      input:
        vcf = vcf,
        sample_table = sample_table,
        min_size = min_size,
        variants_per_shard = variants_per_shard,
        base_docker = base_docker
    }
  }

  call MergeIntervals {
    input:
      intervals = flatten(ExtractVariants.intervals),
      base_docker = base_docker
  }

  call MapFromArrays as MakeBincovMap {
    input:
      keys = sample_set_ids,
      values = bincov_files,
      base_docker = base_docker
  }

  call MapFromArrays as MakeBincovIndexMap {
    input:
      keys = sample_set_ids,
      values = bincov_index_files,
      base_docker = base_docker
  }

  scatter (f in MergeIntervals.merged_intervals) {
    String batch = basename(f)
    call SubsetBincovMatrix {
      input:
        intervals = f,
        bincov = MakeBincovMap.m[batch],
        bincov_index = MakeBincovIndexMap.m[batch],
        base_docker = base_docker
    }
  }

  call MapFromArrays as MakeBincovSubsetMap {
    input:
      keys = SubsetBincovMatrix.batch,
      values = SubsetBincovMatrix.bincov_subset,
      base_docker = base_docker
  }

  call MapFromArrays as MakeBincovSubsetIndexMap {
    input:
      keys = SubsetBincovMatrix.batch,
      values = SubsetBincovMatrix.bincov_subset_index,
      base_docker = base_docker
  }

  call MapFromArrays as MakeMediansMap {
    input:
      keys = sample_set_ids,
      values = medians_files,
      base_docker = base_docker
  }

  call ShardVariants {
    input:
      variants = ExtractVariants.variants,
      sample_table = sample_table,
      variants_per_shard = variants_per_shard,
      bincov_subset_map = write_map(MakeBincovSubsetMap.m),
      bincov_subset_index_map = write_map(MakeBincovSubsetIndexMap.m),
      medians_map = write_map(MakeMediansMap.m),
      base_docker = base_docker
  }

  scatter (i in range(length(ShardVariants.shards))) {
    call MakePlots {
      input:
        variants = ShardVariants.shards[i],
        sample_set_ids = sample_set_ids,
        bincov_files = read_lines(ShardVariants.bincov_paths[i]),
        bincov_index_files= read_lines(ShardVariants.bincov_index_paths[i]),
        medians_files = read_lines(ShardVariants.medians_paths[i]),
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
    File sample_table
    String base_docker
  }

  Float disk_size = size([vcf, sample_table], "GB") * 3 + 16

  runtime {
    bootDiskSizeGb: 8
    cpu: 2
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
    sample_table='~{sample_table}'

    bcftools query --include "(SVTYPE == \"DEL\" || SVTYPE == \"DUP\") & FILTER ~ \"PASS\" & SVLEN >= ${min_size} & GT ~ \"1\"" \
      --format '%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVTYPE\t[%SAMPLE,]\n' "${vcf}" \
      | awk -F'\t' '{sub(/,$/, "", $6); print}' OFS='\t' > cnvs.tsv

    mkdir intervals
    awk -F'\t' 'BEGIN{OFS="\t"}
    ARGIND == 1{a[$1]=$2}
    ARGIND == 2{
      split($6, b, /,/);
      for (i in b) {
        print $1,$2 - 1,$3 > ("intervals/" a[b[i]])
      }
    }' "${sample_table}" cnvs.tsv
  >>>

  output {
    File variants = "cnvs.tsv"
    Array[File] intervals = glob("intervals/*")
  }
}

task MergeIntervals {
  input {
    Array[File] intervals
    String base_docker
  }

  Float disk_size = size(intervals, "GB") * 3 + 16

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

    intervals='~{write_lines(intervals)}'

    mkdir intervals
    while read -r f; do
      batch="$(basename "${f}")"
      cat "${f}" >> "intervals/${batch}"
    done < "${intervals}"

    mkdir merged_intervals
    for f in intervals/*; do
      batch="$(basename "${f}")"
      awk -F'\t' '{s=$3 - $2; pad=int(s / 2) + 1; min=$2 - pad; $2=min < 0 ? 0 : min; $3+=pad} 1' OFS='\t' "${f}" \
        | sort -k1,1 -k2,2n \
        | bedtools merge -i stdin -d 101 \
        | awk -F'\t' '{$2+=1} 1' OFS='\t' > "merged_intervals/${batch}"
    done
  >>>

  output {
    Array[File] merged_intervals = glob("merged_intervals/*")
  }
}

task SubsetBincovMatrix {
  input {
    File bincov
    File bincov_index
    File intervals
    String base_docker
  }

  Float disk_size = size(bincov, "GB") * 2 + 16

  runtime {
    bootDiskSizeGb: 8
    cpu: 2
    disks: "local-disk ${ceil(disk_size)} HDD"
    docker: base_docker
    maxRetries: 1
    memory: "4 GiB"
    preemptible: 2
  }

  String intervals_bn = basename(intervals)
  String output_name = intervals_bn + ".gz"

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    intervals='~{intervals}'
    bincov='~{bincov}'
    bincov_subset='~{output_name}'

    tabix --print-header --regions "${intervals}" "${bincov}" \
      | bgzip > "${bincov_subset}"
    tabix --zero-based --begin 2 --comment '#' --end 3 --sequence 1 "${bincov_subset}"
  >>>

  output {
    String batch = intervals_bn
    File bincov_subset = output_name
    File bincov_subset_index = output_name + ".tbi"
  }
}

task ShardVariants {
  input {
    Array[File] variants
    File sample_table
    Int variants_per_shard
    File bincov_subset_map
    File bincov_subset_index_map
    File medians_map
    String base_docker
  }

  Float disk_size = size(variants, "GB") * 3
    + size(sample_table, "GB")
    + size(bincov_subset_map, "GB")
    + size(bincov_subset_index_map, "GB")
    + size(medians_map, "GB")

  runtime {
    bootDiskSizeGb: 8
    cpu: 1
    disks: "local-disk ${ceil(disk_size)} HDD"
    docker: base_docker
    maxRetries: 1
    memory: "2 GiB"
    preemptible: 3
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    variants='~{write_lines(variants)}'
    sample_table='~{sample_table}'
    variants_per_shard=~{variants_per_shard}
    bincov_subset_map='~{bincov_subset_map}'
    bincov_subset_index_map='~{bincov_subset_index_map}'
    medians_map='~{medians_map}'

    cat "${variants}" | xargs cat > merged_cnvs.tsv
    mkdir shards bincov bincov_index medians
    split -l "${variants_per_shard}" merged_cnvs.tsv shards/cnvs_

    gawk -F'\t' '
    ARGIND == 1 {
      sample_arr[$1] = $2
    }
    ARGIND == 2 {
      bincov_arr[$1] = $2
    }
    ARGIND == 3 {
      bincov_index_arr[$1] = $2
    }
    ARGIND == 4 {
      medians_arr[$1] = $2
    }
    ARGIND > 4 && FNR == 1 {
      if (bincov_out) {
        close(bincov_out)
      }
      if (bincov_index_out) {
        close(bincov_index_out)
      }
      if (medians_out) {
        close(medians_out)
      }
      n = split(FILENAME, p, /\//)
      bn = p[n]
      bincov_out = "sort -u > bincov/" bn
      bincov_index_out = "sort -u > bincov_index/" bn
      medians_out = "sort -u > medians/" bn
    }
    ARGIND > 4 {
      split($6, b, /,/);
      for (i in b) {
        batch = sample_arr[b[i]]
        print bincov_arr[batch] | bincov_out
        print bincov_index_arr[batch] | bincov_index_out
        print medians_arr[batch] | medians_out
      }
    }' "${sample_table}" "${bincov_subset_map}" "${bincov_subset_index_map}" \
      "${medians_map}" shards/cnvs_*
  >>>

  output {
    Array[File] shards = glob("shards/cnvs_*")
    Array[File] bincov_paths = glob("bincov/cnvs_*")
    Array[File] bincov_index_paths = glob("bincov_index/cnvs_*")
    Array[File] medians_paths = glob("medians/cnvs_*")
  }
}

task MakePlots {
  input {
    File variants
    Array[String] sample_set_ids
    Array[File] bincov_files
    Array[File] bincov_index_files
    Array[File] medians_files
    File sample_table
    String r_docker
  }

  Int variant_count = length(read_lines(variants))
  Float input_size = size(bincov_files, "GB")
    + size(bincov_index_files, "GB")
    + size(medians_files, "GB")
    + size(sample_table, "GB")
    + size(variants, "GB")
  Float disk_size = input_size + variant_count * 0.01 + 16

  runtime {
    bootDiskSizeGb: 8
    cpu: 4
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
    medians_files='~{write_lines(medians_files)}'
    variants='~{variants}'
    sample_table='~{sample_table}'

    paste "${batch_ids}" "${bincov_files}" > bincov_map.tsv
    paste "${batch_ids}" "${medians_files}" > medians_map.tsv

    Rscript /opt/gatk-sv-utils/scripts/visualize_cnvs.R \
      --cnvs "${variants}" \
      --sample-batches "${sample_table}" \
      --coverage-paths bincov_map.tsv \
      --medians-paths medians_map.tsv \
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

task MapFromArrays {
  input {
    Array[String] keys
    Array[String] values
    String base_docker
  }

  runtime {
    bootDiskSizeGb: 8
    cpu: 1
    disks: "local-disk 16 HDD"
    docker: base_docker
    maxRetries: 1
    memory: "1 GiB"
    preemptible: 3
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    paste '~{write_lines(keys)}' '~{write_lines(values)}' > map.tsv
  >>>

  output {
    Map[String, String] m = read_map("map.tsv")
  }
}
