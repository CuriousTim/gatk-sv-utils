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

    String runtime_docker
  }

  scatter (vcf in vcfs) {
    call ExtractVariants {
      input:
        vcf = vcf,
        min_size = min_size,
        variants_per_shard = variants_per_shard,
        runtime_docker = runtime_docker
    }
  }

  scatter (i in range(length(bincov_files))) {
    call SubsetBincovMatrix {
      input:
        bincov = bincov_files[i],
        bincov_index = bincov_index_files[i],
        intervals = ExtractVariants.intervals,
        runtime_docker = runtime_docker
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
        runtime_docker = runtime_docker,
    }
  }

  call MergePlotTars {
    input:
      plot_tars = MakePlots.plots,
      plot_prefix = plot_prefix,
      runtime_docker = runtime_docker
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
    String runtime_docker
  }

  Float disk_size = size(vcf, "GB") * 2 + 16

  runtime {
    bootDiskSizeGb: 8
    cpu: 1
    disks: "local-disk ${ceil(disk_size)} HDD"
    docker: runtime_docker
    maxRetries: 1
    memory: "4 GiB"
    preemptible: 3
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    /opt/task_scripts/VisualizeCnvs/ExtractVariants '~{vcf}' '~{min_size}' \
      'cnvs_' '~{variants_per_shard}'  merged_intervals.bed
  >>>

  output {
    Array[File] variants = glob("cnvs_")
    File intervals = "merged_intervals.bed"
  }
}

task SubsetBincovMatrix {
  input {
    File bincov
    File bincov_index
    Array[File] intervals
    String runtime_docker
  }

  Float disk_size = size(bincov, "GB") * 2 + 16

  runtime {
    bootDiskSizeGb: 8
    cpu: 1
    disks: "local-disk ${ceil(disk_size)} HDD"
    docker: runtime_docker
    maxRetries: 1
    memory: "4 GiB"
    preemptible: 3
  }

  String bincov_bn = basename(bincov)

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    /opt/task_scripts/VisualizeCnvs/SubsetBincovMatrix '~{write_lines(intervals)}' \
      '~{bincov}' '~{bincov_bn}'
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
    String runtime_docker
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
    docker: runtime_docker
    maxRetries: 1
    memory: "8 GiB"
    preemptible: 3
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    /opt/task_scripts/VisualizeCnvs/MakePlots \
      '~{variants}' \
      '~{write_lines(sample_set_ids)}' \
      '~{write_lines(bincov_files)}' \
      '~{write_lines(median_files)}' \
      '~{sample_table}' \
      plots
  >>>

  output {
    File plots = "plots.tar.gz"
  }
}

task MergePlotTars {
  input {
    Array[File] plot_tars
    String plot_prefix
    String runtime_docker
  }

  Float disk_size = size(plot_tars, "GB") * 3 + 16

  runtime {
    bootDiskSizeGb: 8
    cpu: 1
    disks: "local-disk ${ceil(disk_size)} HDD"
    docker: runtime_docker
    maxRetries: 1
    memory: "1 GiB"
    preemptible: 3
  }

  command <<<
    /opt/task_scripts/VisualizeCnvs/MergePlotTars '~{write_lines(plot_tars)}' \
      '~{plot_prefix}'
  >>>

  output {
    File plots = "${plot_prefix}.tar.gz"
  }
}
