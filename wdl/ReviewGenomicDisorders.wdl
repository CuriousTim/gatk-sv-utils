version 1.0

# Review potential genomic disorder CNVs
workflow ReviewGenomicDisorders {
  input {
    String tar_prefix
    File gd_regions
    Float? padding
    Float? min_rd_deviation
    Int? max_calls_per_sample
    File sample_table
    File segdups
    File pedigree
    Array[String] sample_set_ids
    Array[File] bincov_matrices
    Array[File] bincov_matrix_indices
    Array[File] medians_files

    String base_docker
    String r_docker
  }

  output {
    File gd_plots = MergePlotTars.gd_plots
    File samples_with_excess_calls = MergePlotTars.excess_call_samples
  }

  scatter (i in range(length(bincov_matrices))) {
    call VisualizeGenomicDisorders {
      input:
        batch_id = sample_set_ids[i],
        bincov = bincov_matrices[i],
        bincov_index = bincov_matrix_indices[i],
        medians_file = medians_files[i],
        gd_regions = gd_regions,
        sample_table = sample_table,
        segdups = segdups,
        pedigree = pedigree,
        min_rd_deviation = min_rd_deviation,
        padding = padding,
        max_calls_per_sample = max_calls_per_sample,
        r_docker = r_docker
    }
  }

  call MergePlotTars {
    input:
      plot_tars = VisualizeGenomicDisorders.plots,
      tar_prefix = tar_prefix,
      excess_call_samples = VisualizeGenomicDisorders.excess_call_samples,
      base_docker = base_docker
  }
}

task VisualizeGenomicDisorders {
  input {
    String batch_id
    File sample_table
    File bincov
    File bincov_index
    File medians_file
    File gd_regions
    File pedigree
    File segdups
    Float? min_rd_deviation
    Float? padding
    Float? max_calls_per_sample
    String r_docker
  }

  parameter_meta {
    batch_id: "Batch ID."
    sample_table: "Two-column table of batches and samples."
    bincov: "Binned coverage matrix file."
    bincov_index: "Binned coverage matrix index file."
    medians_file: "Genome-wide median coverage file."
    gd_regions: "Genomic disorder regions coordinates file."
    pedigree: "Pedigree."
    segdups: "Segmental duplication coordinates file."
    min_rd_deviation: "Minimum read depth ratio deviation from 1 required to make a plot."
    padding: "Fraction of GD region to add as padding."
    max_calls_per_sample: "Maximum number of call-regions per sample."
    r_docker: "Docker image."
  }

  output {
    File plots = "${batch_id}.tar"
    File excess_call_samples = "excess_calls_samples_${batch_id}.tsv"
  }

  Float input_size = size([bincov, medians_file, gd_regions, segdups, pedigree], "GB")
  Int disk_size = ceil(input_size) + 50
  runtime {
    bootDiskSizeGb: 8
    cpu: 1
    disks: "local-disk ${disk_size} SSD"
    docker: r_docker
    maxRetries: 1
    memory: "16 GiB"
    preemptible: 3
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    batch_id='~{batch_id}'
    sample_table='~{sample_table}'
    bincov='~{bincov}'
    medians_file='~{medians_file}'
    gd_regions='~{gd_regions}'
    pedigree='~{pedigree}'
    segdups='~{segdups}'
    min_rd_deviation='~{if defined(min_rd_deviation) then "--min-shift ${min_rd_deviation}" else ""}'
    padding='~{if defined(padding) then "--pad ${padding}" else ""}'
    max_calls_per_sample='~{if defined(max_calls_per_sample) then "--max-calls-per-sample ${max_calls_per_sample}" else ""}'

    awk '$1 == bid {print $2}' bid="${batch_id}" "${sample_table}" > samples.list
    cut -f2,5 "${pedigree}" > ploidy.tsv

    Rscript /opt/gatk-sv-utils/scripts/visualize_gd.R \
      ${min_rd_deviation} \
      ${padding} \
      ${max_calls_per_sample} \
      --violators "excess_calls_samples_${batch_id}.tsv" \
      "${gd_regions}" \
      "${segdups}" \
      "${bincov}" \
      "${medians_file}" \
      samples.list \
      ploidy.tsv \
      plots

    touch "excess_calls_samples_${batch_id}.tsv"
    tar -cf "${batch_id}.tar" plots
  >>>
}

task MergePlotTars {
  input {
    Array[File] plot_tars
    Array[File] excess_call_samples
    String tar_prefix
    String base_docker
  }

  parameter_meta {
    plot_tars: "Archives of plots."
    excess_call_samples: "Samples with excess number of calls and their counts."
    tar_prefix: "What to name the output archive."
    base_docker: "Docker image."
  }

  output {
    File gd_plots = "${tar_prefix}.tar"
    File excess_call_samples = "${tar_prefix}-excess_gds.tsv"
  }

  Float input_size = size(plot_tars, "GB")
  Int disk_size = ceil(input_size) * 2 + 50
  runtime {
    bootDiskSizeGb: 8
    cpu: 1
    disks: "local-disk ${disk_size} SSD"
    docker: base_docker
    maxRetries: 1
    memory: "4 GiB"
    preemptible: 3
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    plot_tars='~{write_lines(plot_tars)}'
    excess_call_samples='~{write_lines(excess_call_samples)}'
    tar_prefix='~{tar_prefix}'

    mkdir temp
    while read -r f; do
      tar -xf "${f}"
      find plots -type f -name '*.jpg' -exec mv '{}' temp \;
      rm -r plots
    done < "${plot_tars}"

    mv temp "${tar_prefix}"
    tar -cf "${tar_prefix}.tar" "${tar_prefix}"

    cat "${excess_call_samples}" \
      | xargs cat > "${tar_prefix}-excess_gds.tsv"
  >>>
}
