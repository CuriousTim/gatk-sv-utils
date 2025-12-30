version 1.0

workflow MakeGDRevisionTable {
  input {
    File gd_table
    Array[File] vcfs
    Array[File] vcf_indices
    File sample_table
    File segdups
    File gd_regions
    Array[String] sample_set_ids
    Array[File] bincov_matrices
    Array[File] bincov_matrix_indices
    Array[File] medians_files
    String output_prefix
    String base_docker
    String r_docker
    String python_docker
  }

  output {
    File remove_vids = MergeVCFChecks.merged_remove_vids
    File missing = MergeVCFChecks.merged_missing
    File missing_carrier_plots = MergePlotTars.plots
  }

  call MakeMaps {
    input:
      sample_set_ids = sample_set_ids,
      bincov_matrix_paths = bincov_matrices,
      bincov_matrix_index_paths = bincov_matrix_indices,
      medians_paths = medians_files,
      base_docker = base_docker
  }

  scatter (i in range(length(vcfs))) {
    call CheckVCFForGDs {
      input:
        gd_table = gd_table,
        vcf = vcfs[i],
        vcf_index = vcf_indices[i],
        shard_index = i,
        python_docker = python_docker
    }
  }

  call MergeVCFChecks {
    input:
      remove_vids = CheckVCFForGDs.remove_vids,
      missing = CheckVCFForGDs.missing,
      output_prefix = output_prefix,
      sample_table = sample_table,
      base_docker = base_docker
  }

  scatter (batch in MergeVCFChecks.batches_to_plot) {
    String batch_id = basename(batch)
    call MakePlots {
      input:
        manifest = batch,
        gd_regions = gd_regions,
        bincov_matrix = MakeMaps.bincov_matrix_map[batch_id],
        bincov_matrix_index = MakeMaps.bincov_matrix_index_map[batch_id],
        medians_file = MakeMaps.medians_map[batch_id],
        segdups = segdups,
        batch_id = batch_id,
        r_docker = r_docker
    }
  }

  call MergePlotTars {
    input:
      plot_tars = MakePlots.plots,
      tar_prefix = output_prefix,
      base_docker = base_docker
  }
}

task MakeMaps {
  input {
    Array[String] sample_set_ids
    Array[String] bincov_matrix_paths
    Array[String] bincov_matrix_index_paths
    Array[String] medians_paths
    String base_docker
  }

  parameter_meta {
    sample_set_ids: "Batch IDs."
    bincov_matrix_paths: "Paths to the binned coverage matrix files."
    bincov_matrix_index_paths: "Paths to the binnned coverage matrix index files."
    medians_paths: "Paths to the median coverages files."
    base_docker: "Docker image."
  }

  runtime {
    bootDiskSizeGb: 8
    cpus: 1
    disks: "local-disk 50 HDD"
    docker: base_docker
    maxRetries: 1
    memory: "1 GiB"
    preemptible: 3
  }

  output {
    Map[String, String] bincov_matrix_map = read_map("bincov_matrix_map")
    Map[String, String] bincov_matrix_index_map = read_map("bincov_matrix_index_map")
    Map[String, String] medians_map = read_map("medians_map")
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    sample_set_ids='~{write_lines(sample_set_ids)}'
    bincov_matrix_paths='~{write_lines(bincov_matrix_paths)}'
    bincov_matrix_index_paths='~{write_lines(bincov_matrix_index_paths)}'
    medians_paths='~{write_lines(medians_paths)}'

    paste "${sample_set_ids}" "${bincov_matrix_paths}" > bincov_matrix_map
    paste "${sample_set_ids}" "${bincov_matrix_index_paths}" > bincov_matrix_index_map
    paste "${sample_set_ids}" "${medians_paths}" > medians_map
  >>>
}

task CheckVCFForGDs {
  input {
    File gd_table
    File vcf
    File vcf_index
    Int shard_index
    String python_docker
  }

  parameter_meta {
    gd_table: "File with genomic disorder regions."
    vcf: "VCF to check."
    vcf_index: "Index of VCF to check."
    shard_index: "Index of this shard."
    python_docker: "Docker image."
  }

  output {
    File remove_vids = "remove_vids-${shard_index}.txt"
    File missing = "missing-${shard_index}.tsv"
  }

  Float disk_size = size([gd_table, vcf, vcf_index], "GB") * 2 + 16

  runtime {
    bootDiskSizeGb: 8
    cpus: 1
    disks: "local-disk ${ceil(disk_size)} HDD"
    docker: python_docker
    maxRetries: 1
    memory: "4 GiB"
    preemptible: 3
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    gd_table='~{gd_table}'
    vcf='~{vcf}'
    shard_index='~{shard_index}'

    # the script silently skips GDs on contigs not in the VCF
    python3 /opt/gatk-sv-utils/scripts/compare_gds_to_vcf.py "${gd_table}" \
      "${vcf}" "remove_vids-${shard_index}.txt" "missing-${shard_index}.tsv"
  >>>
}

task MergeVCFChecks {
  input {
    String output_prefix
    Array[File] remove_vids
    Array[File] missing
    File sample_table
    String base_docker
  }

  parameter_meta {
    output_prefix: "The prefix of the merged files."
    remove_vids: "The files with the VCF IDs to remove."
    missing: "The files with discrepencies between the GDs and the VCF."
    sample_table: "Two-column table of batches and samples."
    base_docker: "Docker image."
  }

  output {
    File merged_remove_vids = "${output_prefix}-gd_remove_vids.txt"
    File merged_missing = "${output_prefix}-gd_missing.tsv"
    Array[File] batches_to_plot = glob("batches_to_plot/*")
  }

  Float disk_size = ((size(remove_vids, "GB") + size(missing, "GB")) * 2) + size(sample_table, "GB") + 50

  runtime {
    bootDiskSizeGb: 8
    cpus: 1
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

    remove_vids='~{write_lines(remove_vids)}'
    missing='~{write_lines(missing)}'
    output_prefix='~{output_prefix}'
    sample_table='~{sample_table}'

    cat "${remove_vids}" | xargs cat > "${output_prefix}-gd_remove_vids.txt"
    cat "${missing}" | xargs cat > "${output_prefix}-gd_missing.tsv"

    mkdir batches_to_plot
    cut -f1,6,11 "${output_prefix}-gd_missing.tsv" \
      | gawk -F'\t' 'NR==FNR{a[$2]=$1} NR>FNR{split($3, b, /,/); for (i in b) {if (b[i] in a) {print ($1"\t"$2"\t"b[i]) > ("batches_to_plot/" a[b[i]])}}}' "${sample_table}" -
  >>>
}

task MakePlots {
  input {
    String batch_id
    File manifest
    File segdups
    File gd_regions
    File bincov_matrix
    File bincov_matrix_index
    File medians_file
    String r_docker
  }

  output {
    File plots = "${batch_id}.tar"
  }

  Float disk_size = size(manifest, "GB") + size(gd_regions, "GB") + size(bincov_matrix, "GB") + size(medians_file, "GB") + 50

  runtime {
    bootDiskSizeGb: 8
    cpus: 1
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

    batch_id='~{batch_id}'
    manifest='~{manifest}'
    segdups='~{segdups}'
    gd_regions='~{gd_regions}'
    bincov_matrix='~{bincov_matrix}'
    medians_file='~{medians_file}'

    Rscript /opt/gatk-sv-utils/scripts/visualize_gds.R \
      ${manifest} \
      ${gd_regions} \
      ${segdups} \
      "${bincov_matrix}" \
      "${medians_file}" \
      plots

    tar -cf "${batch_id}.tar" plots
  >>>
}

task MergePlotTars {
  input {
    Array[File] plot_tars
    String tar_prefix
    String base_docker
  }

  parameter_meta {
    plot_tars: "Archives of plots."
    tar_prefix: "Prefix of the output archive."
    base_docker: "Docker image."
  }

  output {
    File plots = "${tar_prefix}.tar"
  }

  Int disk_size = size(plot_tars, "GB") * 2 + 50

  runtime {
    bootDiskSizeGb: 8
    cpu: 1
    disks: "local-disk ${disk_size} SSD"
    docker: base_docker
    maxRetries: 1
    memory: "2 GiB"
    preemptible: 3
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    plot_tars='~{write_lines(plot_tars)}'
    tar_prefix='~{tar_prefix}'

    mkdir temp
    while read -r f; do
      tar -xf "${f}"
      find plots -type f -name '*.jpg' -exec mv '{}' temp \;
      rm -r plots
    done < "${plot_tars}"

    mv temp "${tar_prefix}"
    tar -cf "${tar_prefix}.tar" "${tar_prefix}"
  >>>
}
