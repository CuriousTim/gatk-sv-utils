version 1.0

# Create visualizations for de novo SVs
workflow VisualizeDeNovoSVs {
  input {
    File variants

    # TSV with sample ID, sample set ID
    File sample_table
    File pedigree
    Array[File] merged_pe
    Array[File] merged_pe_index
    Array[File] merged_sr
    Array[File] merged_sr_index
    Array[File] merged_bincov
    Array[File] merged_bincov_index
    Array[File] median_cov
    Array[String] sample_set_id

    File sequence_dict
    String output_tar_prefix

    String base_docker
    String gatk_docker
    String r_docker
  }

  output {
    File denovo_plots_tar = MergePlotsTars.merged_plots_tar
  }

  call BatchVariants {
    input:
      variants = variants,
      sample_table = sample_table,
      sample_set_id = sample_set_id,
      base_docker = base_docker
  }

  scatter (i in range(length(BatchVariants.batched_variants))) {
    if (size(BatchVariants.batched_variants[i]) > 0) {
      call SubsetEvidence {
        input:
          variants = BatchVariants.batched_variants[i],
          merged_pe = merged_pe[i],
          merged_pe_index = merged_pe_index[i],
          merged_sr = merged_sr[i],
          merged_sr_index = merged_sr_index[i],
          merged_bincov = merged_bincov[i],
          merged_bincov_index = merged_bincov_index[i],
          sequence_dict = sequence_dict,
          gatk_docker = gatk_docker
      }

      call MakePlots {
        input:
          variants = BatchVariants.batched_variants[i],
          merged_pe = SubsetEvidence.subset_pe,
          merged_pe_index = SubsetEvidence.subset_pe_index,
          merged_sr = SubsetEvidence.subset_sr,
          merged_sr_index = SubsetEvidence.subset_sr_index,
          merged_bincov = SubsetEvidence.subset_bincov,
          merged_bincov_index = SubsetEvidence.subset_bincov_index,
          median_cov = median_cov[i],
          sample_set_id = sample_set_id[i],
          r_docker = r_docker
      }
    }
  }

  call MergePlotsTars {
    input:
      plots_tars = select_all(MakePlots.plots_tar),
      merged_tar_prefix = output_tar_prefix,
      base_docker = base_docker
  }
}

task BatchVariants {
  input {
    File variants
    File sample_table
    Array[String] sample_set_id
    String base_docker
  }

  output {
    Array[File] batched_variants = glob("batches/*.tsv")
  }

  Float disk_size = size([variants, sample_table], "GB") * 2 + 16

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

    variants='~{variants}'
    sample_table='~{sample_table}'
    sample_set_id='~{write_lines(sample_set_id)}'

    mkdir batches
    gawk -F'\t' '
      ARGIND == 1 {
        path = sprintf("%06d.tsv", FNR)
        printf "" > path
        outpaths[$1] = path
      }
      ARGIND == 2 { sample_map[$1] = $2 }
      ARGIND == 3 && FNR == 1 { next }
      ARGIND == 3 && !($7 in sample_map) {
        printf "%s does not have a sample set\n" > "/dev/stderr"
        exit 1
      }
      ARGIND == 3 {
        sample_set = sample_map[$7]
        path = outpaths[sample_set]
        print > path
      }
    ' "${sample_set_id}" "${sample_table}" "${variants}"
  >>>
}

task SubsetEvidence {
  input {
    File variants
    File merged_pe
    File merged_pe_index
    File merged_sr
    File merged_sr_index
    File merged_bincov
    File merged_bincov_index
    File sequence_dict
    String gatk_docker
  }

  parameter_meta {
    merged_pe: {
      description: "Sample set PE matrix",
      localization_optional: true
    }
    merged_sr: {
      description: "Sample set SR matrix",
      localization_optional: true
    }
    merged_bincov: {
      description: "Sample set RD matrix",
      localization_optional: true
    }
  }

  output {
    File subset_pe = "subset.PE.txt.gz"
    File subset_pe_index = "subset.PE.txt.gz.tbi"
    File subset_sr = "subset.SR.txt.gz"
    File subset_sr_index = "subset.SR.txt.gz.tbi"
    File subset_bincov = "subset.RD.txt.gz"
    File subset_bincov_index = "subset.RD.txt.gz.tbi"
  }

  runtime {
    bootDiskSizeGb: 8
    cpu: 4
    disks: "local-disk 256 HDD"
    docker: gatk_docker
    maxRetries: 1
    memory: "16 GiB"
    preemptible: 3
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    variants="~{variants}"
    merged_pe="~{merged_pe}"
    merged_sr="~{merged_sr}"
    merged_bincov="~{merged_bincov}"
    sequence_dict="~{sequence_dict}"

    # expand all ranges by 50% upstream and downstream so the visualizations
    # can have padding
    awk -F'\t' '{size=$3-$2+1;pad=size / 2;a=$2-pad;b=$3+pad;a=a<1?1:a;print $1"\t"a-1"\t"b}' \
      "${variants}" > padded_coords.bed

    gatk --java-options "Xmx6G" PrintSVEvidence \
      --evidence-file  "${merged_pe}" \
      --interval-merging-rule ALL \
      --intervals padded_coords.bed \
      --sequence-dictionary "${sequence_dict}" \
      --output "subset.PE.txt.gz"
    gatk --java-options "Xmx6G" PrintSVEvidence \
      --evidence-file  "${merged_sr}" \
      --interval-merging-rule ALL \
      --intervals padded_coords.bed \
      --sequence-dictionary "${sequence_dict}" \
      --output "subset.SR.txt.gz"
    gatk --java-options "Xmx6G" PrintSVEvidence \
      --evidence-file  "${merged_bincov}" \
      --interval-merging-rule ALL \
      --intervals padded_coords.bed \
      --sequence-dictionary "${sequence_dict}" \
      --output "subset.RD.txt.gz"
  >>>
}

task MakePlots {
  input {
    File variants
    File merged_pe
    File merged_pe_index
    File merged_sr
    File merged_sr_index
    File merged_bincov
    File merged_bincov_index
    File median_cov
    String sample_set_id
    String r_docker
  }

  output {
    File plots_tar = "${sample_set_id}.tar"
  }

  Float disk_size = size([merged_pe, merged_sr, merged_bincov, median_cov], "GB") * 1.5 + 32

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

    variants="~{variants}"
    merged_pe="~{merged_pe}"
    merged_sr="~{merged_sr}"
    merged_bincov="~{merged_bincov}"
    median_cov="~{median_cov}"
    sample_set_id="~{sample_set_id}"

    Rscript /opt/gatk-sv-utils/script/visualize_denovos.R \
      "${variants}" pe.txt.gz sr.txt.gz rd.txt.gz \
      "${median_cov}" "${sample_set_id}"

    tar -cf "${sample_set_id}.tar" "${sample_set_id}"
  >>>
}

task MergePlotsTars {
  input {
    Array[File] plots_tars
    String merged_tar_prefix
    String base_docker
  }

  output {
    File merged_plots_tar = "${merged_tar_prefix}"
  }

  Float disk_size = size(plots_tars, "GB") * 3 + 16

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

    plots_tars="~{write_lines(plots_tars)}"
    merged_tar_prefix="~{merged_tar_prefix}"

    touch "${merged_tar_prefix}.tar"
    cat "${plots_tars}" \
      | xargs tar --concatenate --file="${merged_tar_prefix}.tar"
  >>>
}
