version 1.0

# Benchmark de novo SV calling.
workflow BenchmarkDenovo {
  input {
    # See benchmark_denovo.R
    File denovos
    File sample_table
    # The "true" de novo calls
    File truth_vcf
    File truth_vcf_index
    # VCFs that were run through the de novo pipeline split by contig
    Array[File]+ start_vcfs
    Array[File]+ start_vcf_indices
    File contigs
    File primary_contigs_fai
    File reference_dict

    String base_docker
    String gatk_docker
  }

  Array[String] contigs_arr = read_lines(contigs)

  call GetSharedSamples {
    input:
      start_vcf = start_vcfs[0],
      truth_vcf = truth_vcf,
      sample_table = sample_table,
      base_docker = base_docker
  }

  scatter (i in range(length(contigs_arr))) {
    call SubsetVcf as subset_start {
      input:
        vcf = start_vcfs[i],
        vcf_index = start_vcf_indices[i],
        sample_ids = GetSharedSamples.shared_samples,
        primary_contigs_fai = primary_contigs_fai,
        base_docker = base_docker
    }

    call SubsetVcf as subset_truth {
      input:
        vcf = truth_vcf,
        vcf_index = truth_vcf_index,
        sample_ids = GetSharedSamples.shared_samples,
        contig = contigs_arr[i],
        primary_contigs_fai = primary_contigs_fai,
        base_docker = base_docker
    }

    call MakeDenovoVcf {
      input:
        subset_start_vcf = subset_start.subset_vcf,
        denovos = denovos,
        base_docker = base_docker
    }

    call SVConcordance {
      input:
        eval_vcf = MakeDenovoVcf.denovo_vcf,
        eval_vcf_index = MakeDenovoVcf.denovo_vcf_index,
        truth_vcf = subset_truth.subset_vcf,
        truth_vcf_index = subset_truth.subset_vcf_index,
        start_vcf = subset_start.subset_vcf,
        start_vcf_index = subset_start.subset_vcf_index,
        reference_dict = reference_dict,
        gatk_docker = gatk_docker
    }
  }

  call ConcatVcfs as concat_eval_in_truth {
    input:
      vcfs = SVConcordance.eval_in_truth_vcf,
      concat_prefix = "eval_in_truth",
      base_docker = base_docker
  }

  call ConcatVcfs as concat_truth_in_eval {
    input:
      vcfs = SVConcordance.truth_in_eval_vcf,
      concat_prefix = "truth_in_eval",
      base_docker = base_docker
  }

  call ConcatVcfs as concat_truth_in_start {
    input:
      vcfs = SVConcordance.truth_in_start_vcf,
      concat_prefix = "truth_in_start",
      base_docker = base_docker
  }

  output {
    File eval_in_truth_vcf = concat_eval_in_truth.concat_vcf
    File eval_in_truth_vcf_index = concat_eval_in_truth.concat_vcf_index
    File truth_in_eval_vcf = concat_truth_in_eval.concat_vcf
    File truth_in_eval_vcf_index = concat_truth_in_eval.concat_vcf_index
    File truth_in_start_vcf = concat_truth_in_start.concat_vcf
    File truth_in_start_vcf_index = concat_truth_in_start.concat_vcf_index
  }
}

task GetSharedSamples {
  input {
    File start_vcf
    File truth_vcf
    File sample_table
    String base_docker
  }

  Float disk_size = size([start_vcf, truth_vcf, sample_table], "GB") * 2 + 16

  runtime {
    bootDiskSizeGb: 8
    cpus: 1
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

    start_vcf='~{start_vcf}'
    truth_vcf='~{truth_vcf}'
    sample_table='~{sample_table}'

    mv "${sample_table}" samples.tsv
    duckdb ':memory:' "COPY (SELECT DISTINCT \"entity:sample_id\" FROM 'samples.tsv' WHERE cohort_short = 'SSC') TO 'ssc' (HEADER false);"
    bcftools query --list-samples "${start_vcf}" | sort > start_samples
    bcftools query --list-samples "${truth_vcf}" | sort > truth_samples
    comm -12 ssc start_samples > temp
    comm -12 temp truth_samples > shared_samples.list
    if [[ ! -s 'shared_samples.list' ]]; then
      printf 'Start VCF and truth VCF have 0 shared samples\n' >&2
      exit 1
    fi
  >>>

  output {
   File shared_samples = "shared_samples.list"
  }
}

task SubsetVcf {
  input {
    File vcf
    File vcf_index
    File sample_ids
    File primary_contigs_fai
    String? contig
    String base_docker
  }

  Float disk_size = size([vcf, sample_ids, primary_contigs_fai], "GB") * 2 + 16

  runtime {
    bootDiskSizeGb: 8
    cpus: 1
    disks: "local-disk ${ceil(disk_size)} HDD"
    docker: base_docker
    maxRetries: 1
    memory: "4 GiB"
    preemptible: 3
  }

  String subset_vcf_name = "subset-${basename(vcf)}"

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    vcf='~{vcf}'
    sample_ids='~{sample_ids}'
    primary_contigs_fai='~{primary_contigs_fai}'
    contig='~{if defined(contig) then "--regions ${contig}" else ""}'
    subset_vcf_name='~{subset_vcf_name}'

    bcftools view --samples-file "${sample_ids}" --output "${subset_vcf_name}" \
      --output-type z ${contig} "${vcf}"
    bcftools reheader --fai "${primary_contigs_fai}" "${subset_vcf_name}"
    bcftools index --tbi "${subset_vcf_name}"
  >>>

  output {
    File subset_vcf = subset_vcf_name
    File subset_vcf_index = "${subset_vcf_name}.tbi"
  }
}

task MakeDenovoVcf {
  input {
    File subset_start_vcf
    File denovos
    String base_docker
  }

  Float disk_size = size([subset_start_vcf, denovos], "GB") * 2 + 16

  runtime {
    bootDiskSizeGb: 8
    cpus: 1
    disks: "local-disk ${ceil(disk_size)} HDD"
    docker: base_docker
    maxRetries: 1
    memory: "4 GiB"
    preemptible: 3
  }

  String denovo_vcf_name = "denovo-${basename(subset_start_vcf)}"

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    subset_start_vcf='~{subset_start_vcf}'
    denovos='~{denovos}'
    denovo_vcf_name='~{denovo_vcf_name}'

    mv "${denovos}" denovos.tsv.gz
    cat > commands.sql <<EOF
    COPY (
      SELECT DISTINCT "name", "sample"
      FROM 'denovos.tsv.gz'
      WHERE is_de_novo = 'TRUE'
    ) TO 'denovo_vids' (HEADER false, DELIMITER '\t');
    EOF
    duckdb ':memory:' < commands.sql
    if [[ ! -s denovo_vids ]]; then
      printf '%s\n' 'de novo output does not have any de novos' >&2
      exit 1
    fi

    set_denovo_gt "${subset_start_vcf}" denovo_vids "${denovo_vcf_name}"
    bcftools index --tbi "${denovo_vcf_name}"
  >>>

  output {
    File denovo_vcf = denovo_vcf_name
    File denovo_vcf_index = "${denovo_vcf_name}.tbi"
  }
}

task SVConcordance {
  input {
    File eval_vcf
    File eval_vcf_index
    File truth_vcf
    File truth_vcf_index
    File start_vcf
    File start_vcf_index
    File reference_dict
    String gatk_docker
  }

  Float disk_size = size([eval_vcf, truth_vcf, start_vcf], "GB") * 3 + 16

  runtime {
    bootDiskSizeGb: 8
    cpus: 1
    disks: "local-disk ${ceil(disk_size)} HDD"
    docker: gatk_docker
    maxRetries: 1
    memory: "8 GiB"
    preemptible: 3
  }

  String eval_in_truth_name = "${basename(eval_vcf)}"
  String truth_in_eval_name = "${basename(truth_vcf)}"
  String truth_in_start_name = "${basename(truth_vcf)}"

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    eval_vcf='~{eval_vcf}'
    truth_vcf='~{truth_vcf}'
    start_vcf='~{start_vcf}'
    reference_dict='~{reference_dict}'

    gatk --java-options '-Xmx8000M' SVConcordance \
      --keep-all \
      --sequence-dictionary '~{reference_dict}' \
      --eval '~{eval_vcf}' \
      --truth '~{truth_vcf}' \
      --output '~{eval_in_truth_name}'
    gatk --java-options '-Xmx8000M' SVConcordance \
      --keep-all \
      --sequence-dictionary '~{reference_dict}' \
      --eval '~{truth_vcf}' \
      --truth '~{eval_vcf}' \
      --output '~{truth_in_eval_name}'
    gatk --java-options '-Xmx8000M' SVConcordance \
      --keep-all \
      --sequence-dictionary '~{reference_dict}' \
      --eval '~{truth_vcf}' \
      --truth '~{start_vcf}' \
      --output '~{truth_in_start_name}'
  >>>

  output {
    File eval_in_truth_vcf = "${eval_in_truth_name}"
    File truth_in_eval_vcf = "${truth_in_eval_name}"
    File truth_in_start_vcf = "${truth_in_start_name}"
  }
}

task ConcatVcfs {
  input {
    Array[File] vcfs
    String concat_prefix
    String base_docker
  }

  Float disk_size = size(vcfs, "GB") * 2 + 16

  runtime {
    bootDiskSizeGb: 8
    cpus: 1
    disks: "local-disk ${ceil(disk_size)} HDD"
    docker: base_docker
    maxRetries: 1
    memory: "4 GiB"
    preemptible: 3
  }

  String concat_vcf_name = "${concat_prefix}.vcf.gz"

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    vcfs='~{write_lines(vcfs)}'

    while read -r f; do
      bcftools index "${f}"
    done < "${vcfs}"

    bcftools concat --file-list "${vcfs}" --output "${concat_vcf_name}" \
      --output-type z --write-index=tbi
  >>>

  output {
    File concat_vcf = concat_vcf_name
    File concat_vcf_index = "${concat_vcf}.tbi"
  }
}
