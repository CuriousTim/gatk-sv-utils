version 1.0

# Benchmark de novo SV calling.
workflow BenchmarkDenovo {
  input {
    # See benchmark_denovo.R
    File denovos
    # The "true" de novo calls. Must be subset to the samples of interest.
    File truth_vcf
    File truth_vcf_index
    # VCFs that were run through the de novo pipeline
    Array[File] start_vcfs
    File reference_dict

    String base_docker
    String gatk_docker
  }

  scatter (vcf in start_vcfs) {
    call SubsetStartVcf {
      input:
        start_vcf = vcf,
        truth_vcf = truth_vcf,
        denovos = denovos,
        base_docker = base_docker
    }

    call SVConcordance {
      input:
        eval_vcf = SubsetStartVcf.denovo_subset_vcf,
        eval_vcf_index = SubsetStartVcf.denovo_subset_vcf_index,
        truth_vcf = truth_vcf,
        truth_vcf_index = truth_vcf_index,
        start_vcf = SubsetStartVcf.sample_subset_vcf,
        start_vcf_index = SubsetStartVcf.sample_subset_vcf_index,
        reference_dict = reference_dict,
        gatk_docker = gatk_docker
    }
  }

  output {
    Array[File] truth_in_eval_vcfs = SVConcordance.truth_in_eval_vcf
    Array[File] eval_in_truth_vcfs = SVConcordance.eval_in_truth_vcf
    Array[File] truth_in_start_vcfs = SVConcordance.truth_in_start_vcf
  }
}

task SubsetStartVcf {
  input {
    File start_vcf
    File truth_vcf
    File denovos
    String base_docker
  }

  Float disk_size = size([start_vcf, truth_vcf, denovos], "GB") * 2 + 16

  runtime {
    bootDiskSizeGb: 8
    cpus: 1
    disks: "local-disk ${ceil(disk_size)} HDD"
    docker: base_docker
    maxRetries: 1
    memory: "4 GiB"
    preemptible: 3
  }

  String sample_subset_name = "start-${basename(start_vcf)}"
  String denovo_subset_name = "denovo-${basename(start_vcf)}"

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    bcftools query --list-samples '~{truth_vcf}' > truth_samples
    gzip -cd '~{denovos}' \
      | awk -F'\t' 'function check(col) {
          if (!(col in Header)) {
            print "required column \047" col "\047 is missing" > "/dev/stderr"
            exit 1
          }
        }
        NR==1{for(i=1;i<=NF;++i){Header[$i]=i}}
        NR==2{check("name");check("sample");check("is_de_novo")}
        NR>1 && $(Header["is_de_novo"]) == "TRUE"{print $(Header["name"])"\t"$(Header["sample"])}' \
      | LC_ALL=C sort -u > denovo_vids
    if [[ ! -s denovo_vids ]]; then
      printf '%s\n' 'de novo output does not have any de novos' >&2
      exit 1
    fi

    bcftools view --samples-file truth_samples --output '~{sample_subset_name}' \
      --output-type z '~{start_vcf}'
    bcftools index --tbi '~{sample_subset_name}'
    set_denovo_gt '~{sample_subset_name}' denovo_vids '~{denovo_subset_name}'
    bcftools index --tbi '~{denovo_subset_name}'
  >>>

  output {
    File denovo_subset_vcf = denovo_subset_name
    File denovo_subset_vcf_index = "${denovo_subset_name}.tbi"
    File sample_subset_vcf = sample_subset_name
    File sample_subset_vcf_index = "${sample_subset_name}.tbi"
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

  Float disk_size = size([eval_vcf, truth_vcf], "GB") * 2 + 16

  runtime {
    bootDiskSizeGb: 8
    cpus: 1
    disks: "local-disk ${ceil(disk_size)} HDD"
    docker: gatk_docker
    maxRetries: 1
    memory: "8 GiB"
    preemptible: 3
  }

  String truth_in_eval_name = "truth_in_eval-${basename(truth_vcf)}"
  String eval_in_truth_name = "eval_in_truth-${basename(eval_vcf)}"
  String truth_in_start_name = "truth_in_start-${basename(truth_vcf)}"

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

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
    File truth_in_eval_vcf = "${truth_in_eval_name}"
    File eval_in_truth_vcf = "${eval_in_truth_name}"
    File truth_in_start_vcf = "${truth_in_start_name}"
  }
}
