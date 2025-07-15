version 1.0

# Benchmark de novo SV calling.
workflow BenchmarkDenovo {
  input {
    # See benchmark_denovo.R
    File denovo
    File truth
    File cleanvcf_tsv
    # Basically the sample data table, but in a convenient file
    File sample_table
    String r_docker
  }

  call Benchmark {
    input:
      denovo = denovo,
      truth = truth,
      cleanvcf_tsv = cleanvcf_tsv,
      sample_table = sample_table,
      r_docker = r_docker
  }

  output {
    File truth_in_cleanvcf = Benchmark.truth_in_cleanvcf
    File truth_not_in_cleanvcf = Benchmark.truth_not_in_cleanvcf
    File denovo_in_truth = Benchmark.denovo_in_truth
    File denovo_not_in_truth = Benchmark.denovo_not_in_truth
    File not_denovo_in_truth = Benchmark.not_denovo_in_truth
    File not_denovo_not_in_truth = Benchmark.not_denovo_not_in_truth
    File truth_in_denovo = Benchmark.truth_in_denovo
    File truth_not_in_denovo = Benchmark.truth_not_in_denovo
  }
}

task Benchmark {
  input {
    File denovo
    File truth
    File cleanvcf_tsv
    File sample_table
    String r_docker
  }

  Float disk_size = size([denovo, truth, cleanvcf_tsv], "GB") * 2 + size(sample_table, "GB") + 16

  runtime {
    bootDiskSizeGb: 8
    cpus: 2
    disks: "local-disk ${disk_size} HDD"
    docker: r_docker
    maxRetries: 1
    memory: "4 GiB"
    preemptible: 3
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    denovo="~{denovo}"
    truth="~{truth}"
    cleanvcf_tsv="~{cleanvcf_tsv}"
    sample_table="~{sample_table}"

    Rscript /opt/gatk-sv-utils/scripts/benchmark_denovo.R \
      "${denovo}" \
      "${truth}" \
      "${cleanvcf_tsv}" \
      "${sample_table}"
  >>>

  output {
    File truth_in_cleanvcf = "truth_in_cleanvcf.tsv"
    File truth_not_in_cleanvcf = "truth_not_in_cleanvcf.tsv"
    File denovo_in_truth = "denovo_in_truth.tsv"
    File denovo_not_in_truth = "denovo_not_in_truth.tsv"
    File not_denovo_in_truth = "not_denovo_in_truth.tsv"
    File not_denovo_not_in_truth = "not_denovo_not_in_truth.tsv"
    File truth_in_denovo = "truth_in_denovo.tsv"
    File truth_not_in_denovo = "truth_not_in_denovo.tsv"
  }
}
