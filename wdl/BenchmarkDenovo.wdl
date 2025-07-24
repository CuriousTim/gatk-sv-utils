version 1.0

# Benchmark de novo SV calling.
workflow BenchmarkDenovo {
  input {
    # See benchmark_denovo.R
    File denovo
    File truth
    # VCFs that were run through the de novo pipeline
    Array[File] vcfs
    # Basically the sample data table, but in a convenient file
    File sample_table
    String base_docker
    String r_docker
  }

  scatter (vcf in vcfs) {
    call ReformatVcf {
      input:
        vcf = vcf,
        base_docker = base_docker
    }
  }

  call Benchmark {
    input:
      denovo = denovo,
      truth = truth,
      input_vcf_tsvs = ReformatVcf.vcf_tsv,
      sample_table = sample_table,
      r_docker = r_docker
  }

  output {
    File truth_in_input_vcf = Benchmark.truth_in_input_vcf
    File truth_not_in_input_vcf = Benchmark.truth_not_in_input_vcf
    File denovo_in_truth = Benchmark.denovo_in_truth
    File denovo_not_in_truth = Benchmark.denovo_not_in_truth
    File not_denovo_in_truth = Benchmark.not_denovo_in_truth
    File not_denovo_not_in_truth = Benchmark.not_denovo_not_in_truth
    File truth_in_denovo = Benchmark.truth_in_denovo
    File truth_not_in_denovo = Benchmark.truth_not_in_denovo
    File overlaps_plot = Benchmark.overlaps_plot
    File benchmark_log = Benchmark.benchmark_log
  }
}

task ReformatVcf {
  input {
    File vcf
    String base_docker
  }

  Float disk_size = size(vcf, "GB") * 2 + 16

  runtime {
    bootDiskSizeGb: 8
    cpus: 2
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

    bcftools query -i 'GT="alt" & INFO/SVTYPE != "CNV" & INFO/SVTYPE != "BND"' \
      -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t[%SAMPLE,]\n' \
      "${vcf}" \
      | gawk -F'\t' '{sub(/,$/, "", $5); split($5, a, /,/); for(i in a){$5=a[i]; print}}' OFS='\t' \
      | gzip -c > vcf.tsv.gz
  >>>

  output {
    File vcf_tsv = "vcf.tsv.gz"
  }
}

task Benchmark {
  input {
    File denovo
    File truth
    Array[File] input_vcf_tsvs
    File sample_table
    String r_docker
  }

  Float disk_size = size([denovo, truth], "GB") * 2
    + (size(input_vcf_tsvs, "GB") * 2)
    + size(sample_table, "GB") + 16

  runtime {
    bootDiskSizeGb: 8
    cpus: 2
    disks: "local-disk ${ceil(disk_size)} HDD"
    docker: r_docker
    maxRetries: 1
    memory: "4 GiB"
    preemptible: 3
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    denovo='~{denovo}'
    truth='~{truth}'
    input_vcf_tsvs='~{write_lines(input_vcf_tsvs)}'
    sample_table='~{sample_table}'

    cat "${input_vcf_tsvs}" \
      | xargs cat > vcf.tsv.gz

    Rscript /opt/gatk-sv-utils/scripts/benchmark_denovo.R \
      "${denovo}" \
      "${truth}" \
      vcf.tsv.gz \
      "${sample_table}" \
        |& tee benchmark.log
  >>>

  output {
    File truth_in_input_vcf = "truth_in_input_vcf.tsv"
    File truth_not_in_input_vcf = "truth_not_in_input_vcf.tsv"
    File denovo_in_truth = "denovo_in_truth.tsv"
    File denovo_not_in_truth = "denovo_not_in_truth.tsv"
    File not_denovo_in_truth = "not_denovo_in_truth.tsv"
    File not_denovo_not_in_truth = "not_denovo_not_in_truth.tsv"
    File truth_in_denovo = "truth_in_denovo.tsv"
    File truth_not_in_denovo = "truth_not_in_denovo.tsv"
    File overlaps_plot = "benchmark_venn.jpg"
    File benchmark_log = "benchmark.log"
  }
}
