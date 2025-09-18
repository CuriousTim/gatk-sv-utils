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
    # The contigs must be in the same order as the VCFs
    Array[String] contigs = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
          "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
          "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"]
    File primary_contigs_fai
    File reference_dict

    Float small_cnv_reciprocal_ovp = 0.1
    Float small_cnv_size_sim = 0
    Int small_cnv_breakend_win = 300
    Float large_cnv_reciprocal_ovp = 0.8
    Float large_cnv_size_sim = 0
    Int large_cnv_breakend_win = 10000000
    Float inv_reciprocal_ovp = 0.1
    Float inv_size_sim = 0
    Int inv_breakend_win = 300
    Float ins_reciprocal_ovp = 0
    Float ins_size_sim = 0
    Int ins_breakend_win = 300

    String base_docker
    String gatk_docker
  }

  call GetSharedSamples {
    input:
      start_vcf = start_vcfs[0],
      truth_vcf = truth_vcf,
      sample_table = sample_table,
      base_docker = base_docker
  }

  scatter (i in range(length(contigs))) {
    call SubsetVcf as subset_start {
      input:
        vcf = start_vcfs[i],
        vcf_index = start_vcf_indices[i],
        sample_ids = GetSharedSamples.shared_samples,
        contig = contigs[i],
        primary_contigs_fai = primary_contigs_fai,
        base_docker = base_docker
    }

    call SubsetVcf as subset_truth {
      input:
        vcf = truth_vcf,
        vcf_index = truth_vcf_index,
        sample_ids = GetSharedSamples.shared_samples,
        contig = contigs[i],
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
        small_cnv_reciprocal_ovp = small_cnv_reciprocal_ovp,
        small_cnv_size_sim = small_cnv_size_sim,
        small_cnv_breakend_win = small_cnv_breakend_win,
        large_cnv_reciprocal_ovp = large_cnv_reciprocal_ovp,
        large_cnv_size_sim = large_cnv_size_sim,
        large_cnv_breakend_win = large_cnv_breakend_win,
        inv_reciprocal_ovp = inv_reciprocal_ovp,
        inv_size_sim = inv_size_sim,
        inv_breakend_win = inv_breakend_win,
        ins_reciprocal_ovp = ins_reciprocal_ovp,
        ins_size_sim = ins_size_sim,
        ins_breakend_win = ins_breakend_win,
        gatk_docker = gatk_docker
    }

    call CountConcordance {
      input:
        eval_in_truth_vcf = SVConcordance.eval_in_truth_vcf,
        truth_in_eval_vcf = SVConcordance.truth_in_eval_vcf,
        truth_in_start_vcf =  SVConcordance.truth_in_start_vcf,
        start_vcf = subset_start.subset_vcf,
        base_docker = base_docker
    }
  }

  output {
    Array[File] eval_in_truth_counts = CountConcordance.eval_in_truth_counts
    Array[File] truth_in_eval_counts = CountConcordance.truth_in_eval_counts
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
    String contig
    File primary_contigs_fai
    String base_docker
  }

  Float disk_size = size([vcf, sample_ids], "GB") * 3 + 16

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
    contig='~{contig}'
    primary_contigs_fai='~{primary_contigs_fai}'
    subset_vcf_name='~{subset_vcf_name}'

    # CPX, CTX, CNV, and BND are excluded from benchmarking
    # INFO/ALGORITHMS field is set to pesr for all sites so sites will not be
    # matched according to algorithm, but SVConcordance will accept the VCF
    bcftools view --samples-file "${sample_ids}" --regions "${contig}" \
      --exclude 'INFO/SVTYPE == "CPX" || INFO/SVTYPE == "CTX" || INFO/SVTYPE == "CNV" || INFO/SVTYPE == "BND"' \
      "${vcf}" \
      | gawk -f /opt/gatk-sv-utils/scripts/set_vcf_algorithms.awk - \
      | bgzip -c > temp.vcf.gz

    bcftools head temp.vcf.gz \
      | awk '/^#CHROM/{print > "samples"; next} /^##contig/{next} {print > "newheader"}'

    gawk -F'\t' '{print "##contig=<ID="$1",length="$2">"}' \
      "${primary_contigs_fai}" >> newheader
    cat samples >> newheader

    bcftools reheader --header newheader --output "${subset_vcf_name}" temp.vcf.gz
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

    Float small_cnv_reciprocal_ovp
    Float small_cnv_size_sim
    Int small_cnv_breakend_win
    Float large_cnv_reciprocal_ovp
    Float large_cnv_size_sim
    Int large_cnv_breakend_win
    Float inv_reciprocal_ovp
    Float inv_size_sim
    Int inv_breakend_win
    Float ins_reciprocal_ovp
    Float ins_size_sim
    Int ins_breakend_win

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

    printf 'NAME\tSVTYPE\tMIN_SIZE\tMAX_SIZE\tTRACKS\n' > stratify.tsv
    printf 'DEL_small\tDEL\t-1\t5000\tNULL\n' >> stratify.tsv
    printf 'DUP_small\tDUP\t-1\t5000\tNULL\n' >> stratify.tsv
    printf 'DEL_large\tDEL\t5000\t-1\tNULL\n' >> stratify.tsv
    printf 'DUP_large\tDUP\t5000\t-1\tNULL\n' >> stratify.tsv
    printf 'INV\tINV\t-1\t-1\tNULL\n' >> stratify.tsv
    printf 'INS\tINS\t-1\t-1\tNULL\n' >> stratify.tsv

    printf 'NAME\tRECIPROCAL_OVERLAP\tSIZE_SIMILARITY\tBREAKEND_WINDOW\tSAMPLE_OVERLAP\n' > cluster.tsv
    printf 'DEL_small\t%0.1f\t%0.1f\t%d\t0\n' ~{small_cnv_reciprocal_ovp} ~{small_cnv_size_sim} ~{small_cnv_breakend_win} >> cluster.tsv
    printf 'DUP_small\t%0.1f\t%0.1f\t%d\t0\n' ~{small_cnv_reciprocal_ovp} ~{small_cnv_size_sim} ~{small_cnv_breakend_win} >> cluster.tsv
    printf 'DEL_large\t%0.1f\t%0.1f\t%d\t0\n' ~{large_cnv_reciprocal_ovp} ~{large_cnv_size_sim} ~{large_cnv_breakend_win} >> cluster.tsv
    printf 'DUP_large\t%0.1f\t%0.1f\t%d\t0\n' ~{large_cnv_reciprocal_ovp} ~{large_cnv_size_sim} ~{large_cnv_breakend_win} >> cluster.tsv
    printf 'INV\t%0.1f\t%0.1f\t%d\t0\n' ~{inv_reciprocal_ovp} ~{inv_size_sim} ~{inv_breakend_win} >> cluster.tsv
    printf 'INS\t%0.1f\t%0.1f\t%d\t0\n' ~{ins_reciprocal_ovp} ~{ins_size_sim} ~{ins_breakend_win} >> cluster.tsv

    gatk --java-options '-Xmx8000M' SVConcordance \
      --keep-all \
      --sequence-dictionary "${reference_dict}" \
      --eval "${eval_vcf}" \
      --truth "${truth_vcf}" \
      --output '~{eval_in_truth_name}'
      --stratify-config stratify.tsv \
      --clustering-config cluster.tsv
    gatk --java-options '-Xmx8000M' SVConcordance \
      --eval "${truth_vcf}" \
      --truth "${eval_vcf}" \
      --output '~{truth_in_eval_name}' \
      --keep-all \
      --sequence-dictionary "${reference_dict}" \
      --stratify-config stratify.tsv \
      --clustering-config cluster.tsv
    gatk --java-options '-Xmx8000M' SVConcordance \
      --eval "${truth_vcf}" \
      --truth "${start_vcf}" \
      --output '~{truth_in_start_name}'
      --keep-all \
      --sequence-dictionary "${reference_dict}" \
      --stratify-config stratify.tsv \
      --clustering-config cluster.tsv
  >>>

  output {
    File eval_in_truth_vcf = eval_in_truth_name
    File truth_in_eval_vcf = truth_in_eval_name
    File truth_in_start_vcf = truth_in_start_name
  }
}

task CountConcordance {
  input {
    File eval_in_truth_vcf
    File truth_in_eval_vcf
    File truth_in_start_vcf
    File start_vcf
    String base_docker
  }

  Float disk_size = size([eval_in_truth_vcf, truth_in_eval_vcf, truth_in_start_vcf, start_vcf], "GB") * 2 + 16

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

    eval_in_truth_vcf='~{eval_in_truth_vcf}'
    truth_in_eval_vcf='~{truth_in_eval_vcf}'
    truth_in_start_vcf='~{truth_in_start_vcf}'
    start_vcf='~{start_vcf}'

    bcftools query --include 'GT="alt"' \
      --format '%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/TRUTH_VID[\t%SAMPLE]\n' \
      "${eval_in_truth_vcf}" > 'eval.tsv'
    bcftools query --include 'GT="alt"' \
      --format '%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/TRUTH_VID[\t%SAMPLE]\n' \
      "${truth_in_eval_vcf}" > 'truth.tsv'
    bcftools query --include 'INFO/TRUTH_VID != "."' \
      --format '%ID\t%TRUTH_VID\n' \
      "${truth_in_start_vcf}" > 'truth_in_start.tsv'

    cut -f 2 'truth_in_start.tsv' | tr ',' '\n' | LC_ALL=C sort -u > start_vcf_vids

    bcftools query --include 'ID=@start_vcf_vids & GT="alt"' \
      --format '%ID[\t%SAMPLE]\n' "${start_vcf}" > 'start.tsv'

    gawk -f /opt/gatk-sv-utils/scripts/benchmark_denovo.awk \
      'eval.tsv' 'truth_in_start.tsv' 'truth.tsv' 'start.tsv'
  >>>

  output {
    File eval_in_truth_counts = "eval_vs_truth.tsv"
    File truth_in_eval_counts = "truth_vs_eval.tsv"
  }
}
