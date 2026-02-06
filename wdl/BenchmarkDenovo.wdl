version 1.0

# Benchmark de novo SV calling.
workflow BenchmarkDenovo {
  input {
    # See benchmark_denovo.R
    File denovos
    File sample_table
    File pedigree
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
    File genomic_disorders_bed
    File blacklist_bed
    File sr_bed
    File rm_bed
    File sd_bed
    File gencode_gff

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
    String r_docker
  }

  call GetSharedSamples {
    input:
      start_vcf = start_vcfs[0],
      truth_vcf = truth_vcf,
      pedigree = pedigree,
      sample_table = sample_table,
      base_docker = base_docker
  }

  call FilterTruthVcf {
    input:
      vcf = truth_vcf,
      gd_bed = genomic_disorders_bed,
      bl_bed = blacklist_bed,
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
        vcf = FilterTruthVcf.filtered_vcf,
        vcf_index = FilterTruthVcf.filtered_vcf_index,
        sample_ids = GetSharedSamples.shared_samples,
        contig = contigs[i],
        primary_contigs_fai = primary_contigs_fai,
        base_docker = base_docker
    }

    call MakeDenovoVcf {
      input:
        subset_start_vcf = subset_start.for_concordance_vcf,
        denovos = denovos,
        base_docker = base_docker
    }

    call SVConcordance {
      input:
        eval_vcf = MakeDenovoVcf.denovo_vcf,
        eval_vcf_index = MakeDenovoVcf.denovo_vcf_index,
        truth_vcf = subset_truth.for_concordance_vcf,
        truth_vcf_index = subset_truth.for_concordance_vcf_index,
        start_vcf = subset_start.for_concordance_vcf,
        start_vcf_index = subset_start.for_concordance_vcf_index,
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
        start_vcf = subset_start.for_concordance_vcf,
        base_docker = base_docker
    }
  }

  call MakePlots {
    input:
      eval_in_truth_counts = CountConcordance.eval_in_truth_counts,
      truth_in_eval_counts = CountConcordance.truth_in_eval_counts,
      sr_bed = sr_bed,
      rm_bed = rm_bed,
      sd_bed = sd_bed,
      gencode_gff = gencode_gff,
      r_docker = r_docker
  }

  output {
    File benchmark_plots = MakePlots.benchmark_plots
    Array[File] false_positives = CountConcordance.false_positives
    Array[File] false_negatives_type1 = CountConcordance.false_negatives_type1
    Array[File] false_negatives_type2 = CountConcordance.false_negatives_type2
    Array[File] subset_vcfs = subset_truth.subset_vcf
    Array[File] subset_vcf_indicies = subset_truth.subset_vcf_index
  }
}

task GetSharedSamples {
  input {
    File start_vcf
    File truth_vcf
    File sample_table
    File pedigree
    String base_docker
  }

  Float disk_size = size([start_vcf, truth_vcf, sample_table, pedigree], "GB") * 2 + 16

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
    pedigree='~{pedigree}'

    mv "${sample_table}" samples.tsv
    duckdb ':memory:' "COPY (SELECT DISTINCT \"entity:sample_id\" FROM 'samples.tsv' WHERE cohort_short = 'SSC') TO 'ssc' (HEADER false);"
    bcftools query --list-samples "${start_vcf}" | sort > start_samples
    bcftools query --list-samples "${truth_vcf}" | sort > truth_samples

    # filter by pedigree in case the starting VCF and pedigree are not in sync
    awk -F'\t' '$2 && $3 && $4 {print $2}' "${pedigree}" \
      sort > pedigree_samples
    comm -12 ssc pedigree_samples > step1
    comm -12 step1 start_samples > step2
    comm -12 step2 truth_samples > shared_samples.list
    if [[ ! -s 'shared_samples.list' ]]; then
      printf 'Start VCF and truth VCF have 0 shared samples\n' >&2
      exit 1
    fi
  >>>

  output {
   File shared_samples = "shared_samples.list"
  }
}

task FilterTruthVcf {
  input {
    File vcf
    File gd_bed
    File bl_bed
    String base_docker
  }

  Float disk_size = size(vcf, "GB") * 3 + size([gd_bed, bl_bed], "GB") + 16

  runtime {
    bootDiskSizeGb: 8
    cpus: 1
    disks: "local-disk ${ceil(disk_size)} HDD"
    docker: base_docker
    maxRetries: 1
    memory: "4 GiB"
    preemptible: 3
  }

  String filtered_truth_vcf = "filtered-${basename(vcf)}"

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    vcf='~{vcf}'
    gd_bed='~{gd_bed}'
    bl_bed='~{bl_bed}'
    filtered_truth_vcf='~{filtered_truth_vcf}'

    bcftools query --format '%CHROM\t%POS0\t%END\t%ID\n' \
      --include 'SVTYPE = "DEL" || SVTYPE = "DUP"' "${vcf}" > cnvs.bed
    bedtools intersect -a cnvs.bed -b "${gd_bed}" -r 0.5 -u \
      | cut -f 4 > blacklist
    bcftools query --format '%CHROM\t%POS0\t%END\t%ID\n' "${vcf}" > sites.bed

    bedtools coverage -a sites.bed -b bl_bed \
      | awk -F'\t' '%8 >= 0.5 {print $4}' >> blacklist

    bcftools view --exclude 'ID = @blacklist || SVLEN >= 1000000' --output-type z \
      --output "${filtered_truth_vcf}" --write-index=tbi "${vcf}"
  >>>

  output {
    File filtered_vcf = filtered_truth_vcf
    File filtered_vcf_index = "${filtered_truth_vcf}.tbi"
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
  String for_concordance_vcf_name = "safe-${basename(vcf)}"

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    vcf='~{vcf}'
    sample_ids='~{sample_ids}'
    contig='~{contig}'
    primary_contigs_fai='~{primary_contigs_fai}'
    subset_vcf_name='~{subset_vcf_name}'
    for_concordance_vcf_name='~{for_concordance_vcf_name}'

    # CPX, CTX, CNV, and BND are excluded from benchmarking
    # INFO/ALGORITHMS field is set to pesr for all sites so sites will not be
    # matched according to algorithm, but SVConcordance will accept the VCF
    bcftools view --samples-file "${sample_ids}" --regions "${contig}" \
      --exclude 'INFO/SVTYPE == "BND"' --output "${subset_vcf_name}" \
      --output-type z "${vcf}"
    bcftools index --tbi "${subset_vcf_name}"
    bcftools view \
      --exclude 'INFO/SVTYPE == "CPX" || INFO/SVTYPE == "CTX" || INFO/SVTYPE == "CNV"' \
      "${subset_vcf_name}" \
      | gawk -f /opt/gatk-sv-utils/scripts/set_vcf_algorithms.awk - \
      | bgzip -c > temp.vcf.gz

    bcftools head temp.vcf.gz \
      | awk '/^#CHROM/{print > "samples"; next} /^##contig/{next} {print > "newheader"}'

    gawk -F'\t' '{print "##contig=<ID="$1",length="$2">"}' \
      "${primary_contigs_fai}" >> newheader
    cat samples >> newheader

    bcftools reheader --header newheader --output "${for_concordance_vcf_name}" temp.vcf.gz
    bcftools index --tbi "${for_concordance_vcf_name}"
  >>>

  output {
    File subset_vcf = subset_vcf_name
    File subset_vcf_index = "${subset_vcf_name}.tbi"
    File for_concordance_vcf = for_concordance_vcf_name
    File for_concordance_vcf_index = "${for_concordance_vcf_name}.tbi"
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
      --output '~{eval_in_truth_name}' \
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
      --output '~{truth_in_start_name}' \
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
      --format '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%ID\t%INFO/TRUTH_VID[\t%SAMPLE]\n' \
      "${eval_in_truth_vcf}" > 'eval.tsv'
    bcftools query --include 'GT="alt"' \
      --format '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%ID\t%INFO/TRUTH_VID[\t%SAMPLE]\n' \
      "${truth_in_eval_vcf}" > 'truth.tsv'
    bcftools query --include 'INFO/TRUTH_VID != "."' \
      --format '%ID\t%TRUTH_VID\n' \
      "${truth_in_start_vcf}" > 'truth_in_start.tsv'

    cut -f 2 'truth_in_start.tsv' | tr ',' '\n' | LC_ALL=C sort -u > start_vcf_vids

    bcftools query --include 'ID=@start_vcf_vids & GT="alt"' \
      --format '%ID[\t%SAMPLE]\n' "${start_vcf}" > 'start.tsv'

    : > eval_vs_truth.tsv
    : > truth_vs_eval.tsv
    : > false_positives.tsv
    : > false_negatives-type1.tsv
    : > false_negatives-type2.tsv

    gawk -f /opt/gatk-sv-utils/scripts/benchmark_denovo.awk \
      'eval.tsv' 'truth_in_start.tsv' 'truth.tsv' 'start.tsv'

    gzip false_positives.tsv
    gzip false_negatives-type1.tsv
    gzip false_negatives-type2.tsv
  >>>

  output {
    File eval_in_truth_counts = "eval_vs_truth.tsv"
    File truth_in_eval_counts = "truth_vs_eval.tsv"
    File false_positives = "false_positives.tsv.gz"
    File false_negatives_type1 = "false_negatives-type1.tsv.gz"
    File false_negatives_type2 = "false_negatives-type2.tsv.gz"
  }
}

task MakePlots {
  input {
    Array[File] eval_in_truth_counts
    Array[File] truth_in_eval_counts
    File sr_bed
    File rm_bed
    File sd_bed
    File gencode_gff
    String r_docker
  }

  Float disk_size = size(eval_in_truth_counts, "GB") + size(truth_in_eval_counts, "GB") + 32

  runtime {
    bootDiskSizeGb: 8
    cpus: 1
    disks: "local-disk ${ceil(disk_size)} HDD"
    docker: r_docker
    maxRetries: 1
    memory: "16 GiB"
    preemptible: 3
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    eval_in_truth_counts='~{write_lines(eval_in_truth_counts)}'
    truth_in_eval_counts='~{write_lines(truth_in_eval_counts)}'
    sr_bed='~{sr_bed}'
    rm_bed='~{rm_bed}'
    sd_bed='~{sd_bed}'
    gencode_gff='~{gencode_gff}'

    mv "${sr_bed}" hg38_SR.bed.gz
    mv "${rm_bed}" hg38_RM.bed.gz
    mv "${sd_bed}" hg38_SD.bed.gz
    mv "${gencode_gff}" gencode.v49.basic.annotation.gff3.gz

    cat "${eval_in_truth_counts}" | xargs cat > eval_vs_truth.tsv
    cat "${truth_in_eval_counts}" | xargs cat > truth_vs_eval.tsv

    Rscript /opt/gatk-sv-utils/scripts/benchmark_denovo.R

    mkdir denovo_benchmark
    mv *.jpg denovo_benchmark
    tar -czf denovo_benchmark.tar.gz denovo_benchmark
  >>>

  output {
    File benchmark_plots = "denovo_benchmark.tar.gz"
  }
}
