version 1.0

# Create visualizations for de novo SVs
workflow VisualizeDeNovoSVs {
  input {
    # The variants file can be a TSV in one of two formats:
    # 1. The file has a header line as the first line and there must be columns named chr, start,
    #    end, svlen, vid, svtype, and sample. These columns can be in any order relative to each
    #    other and there can be other columns.
    # 2. The file has a header line as the first line and there are exactly seven columns. The names
    #    can be anything, but they must describe the variant chromosome, start, end, SV length,
    #    variant ID, SV type and sample ID, in that order.
    File variants
    Int variants_per_batch = 25

    # TSV with sample ID, sample set ID
    File sample_table
    File pedigree
    Array[File] merged_pe
    Array[File] merged_pe_index
    Array[File] merged_sr
    Array[File] merged_sr_index
    Array[File] merged_rd
    Array[File] merged_rd_index
    Array[File] median_cov
    Array[String] sample_set_id

    String output_tar_prefix

    String base_docker
    String r_docker
  }

  output {
    File denovo_plots_tar = MergePlotsTars.merged_plots_tar
  }

  call BatchVariants {
    input:
      variants = variants,
      variants_per_batch = variants_per_batch,
      base_docker = base_docker
  }

  call MakeEvidenceManifest {
    input:
      sample_set_id = sample_set_id,
      merged_pe = merged_pe,
      merged_sr = merged_sr,
      merged_rd = merged_rd,
      median_cov = median_cov,
      base_docker = base_docker
  }

  scatter (i in range(length(BatchVariants.batched_variants))) {
    if (size(BatchVariants.batched_variants[i]) > 0) {
      call MakePlots {
        input:
          variants = BatchVariants.batched_variants[i],
          evidence_manifest = MakeEvidenceManifest.evidence_manifest,
          tar_prefix = "batch_${i}",
          sample_table = sample_table,
          pedigree = pedigree,
          r_docker = r_docker,
          max_svlen_file = BatchVariants.max_svlens[i]
      }
    }
  }

  call MergePlotsTars {
    input:
      plots_tars = select_all(MakePlots.plots_tar),
      variants = BatchVariants.batched_variants,
      merged_tar_prefix = output_tar_prefix,
      base_docker = base_docker
  }
}

task BatchVariants {
  input {
    File variants
    Int variants_per_batch
    String base_docker
  }

  output {
    Array[File] batched_variants = glob("batches/*.tsv")
    Array[File] max_svlens = glob("mems/*")
  }

  Float disk_size = size(variants, "GB") * 2 + 16

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
    variants_per_batch='~{variants_per_batch}'

    cat2() {
      local magic_num
      magic_num="$(od -N 2 -t x1 "$1" | awk 'NF>1{$1=""; gsub(/ /, ""); print}')"
      if [[ "${magic_num}" = '1f8b' ]]; then
        gzip -cd "$1"
      else
        cat "$1"
      fi
    }

    mkdir batches
    mkdir mems
    gawk -F'\t' -v n="${variants_per_batch}" '
      BEGIN { OFS = "\t" }
      FNR == 1 {
        for (i = 1; i <= NF; ++i) {
          if ($i == "sample") {
            sample_field = i
          }

          if ($i == "svlen") {
            svlen_field = i
          }
        }
        if (!sample_field || !svlen_field) {
          if (NF != 7) {
            print "variants file does not have \"sample\" and \"svlen\" columns" > "/dev/stderr"
            exit 1
          } else {
            header = "chr\tstart\tend\tsvlen\tvid\tsvtype\tsample_id"
            svlen_field = 4
            sample_field = 7
          }
        } else {
          $sample_field = "sample_id"
          header = $0
        }

        next
      }
      j % n == 0 {
        if (outpath) {
          close(outpath)
          mem_outpath = sprintf("mems/%06d", k - 1)
          print max_svlen > mem_outpath
          close(mem_outpath)
          max_svlen = 0
        }
        outpath = sprintf("batches/%06d.tsv", k++)
        print header > outpath
      }
      {
        print > outpath
        max_svlen = max_svlen >= $svlen_field ? max_svlen : $svlen_field
        ++j
      }
      END {
        if (k > 0) {
          mem_outpath = sprintf("mems/%06d", k - 1)
          print max_svlen > mem_outpath
          close(mem_outpath)
        }
      }
    ' <(cat2 "${variants}")
  >>>
}

task MakeEvidenceManifest {
  input {
    Array[String] sample_set_id
    Array[String] merged_pe
    Array[String] merged_sr
    Array[String] merged_rd
    Array[String] median_cov
    String base_docker
  }

  output {
    File evidence_manifest = "evidence_manifest.tsv"
  }

  runtime {
    bootDiskSizeGb: 8
    cpu: 1
    disks: "local-disk 64 HDD"
    docker: base_docker
    maxRetries: 1
    memory: "2 GiB"
    preemptible: 3
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    sample_set_id='~{write_lines(sample_set_id)}'
    merged_pe='~{write_lines(merged_pe)}'
    merged_sr='~{write_lines(merged_sr)}'
    merged_rd='~{write_lines(merged_rd)}'
    median_cov='~{write_lines(median_cov)}'

    paste "${sample_set_id}" "${merged_pe}" "${merged_sr}" "${merged_rd}" \
      "${median_cov}" > 'evidence_manifest.tsv'
  >>>
}

task MakePlots {
  input {
    File variants
    File evidence_manifest
    File sample_table
    String tar_prefix
    File pedigree
    String r_docker
    File max_svlen_file

    Float? mem_gib
    Int? disk_gb
    Int? cpu
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
  }

  output {
    File plots_tar = "${tar_prefix}.tar"
  }

  Int max_svlen = read_int(max_svlen_file)
  Float default_mem_gib = if max_svlen < 10000000 then 16 else 32

  runtime {
    bootDiskSizeGb: select_first([boot_disk_gb, 8])
    cpu: select_first([cpu, 2])
    disks: "local-disk " + select_first([disk_gb, 128]) + " HDD"
    docker: r_docker
    maxRetries: select_first([max_retries, 1])
    memory: select_first([mem_gib, default_mem_gib]) + " GiB"
    preemptible: select_first([preemptible_tries, 3])
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    variants="~{variants}"
    evidence_manifest="~{evidence_manifest}"
    sample_table="~{sample_table}"
    tar_prefix="~{tar_prefix}"
    pedigree="~{pedigree}"

    Rscript /opt/gatk-sv-utils/scripts/visualize_denovos.R \
      "${variants}" "${pedigree}" "${evidence_manifest}" "${sample_table}" \
      "${tar_prefix}" 'exclusions.tsv'

    mv exclusions.tsv "${tar_prefix}/exclusions.tsv"
    tar -cf "${tar_prefix}.tar" "${tar_prefix}"
  >>>
}

task MergePlotsTars {
  input {
    Array[File] plots_tars
    Array[File] variants
    String merged_tar_prefix
    String base_docker
  }

  output {
    File merged_plots_tar = "${merged_tar_prefix}.tar.gz"
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
    variants="~{write_lines(variants)}"

    mkdir store

    cat "${plots_tars}" | xargs -I '{}' tar -xf '{}' -C store
    find store -type f -name '*.png' \
      | awk -F'/' '{print $NF "\t" $0}' > manifest.tsv

    mkdir -p "${merged_tar_prefix}/"{INS,small_CNV,large_CNV,INV,other,not_denovo}
    gawk -F'\t' -v dest="${merged_tar_prefix}" '
    function quote(x) {
      return "\047" gensub(/\047/, "\047\\\047\047", "g", x) "\047"
    }
    ARGIND == 1 {
      plots[$1] = $2
    }
    ARGIND == 2 {
      fnr = 0
      while ((getline line < $0) > 0) {
        ++fnr
        split(line, a, /\t/)
        if (fnr == 1) {
          for (i in a) {
            h[a[i]] = i
          }
          continue
        }
        fn = a[h["vid"]] "~~" a[h["sample_id"]] "~~" a[h["chr"]] "_" a[h["start"]] "-" a[h["end"]] ".png"
        if (fn in plots) {
          svtype = a[h["svtype"]]
          svlen = a[h["svlen"]]
          if (("is_de_novo" in h) && (a[h["is_de_novo"]] == "false" || a[h["is_de_novo"]] == "FALSE")) {
            to = "not_denovo"
          } else if (svtype == "INS") {
            to = "INS"
          } else if (svtype ~ /DEL|DUP/ && svlen < 5000) {
            to = "small_CNV"
          } else if (svtype ~ /DEL|DUP/ && svlen >= 5000) {
            to = "large_CNV"
          } else if (svtype == "INV") {
            to = "INV"
          } else {
            to = "other"
          }
          print ("mv " quote(plots[fn]) " " quote(dest "/" to)) | "sh"
        }
      }
      close($0)
    }' manifest.tsv "${variants}"
    find store -type f -name 'exclusions.tsv' -exec cat '{}' \; \
      | awk '/^sample_id/{if(!f){print; f=1}; next} 1' > "${merged_tar_prefix}/exclusions.tsv"

    tar -czf "${merged_tar_prefix}.tar.gz" "${merged_tar_prefix}"
  >>>
}
