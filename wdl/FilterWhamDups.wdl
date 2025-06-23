version 1.0

# Filter Wham-only DUPs
workflow FilterWhamDups {
  input {
    Array[File] vcfs
    File blacklist_regions
    # Wham-only DUPs matching `extra_filters` will be checked for overlap
    # against `blacklist_regions`. The ones overlapping the regions will
    # be removed from the VCF.
    String? extra_filters
    String runtime_docker
  }

  scatter (vcf in vcfs) {
    call FilterSites {
      input:
        vcf = vcf,
        blacklist = blacklist_regions,
        extra_filters = extra_filters,
        runtime_docker = runtime_docker
    }
  }

  output {
    Array[File] filtered_vcfs = FilterSites.filtered_vcf
    Array[File] filtered_vcf_indicies = FilterSites.filtered_vcf_index
  }
}

task FilterSites {
  input {
    File vcf
    File blacklist
    String extra_filters = '(EVIDENCE == "SR" || EVIDENCE ~ "^RD,SR$")'
    String runtime_docker

    Float? memory_gib
    Int? boot_disk_gb
    Int? cpus
    Int? disk_gb
    Int? max_retries
    Int? preemptible_tries
  }

  Float disk_size = size(vcf, "GB") * 2.2  + size(blacklist, "GB") * 2 + 16
  String output_vcf = "filtered-${basename(vcf)}"
  String output_vcf_index = "${output_vcf}.tbi"

  runtime {
    bootDiskSizeGb: select_first([boot_disk_gb, 8])
    cpus: select_first([cpus, 1])
    disks: "local-disk ${select_first([disk_gb, ceil(disk_size)])} HDD"
    docker: runtime_docker
    maxRetries: select_first([max_retries, 1])
    memory: "${select_first([memory_gib, 2])} GiB"
    preemptible: select_first([preemptible_tries, 3])
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    in_vcf='~{vcf}'
    blacklist='~{blacklist}'
    out_vcf='~{output_vcf}'
    extra_filters='~{if defined(extra_filters) then extra_filters else ""}'

    cat2() {
      if [[ "$1" == *.gz ]]; then
        gzip -cd "$1"
      else
        cat "$1"
      fi
    }

    vcf_sites="$(mktemp -p "${PWD}" tmp.XXXXXXXXX)"
    bcftools view --include 'SVTYPE=="DUP" && ALGORITHMS=="wham"' --output-type u "${in_vcf}" \
      | bcftools query "${extra_filters:+--include}" "${extra_filters:-}" \
          --format '%CHROM\t%POS0\t%INFO/END\t%ID\n' \
      | LC_ALL=C sort -k1,1 -k2,2n > "${vcf_sites}"
    read -r wham_only_count _ < <(wc -l "${vcf_sites}")
    printf '%d Wham-only DUPs passing initial filters\n' "${wham_only_count}" >&2

    bl_sorted="$(mktemp -p "${PWD}" tmp.XXXXXXXXX)"
    cat2 "${blacklist}" | LC_ALL=C sort -k1,1 -k2,2n > "${bl_sorted}"

    bad_sites="$(mktemp -p "${PWD}" tmp.XXXXXXXXX)"
    bedtools intersect -a "${bl_sorted}" -b "${vcf_sites}" -wb -F 0.5 -sorted \
      | cut -f 7 \
      | LC_ALL=C sort -u > "${bad_sites}"

    trap 'rm -f "${vcf_sites}" "${bl_sorted}" "${bad_sites}"' EXIT
    read -r bad_sites_count _ < <(wc -l "${bad_sites}")
    if (( bad_sites_count > 0 )); then
      printf 'removing %d Wham-only DUPs overlapping blacklist regions\n' "${bad_sites_count}" >&2
    else
      printf 'no Wham-only DUPs overlapping blacklist regions with >=0.5 coverage of SV\n' >&2
      cp "${in_vcf}" "${out_vcf}"
      if [[ ! -f "${out_vcf}.tbi" ]]; then
        bcftools index --tbi "${out_vcf}"
      fi
      exit 0
    fi

    bcftools filter --exclude "ID=@${bad_sites}" --output "${out_vcf}" \
      --output-type z "${in_vcf}" --write-index=tbi
  >>>

  output {
    File filtered_vcf = "${output_vcf}"
    File filtered_vcf_index = "${output_vcf_index}"
  }
}
