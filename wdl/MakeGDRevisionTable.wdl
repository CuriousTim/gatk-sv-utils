version 1.0

workflow MakeGDRevisionTable {
  input {
    File gd_table
    Array[File] vcfs
    Array[File] vcf_indices
    String output_prefix
    String base_docker
    String python_docker
  }

  output {
    File remove_vids = Merge.merged_remove_vids
    File missing = Merge.merged_missing
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

  call Merge {
    input:
      remove_vids = CheckVCFForGDs.remove_vids,
      missing = CheckVCFForGDs.missing,
      output_prefix = output_prefix,
      base_docker = base_docker
  }
}

task CheckVCFForGDs {
  input {
    File gd_table
    File vcf
    File vcf_index
    Int shard_index
    String python_docker
  }

  Float disk_size = size([gd_table, vcf, vcf_index], "GB") * 2 + 16

  runtime {
    bootDiskSizeGb: 8
    cpus: 1
    disks: "local-disk ${ceil(disk_size)}  HDD"
    docker: python_docker
    maxRetries: 1
    memory: "4 GiB"
    preemptible: 3
  }

  output {
    File remove_vids = "remove_vids-${shard_index}.txt"
    File missing = "missing-${shard_index}.tsv"
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    gd_table='~{gd_table}'
    vcf='~{vcf}'
    shard_index='~{shard_index}'

    # the script silently skips GDs on contigs not in the VCF
    python3 /opt/gatk-sv-utils/script/compare_gds_to_vcf.py "${gd_table}" \
      "${vcf}" "remove_vids-${shard_index}.txt" "missing-${shard_index}.tsv"
  >>>
}

task Merge {
  input {
    String output_prefix
    Array[File] remove_vids
    Array[File] missing
    String base_docker
  }

  Float disk_size = (size(remove_vids, "GB") + size(missing, "GB")) * 2 + 16

  runtime {
    bootDiskSizeGb: 8
    cpus: 1
    disks: "local-disk ${ceil(disk_size)}  HDD"
    docker: base_docker
    maxRetries: 1
    memory: "1 GiB"
    preemptible: 3
  }

  output {
    File merged_remove_vids = "${output_prefix}-gd_remove_vids.txt"
    File merged_missing = "${output_prefix}-gd_missing.txt"
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    remove_vids='~{write_lines(remove_vids)}'
    missing='~{write_lines(missing)}'
    output_prefix='~{output_prefix}'

    cat "${remove_vids}" | xargs cat > "${output_prefix}-gd_remove_vids.txt"
    cat "${missing}" | xargs cat > "${output_prefix}-gd_missing.txt"
  >>>
}
