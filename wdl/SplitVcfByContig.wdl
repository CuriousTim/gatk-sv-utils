version 1.0

# Split a VCF by contig
workflow SplitVcfByContig {
  input {
    File vcf
    File vcf_index
    # One of `contig_list` or `contig_list_file` must be given. `contig_list` overrides
    # `contig_list_file`. There will be one VCF generated for each contig in the list.
    Array[String]? contig_list
    File? contig_list_file
    String runtime_docker
  }

  Array[String] contigs = select_first([contig_list, read_lines(select_first([contig_list_file]))])
  scatter (i in range(length(contigs))) {
    call GetContigFromVcf {
      input:
        vcf = vcf,
        vcf_index = vcf_index,
        contig = contigs[i],
        runtime_docker = runtime_docker
    }
  }

  output {
    Array[File] contig_vcfs = GetContigFromVcf.contig_vcf
    Array[File] contig_vcf_indicies = GetContigFromVcf.contig_vcf_index
  }
}

task GetContigFromVcf {
    input {
      File vcf
      File vcf_index
      String contig
      String runtime_docker

      Float? memory_gib
      Int? disk_gb
      Int? cpus
      Int? preemptible_tries
      Int? max_retries
      Int? boot_disk_gb
    }

    Float disk_size = size([vcf, vcf_index], "GB") * 1.2 + 16
    String output_vcf = "${contig}-${basename(vcf)}"
    String output_vcf_index = "${output_vcf}.tbi"

    runtime {
      memory: "${select_first([memory_gib, 2])} GiB"
      disks: "local-disk ${select_first([disk_gb, ceil(disk_size)])} HDD"
      cpus: select_first([cpus, 1])
      preemptible: select_first([preemptible_tries, 3])
      docker: runtime_docker
      bootDiskSizeGb: select_first([boot_disk_gb, 16])
    }

    command <<<
      bash /opt/task-scripts/SplitVcfByContig/GetContigFromVcf.bash '~{vcf}' '~{contig}' '~{output_vcf}'
    >>>

    output {
      File contig_vcf = "${output_vcf}"
      File contig_vcf_index = "${output_vcf_index}"
    }
}
