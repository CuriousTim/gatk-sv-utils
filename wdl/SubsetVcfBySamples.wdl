version 1.0

# Subset one or more VCFs with a list of samples.
#
# If a pedigree is given, the VCF(s) will be subset by families. The families
# to keep can be given by the input `families`. If that input is not given,
# up to `nfamilies` will be randomly selected from the pedigree and used.
# If a pedigree is not given, but a list of samples is given, the VCF(s) will
# be subset using the list of samples. If neither a pedigree nor a list of
# samples is given, `nsamples` samples will be randomly selected from the
# VCF(s) and used.
workflow SubsetVcfBySamples {
  input {
    Array[File] vcfs
    # FamilyID, IndividualID, PaternalID, MaternalID, Sex, Phenotype
    # The pedigree is tab-separated and can have comment lines starting with
    # '#'. That is the only way to include a header.
    File? pedigree
    File? families
    Int nfamilies = 1000
    File? samples
    Int nsamples = 1000

    Boolean? remove_private_sites
    Boolean? keep_af

    # Same length as `vcfs` to give each output VCF a different prefix
    Array[String]? output_prefix_list
    # File with one line per prefix to give each output VCF a different prefix
    File? output_prefix_file
    # Single prefix to use for all output VCFs
    String? output_prefix

    String runtime_docker
  }

  if (defined(pedigree)) {
    call GetSamplesFromFamilies {
      input:
        pedigree = select_first([pedigree]),
        families = families,
        nfamilies = nfamilies,
        runtime_docker = runtime_docker
    }
  }

  Array[String]? prefix_list = if (!defined(output_prefix_list) && defined(output_prefix_file)) then read_lines(select_first([output_prefix_file])) else output_prefix_list
  String static_prefix = select_first([output_prefix, ""])
  scatter (i in range(length(vcfs))) {
    call SubsetVcf {
      input:
        vcf = vcfs[i],
        samples = samples,
        nsamples = nsamples,
        remove_private_sites = remove_private_sites,
        keep_af = keep_af,
        output_prefix = if defined(prefix_list) then select_first([prefix_list])[i] else static_prefix,
        runtime_docker = runtime_docker
    }
  }

  output {
    Array[File] subset_vcfs = SubsetVcf.subset_vcf
    Array[File] subset_vcf_indicies = SubsetVcf.subset_vcf_index
  }
}

task GetSamplesFromFamilies {
  input {
    File pedigree
    File? families
    Int nfamilies
    String runtime_docker

    Int? boot_disk_gb
    Int? cpus
    Int? disk_gb
    Int? max_retries
    Float? memory_gib
    Int? preemptible_tries
  }

  Float disk_size = size(pedigree, "GB") + 16

  runtime {
    bootDiskSizeGb: select_first([boot_disk_gb, 8])
    cpus: select_first([cpus, 1])
    disks: "local-disk ${select_first([disk_gb, ceil(disk_size)])} HDD"
    docker: runtime_docker
    maxRetries: select_first([max_retries, 1])
    memory: "${select_first([memory_gib, 1])} GiB"
    preemptible: select_first([preemptible_tries, 3])
  }

  command <<<
    /opt/task_scripts/SubsetVcfBySamples/GetSamplesFromFamilies \
      ~{if defined(families) then "--families '" + families + "'" else "--nfamilies " + nfamilies} \
      '~{pedigree}' > "samples.list"
  >>>

  output {
    File samples = "samples.list"
  }
}

task SubsetVcf {
  input {
    File vcf
    File? samples
    Int nsamples
    String output_prefix
    Boolean remove_private_sites = true
    Boolean keep_af = true
    String runtime_docker

    Int? boot_disk_gb
    Int? cpus
    Int? disk_gb
    Int? max_retries
    Float? memory_gib
    Int? preemptible_tries
  }

  Float disk_size = size(vcf, "GB") * 2 + size(samples, "GB") + 16

  runtime {
    bootDiskSizeGb: select_first([boot_disk_gb, 8])
    cpus: select_first([cpus, 1])
    disks: "local-disk ${select_first([disk_gb, ceil(disk_size)])} HDD"
    docker: runtime_docker
    maxRetries: select_first([max_retries, 1])
    memory: "${select_first([memory_gib, 1])} GiB"
    preemptible: select_first([preemptible_tries, 3])
  }

  File output_vcf = "${output_prefix}-subset.vcf.gz"
  File output_vcf_index = "${output_vcf}.tbi"

  command <<<
    /opt/task-scripts/SubsetVcfBySamples/SubsetVcf \
      ~{if defined(samples) then "--samples '" + samples + "'" else "--nsamples " + nsamples} \
      ~{if remove_private_sites then "" else "--keep-private-sites"} \
      ~{if keep_af then "" else "--update-af"} \
      '~{vcf}' '~{output_vcf}'
  >>>

  output {
    File subset_vcf = output_vcf
    File subset_vcf_index = output_vcf_index
  }
}
