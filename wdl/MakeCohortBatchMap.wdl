version 1.0

workflow MakeCohortBatchMap {
  input {
    String sample_set_set_id
    String workspace_namespace
    String workspace_name
    String runtime_docker
  }

  call GatherCohortSamples {
    input:
      workspace_namespace = workspace_namespace,
      workspace_name = workspace_name,
      sample_set_set_id = sample_set_set_id,
      runtime_docker = runtime_docker
  }

  output {
    File batch_map = GatherCohortSamples.batch_map
  }
}

task GatherCohortSamples {
  input {
    String workspace_namespace
    String workspace_name
    String sample_set_set_id
    String runtime_docker

    Float? memory_gib
    Int? boot_disk_gb
    Int? cpus
    Int? disk_gb
    Int? max_retries
    Int? preemptible_tries
  }

  runtime {
    bootDiskSizeGb: select_first([boot_disk_gb, 8])
    cpus: select_first([cpus, 1])
    disks: "local-disk ${select_first([disk_gb, 16])} HDD"
    docker: runtime_docker
    maxRetries: select_first([max_retries, 1])
    memory: "${select_first([memory_gib, 2])} GiB"
    preemptible: select_first([preemptible_tries, 3])
  }

  String sample_batch_map = sample_set_set_id + "-batch_map.tsv"
  command <<<
    /opt/task_scripts/MakeCohortBatchMap/GatherCohortSamples \
      '~{workspace_namespace}' \
      '~{workspace_name}' \
      '~{sample_set_set_id}' \
      '~{sample_batch_map}'
  >>>

  output {
    File batch_map = sample_batch_map
  }
}
