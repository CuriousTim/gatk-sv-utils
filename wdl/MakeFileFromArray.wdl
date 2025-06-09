version 1.0

# Make a File from an Array, with one line per Array element
workflow MakeFileFromArray {
  input {
    Array[String] arr
    String runtime_docker
  }

  call MakeFileFromArrayTask {
    input:
      arr = arr,
      runtime_docker = runtime_docker
  }

  output {
    File file = MakeFileFromArrayTask.file
  }
}

task MakeFileFromArrayTask {
  input {
    Array[String] arr
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
    memory: "${select_first([memory_gib, 1])} GiB"
    preemptible: select_first([preemptible_tries, 3])
  }

  command <<<
    /opt/task_scripts/MakeFileFromArray/MakeFileFromArrayTask \
      '~{write_lines(arr)}' arrayfile.list
  >>>

  output {
    File file = "arrayfile.list"
  }
}
