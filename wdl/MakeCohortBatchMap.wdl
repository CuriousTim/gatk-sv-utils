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
    set -o errexit
    set -o nounset
    set -o pipefail

    workspace_namespace='~{workspace_namespace}'
    workspace_name='~{workspace_name}'
    sample_set_set_id='~{sample_set_set_id}'
    output='~{sample_batch_map}'

    sample_sets="$(mktemp -p "${PWD}" tmp_XXXXXXXXX)"
    samples="$(mktemp -p "${PWD}" tmp_XXXXXXXXX)"
    response="$(mktemp -p "${PWD}" tmp_XXXXXXXXX)"
    trap 'rm -f "${sample_sets}" "${samples}" "${response}"' EXIT

    curl --oauth2-bearer "$(gcloud auth application-default print-access-token)" \
      --request GET --header 'accept: */*' \
      --silent --show-error \
      "https://api.firecloud.org/api/workspaces/${workspace_namespace}/${workspace_name}/entities/sample_set_set/${sample_set_set_id}" \
      | jq --raw-output '.attributes.sample_sets.items[] | .entityName' > "${sample_sets}"
    if [[ ! -s "${sample_sets}" ]]; then
      printf 'no sample sets found\n' >&2
      exit 1
    fi

    read -r id < "${sample_sets}"
    if [[ "${id}" = 'null' ]]; then
      printf 'no sample sets found\n' >&2
      exit 1
    fi

    declare -i i=1
    while true; do
      curl --oauth2-bearer "$(gcloud auth application-default print-access-token)" \
        --request GET --header 'accept: application/json' \
        --url-query "page=${i}" \
        --url-query "pageSize=10" \
        --url-query "fields=samples" \
        --silent --show-error \
        "https://api.firecloud.org/api/workspaces/${workspace_namespace}/${workspace_name}/entityQuery/sample_set" > "${response}"
      declare -i max_pages
      max_pages=$(jq .resultMetadata.filteredPageCount "${response}")
      jq --raw-output \
        '.results[] | .name as $batch | .attributes.samples.items[] | [$batch, .entityName] | @tsv' "${response}"  >> "${samples}"
      if (( i >= max_pages )); then
        break
      fi
      i=$(( i + 1 ))
    done

    awk -F'\t' 'NR==FNR{a[$1]} NR>FNR && ($1 in a)' "${sample_sets}" "${samples}" > "${output}"
  >>>

  output {
    File batch_map = sample_batch_map
  }
}
