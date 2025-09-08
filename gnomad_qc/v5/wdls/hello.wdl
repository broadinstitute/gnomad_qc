version 1.0

workflow hello_workflow {
  input {
    String workspace_bucket
  }

  call helloTask

  output {
    File hello_file = helloTask.hello_txt
  }
}

task helloTask {
  command <<< 
    echo "DEBUG: Workspace bucket is '~{workspace_bucket}'"
    gsutil ls "~{workspace_bucket}" || echo "Bucket path does not exist"
    echo "Hello world" > hello.txt
    gsutil cp hello.txt "~{workspace_bucket}/hello.txt"
>>>

  output {
    File hello_txt = "hello.txt"
  }

  runtime {
    docker: "google/cloud-sdk:latest"
    cpu: 1
    memory: "2 GB"
    bootDiskSizeGb: 20
    stdout: "hello.stdout"
    stderr: "hello.stderr"
  }
}
