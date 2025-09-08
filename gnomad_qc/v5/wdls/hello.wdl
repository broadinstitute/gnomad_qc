version 1.0

workflow hello_workflow {
  input {
    String workspace_bucket
  }

  call helloTask {
    input:
      workspace_bucket_name = workspace_bucket
  }

  output {
    File hello_file = helloTask.hello_txt
  }
}

task helloTask {
    String workspace_bucket_name
  

  command <<< 
    echo "DEBUG: Workspace bucket is '~{workspace_bucket_name}'"
    gsutil ls "~{workspace_bucket_name}" || echo "Bucket path does not exist"
    echo "Hello world" > hello.txt
    gsutil cp hello.txt "~{workspace_bucket_name}/hello.txt"
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
