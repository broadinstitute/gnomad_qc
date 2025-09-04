version 1.0

workflow hello_workflow {
  input {
    String workspace_bucket
  }

  call helloTask {
    input:
      workspace_bucket = workspace_bucket
  }

  output {
    File hello_file = helloTask.hello_txt
  }
}

task helloTask {
  input {
    String workspace_bucket
  }

  command <<< 
    echo "Hello world" > hello.txt
    pwd
    ls -l
    gsutil cp hello.txt ${workspace_bucket}/hello.txt
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
