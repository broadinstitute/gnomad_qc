version 1.0

workflow describe_ht_workflow {
  input {
    File input_ht        # your Hail Table
    String output_name   # base name, e.g. "my_ht_description"
    String workspace_bucket  # e.g. "gs://fc-secure-1234-5678"
  }

  call describeHT {
    input:
      ht = input_ht,
      output_name = output_name,
      workspace_bucket = workspace_bucket
  }

  output {
    File description = describeHT.final_output
  }
}

task describeHT {
  input {
    File ht
    String output_name
    String workspace_bucket
  }

  command <<<
    python <<CODE
import hail as hl

ht = hl.read_table("${ht}")
desc = ht.describe()

outfile = f"{output_name}.txt"
with open(outfile, "w") as f:
    f.write(f"Output name: {output_name}\n")
    f.write("="*40 + "\n")
    f.write(desc)
CODE

    # copy to the workspace bucket so it's easy to find
    gsutil cp ${output_name}.txt ${workspace_bucket}/${output_name}.txt
  >>>

  output {
    File final_output = "${workspace_bucket}/${output_name}.txt"
  }

  runtime {
    docker: "hailgenetics/hail:latest"
  }
}
