# gnomAD v4 production overview:
```mermaid
flowchart LR;
  classDef script_color fill:#2C5D4A,color:#000000
  classDef step_color fill:#98BFFF,color:#000000
  classDef func_color fill:#E0D5E6,color:#000000
  classDef gnomad_methods_color fill:#F9F2CE,color:#000000
  classDef hail_color fill:#FAEACE,color:#000000
  classDef resource_color fill:#D6D6D6,color:#000000
  classDef validity_check_color fill:#2C5D4A,color:#000000
  sample_qc --> annotations;
  click sample_qc "./sample_qc/README.md"
  annotations --> variant_qc;
  click annotations "./annotations/README.md"
  variant_qc --> create_release;
  click variant_qc "./variant_qc/README.md"
  click create_release "./create_release/README.md"
```
