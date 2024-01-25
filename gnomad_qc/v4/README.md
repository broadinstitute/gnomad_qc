# gnomAD v4 production overview:
```mermaid
flowchart LR;
  classDef script_color fill:#2C5D4A
  classDef step_color fill:#2C5D4A
  classDef func_color fill:#2C5D4A
  classDef gnomad_methods_color fill:#2C5D4A
  classDef hail_color fill:#2C5D4A
  classDef resource_color fill:#2C5D4A
  classDef validity_check_color fill:#2C5D4A
  sample_qc --> annotations;
  click sample_qc "./sample_qc/README.md"
  annotations --> variant_qc;
  click annotations "./annotations/README.md"
  variant_qc --> create_release;
  click variant_qc "./variant_qc/README.md"
  click create_release "./create_release/README.md"
```
