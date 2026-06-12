H358 local RNA-seq input format
================================

Copy your lab data into:
  data/raw/h358/expression_matrix.csv
  data/raw/h358/sample_metadata.csv

expression_matrix.csv
  - Rows: gene symbols (HGNC)
  - Columns: sample IDs matching metadata
  - Values: raw counts OR FPKM/TPM (script auto-detects log scale)

sample_metadata.csv (minimum columns)
  - sample_id: column names in expression matrix
  - cell_line: H358
  - genotype: WT or KEAP1_KO
  - group: WT or KEAP1_KO  (used for DEA: KO vs WT)

Use templates/h358_sample_metadata.csv as a template.
Only H358 rows are analyzed if multiple cell lines are present.

Public reference (optional): ArrayExpress E-MTAB-9724 (H358/H292 KEAP1 KO).
Place processed H358 counts in the paths above if using your own export from that study.

Run: Rscript R/06_h358_predatory_matrix.R
