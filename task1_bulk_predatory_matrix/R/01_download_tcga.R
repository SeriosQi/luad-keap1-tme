# =============================================================================
# Step 1: Download TCGA-LUAD RNA-seq counts and MAF (mutation) data
# =============================================================================
# Usage: Rscript R/01_download_tcga.R
# Requires: TCGAbiolinks, SummarizedExperiment
# =============================================================================

source("config.R")
source("R/utils.R")

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
})

log_msg("=== Step 1: TCGA-LUAD data acquisition ===")

counts_rds <- file.path(PATHS$processed, "tcga_luad_counts_se.rds")
maf_rds    <- file.path(PATHS$processed, "tcga_luad_maf.rds")

# --- 1A. RNA-seq raw counts (STAR - Counts) ---
if (!file.exists(counts_rds)) {
  log_msg("Querying GDC for TCGA-LUAD gene expression counts...")
  query_exp <- GDCquery(
    project         = TCGA_PROJECT,
    data.category   = "Transcriptome Profiling",
    data.type       = "Gene Expression Quantification",
    workflow.type   = "STAR - Counts"
  )
  owd <- getwd()
  on.exit(setwd(owd), add = TRUE)
  setwd(PATHS$raw)
  GDCdownload(query_exp, directory = ".")
  setwd(owd)
  se <- GDCprepare(query_exp, directory = PATHS$raw, summarizedExperiment = TRUE)

  # Keep primary tumor samples only (barcode type 01)
  sample_types <- substr(colnames(se), 14, 15)
  se <- se[, sample_types == "01"]

  saveRDS(se, counts_rds)
  log_msg("Saved SummarizedExperiment: ", counts_rds,
          " (", ncol(se), " tumor samples)")
} else {
  se <- readRDS(counts_rds)
  log_msg("Loaded cached counts: ", ncol(se), " samples")
}

# --- 1B. Somatic mutation MAF ---
if (!file.exists(maf_rds)) {
  log_msg("Querying GDC for TCGA-LUAD somatic mutations (MAF)...")
  query_maf <- GDCquery(
    project         = TCGA_PROJECT,
    data.category   = "Simple Nucleotide Variation",
    data.type       = "Masked Somatic Mutation",
    workflow.type   = "Aliquot Ensemble Somatic Variant Merging and Masking"
  )
  owd <- getwd()
  on.exit(setwd(owd), add = TRUE)
  setwd(PATHS$raw)
  GDCdownload(query_maf, directory = ".")
  setwd(owd)
  maf <- GDCprepare(query_maf, directory = PATHS$raw)
  saveRDS(maf, maf_rds)
  log_msg("Saved MAF: ", nrow(maf), " variants")
} else {
  maf <- readRDS(maf_rds)
  log_msg("Loaded cached MAF: ", nrow(maf), " variants")
}

# --- 1C. Build TPM matrix for visualization / correlation ---
tpm_rds <- file.path(PATHS$processed, "tcga_luad_tpm_matrix.rds")
if (!file.exists(tpm_rds)) {
  log_msg("Building TPM matrix from SummarizedExperiment...")
  # GDC SummarizedExperiment: unstranded counts + tpm_unstrand in assays
  if ("tpm_unstrand" %in% assayNames(se)) {
    tpm <- assay(se, "tpm_unstrand")
  } else if ("fpkm_unstrand" %in% assayNames(se)) {
    tpm <- assay(se, "fpkm_unstrand")
    log_msg("TPM not found; using FPKM as proxy.")
  } else {
    stop("No tpm_unstrand or fpkm_unstrand assay found in SE object.")
  }

  gene_ids <- rowData(se)$gene_id
  if (is.null(gene_ids)) gene_ids <- rownames(se)
  symbols <- map_ensembl_to_symbol(gene_ids)
  tpm_sym <- collapse_by_symbol(tpm, symbols)
  colnames(tpm_sym) <- normalize_sample_id(colnames(tpm_sym))

  saveRDS(tpm_sym, tpm_rds)
  log_msg("TPM matrix: ", nrow(tpm_sym), " genes x ", ncol(tpm_sym), " samples")
} else {
  log_msg("TPM matrix already exists: ", tpm_rds)
}

log_msg("Step 1 complete.")
