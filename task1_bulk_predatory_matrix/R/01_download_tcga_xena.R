# =============================================================================
# Step 1 (Scheme A): Download TCGA-LUAD from GDC Xena Hub (2 files)
# =============================================================================
# Usage: Rscript R/01_download_tcga_xena.R
# Output: tcga_luad_tpm_matrix.rds, xena_mutation_dt.rds, data_source.txt
# =============================================================================

source("config.R")
source("R/utils.R")

log_msg("=== Step 1 (Scheme A): GDC Xena Hub download ===")

xena_dir <- file.path(PATHS$raw, "xena")
dir.create(xena_dir, recursive = TRUE, showWarnings = FALSE)

expr_gz  <- file.path(xena_dir, XENA_FILES$expression$dest)
mut_gz   <- file.path(xena_dir, XENA_FILES$mutation$dest)
expr_rds <- file.path(PATHS$processed, "tcga_luad_tpm_matrix.rds")
mut_rds  <- file.path(PATHS$processed, "xena_mutation_dt.rds")
src_flag <- file.path(PATHS$processed, "data_source.txt")

download_xena <- function(url, dest) {
  if (file.exists(dest) && file.info(dest)$size > 1e6) {
    log_msg("Cached: ", basename(dest), " (", round(file.info(dest)$size / 1e6, 1), " MB)")
    return(invisible(dest))
  }
  log_msg("Downloading: ", basename(dest))
  # wget -c supports resume on interrupted transfers
  cmd <- sprintf('wget -c --tries=10 --timeout=60 -O "%s" "%s"', dest, url)
  ret <- system(cmd, intern = FALSE)
  if (ret != 0 || !file.exists(dest) || file.info(dest)$size < 1e3) {
    stop("Download failed: ", basename(dest))
  }
  log_msg("Downloaded: ", basename(dest), " (", round(file.info(dest)$size / 1e6, 1), " MB)")
  invisible(dest)
}

download_xena(XENA_FILES$expression$url, expr_gz)
download_xena(XENA_FILES$mutation$url, mut_gz)

# --- Parse expression matrix (genes x samples, FPKM-UQ) ---
if (!file.exists(expr_rds)) {
  log_msg("Parsing expression matrix...")
  dt <- data.table::fread(expr_gz, data.table = FALSE, check.names = FALSE)
  gene_col <- colnames(dt)[1]
  genes <- dt[[gene_col]]
  expr  <- as.matrix(dt[, -1, drop = FALSE])
  rownames(expr) <- genes

  # Map Ensembl → symbol if needed
  if (any(grepl("^ENSG", genes))) {
    symbols <- map_ensembl_to_symbol(genes)
    expr <- collapse_by_symbol(expr, symbols)
  } else {
    rownames(expr) <- genes
  }

  # Keep primary tumor samples (barcode type 01)
  sample_ids <- normalize_sample_id(colnames(expr))
  tumor_idx  <- substr(sample_ids, 14, 15) == "01"
  expr       <- expr[, tumor_idx, drop = FALSE]
  colnames(expr) <- normalize_sample_id(colnames(expr))

  saveRDS(expr, expr_rds)
  log_msg("Expression matrix: ", nrow(expr), " genes x ", ncol(expr), " tumor samples")
} else {
  log_msg("Expression matrix cached: ", expr_rds)
}

# --- Parse mutation table ---
if (!file.exists(mut_rds)) {
  log_msg("Parsing somatic mutation matrix...")
  mut_dt <- data.table::fread(mut_gz, data.table = TRUE, check.names = FALSE)
  saveRDS(mut_dt, mut_rds)
  log_msg("Mutation table: ", nrow(mut_dt), " rows x ", ncol(mut_dt), " cols")
} else {
  log_msg("Mutation table cached: ", mut_rds)
}

writeLines("xena_gdc_hub", src_flag)
log_msg("Data source flag: ", src_flag)
log_msg("Step 1 (Scheme A) complete.")
