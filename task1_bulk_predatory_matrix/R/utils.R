# =============================================================================
# Utility functions for Task 1 pipeline
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
})

#' Normalize TCGA barcode to 15-char patient ID
normalize_patient_id <- function(x) {
  x <- gsub("\\.", "-", as.character(x))
  substr(x, 1, 15)
}

#' Normalize TCGA barcode to 16-char sample ID (tumor/normal suffix)
normalize_sample_id <- function(x) {
  x <- gsub("\\.", "-", as.character(x))
  substr(x, 1, 16)
}

#' Map Ensembl IDs to HGNC symbols (handles version suffix)
map_ensembl_to_symbol <- function(ensembl_ids) {
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
      !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("Install org.Hs.eg.db: BiocManager::install('org.Hs.eg.db')")
  }
  clean <- sub("\\..*", "", ensembl_ids)
  sym <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = clean,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  unname(sym[match(clean, names(sym))])
}

#' Collapse duplicate gene symbols by max expression
collapse_by_symbol <- function(expr_mat, gene_symbols) {
  dt <- as.data.table(expr_mat, keep.rownames = FALSE)
  dt[, symbol := gene_symbols]
  sample_cols <- colnames(expr_mat)
  dt <- dt[, lapply(.SD, max, na.rm = TRUE), by = symbol, .SDcols = sample_cols]
  mat <- as.matrix(dt[, ..sample_cols])
  rownames(mat) <- dt$symbol
  mat
}

#' Extract predatory matrix genes from expression matrix
extract_predatory_expr <- function(expr_mat, genes = PREDATORY_GENES) {
  present <- intersect(genes, rownames(expr_mat))
  missing <- setdiff(genes, present)
  if (length(missing) > 0) {
    message("Missing genes in expression matrix: ", paste(missing, collapse = ", "))
  }
  expr_mat[present, , drop = FALSE]
}

#' Compute predatory matrix module score (mean z-score)
compute_module_score <- function(expr_mat, genes = PREDATORY_GENES) {
  sub <- extract_predatory_expr(expr_mat, genes)
  z <- t(scale(t(sub)))
  colMeans(z, na.rm = TRUE)
}

#' Pairwise correlation with p-values
cor_test_matrix <- function(mat, method = COR_METHOD) {
  genes <- rownames(mat)
  n <- length(genes)
  r_mat <- matrix(NA, n, n, dimnames = list(genes, genes))
  p_mat <- matrix(NA, n, n, dimnames = list(genes, genes))

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) {
        r_mat[i, j] <- 1
        p_mat[i, j] <- 0
      } else if (i < j) {
        ct <- stats::cor.test(mat[i, ], mat[j, ], method = method, exact = FALSE)
        r_mat[i, j] <- r_mat[j, i] <- ct$estimate
        p_mat[i, j] <- p_mat[j, i] <- ct$p.value
      }
    }
  }
  list(r = r_mat, p = p_mat)
}

#' BH-adjust p-values in correlation matrix (upper triangle)
adjust_cor_pvalues <- function(p_mat) {
  idx <- upper.tri(p_mat)
  p_adj <- p_mat
  p_adj[idx] <- p.adjust(p_mat[idx], method = COR_PADJ_METHOD)
  p_adj[lower.tri(p_adj)] <- t(p_adj)[lower.tri(p_adj)]
  p_adj
}

#' Save session info for reproducibility
save_session_info <- function(out_file) {
  writeLines(capture.output(sessionInfo()), out_file)
}

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)
  message(msg)
  invisible(msg)
}

#' Row annotation vector for predatory matrix heatmaps
predatory_module_annotation <- function(genes) {
  mod <- PREDATORY_GENE_MODULE[genes]
  mod[is.na(mod)] <- "Other"
  factor(mod, levels = names(MODULE_COLORS))
}
