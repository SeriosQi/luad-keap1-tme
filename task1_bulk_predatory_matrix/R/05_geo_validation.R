# =============================================================================
# Step 5: GEO validation — independent cohorts
# =============================================================================
# Usage: Rscript R/05_geo_validation.R
# Validates predatory matrix upregulation in KEAP1-related GEO datasets
# =============================================================================

source("config.R")
source("R/utils.R")

suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(ggplot2)
  library(data.table)
  library(dplyr)
})

log_msg("=== Step 5: GEO validation ===")

geo_results_dir <- file.path(PATHS$results, "geo_validation")
dir.create(geo_results_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 5A. GSE142694 — KEAP1 CRISPR KO in H1299 (LUAD cell line)
# ---------------------------------------------------------------------------
run_gse142694 <- function() {
  gse_id <- "GSE142694"
  log_msg("Processing ", gse_id, "...")

  cache_rds <- file.path(PATHS$geo, paste0(gse_id, "_expr.rds"))
  if (!file.exists(cache_rds)) {
    gse <- getGEO(gse_id, GSEMatrix = TRUE, getGPL = FALSE)[[1]]
    expr <- exprs(gse)
    pdata <- pData(gse)
    saveRDS(list(expr = expr, pdata = pdata, gse = gse), cache_rds)
  } else {
    cached <- readRDS(cache_rds)
    expr <- cached$expr
    pdata <- cached$pdata
  }

  # Platform gene annotation
  gse <- if (exists("cached")) cached$gse else getGEO(gse_id, GSEMatrix = TRUE)[[1]]
  fdata <- fData(gse)

  # Map probes to symbols
  symbol_col <- grep("Symbol|GENE_SYMBOL|Gene Symbol", colnames(fdata), value = TRUE)[1]
  if (!is.na(symbol_col)) {
    symbols <- fdata[[symbol_col]]
  } else {
    symbols <- map_ensembl_to_symbol(rownames(expr))
  }
  expr_sym <- collapse_by_symbol(expr, symbols)

  # Detect groups from title/characteristics
  group_labels <- pdata$title
  if (is.null(group_labels)) group_labels <- rownames(pdata)

  # GSE142694: typically "KEAP1 KO" vs "Control" or sgKEAP1 vs sgControl
  mut_idx <- grepl("KO|knockout|sgKEAP1|KEAP1[-_ ]?del", group_labels, ignore.case = TRUE)
  wt_idx  <- grepl("control|Ctrl|sgCtrl|NT|wild", group_labels, ignore.case = TRUE)

  if (sum(mut_idx) == 0 || sum(wt_idx) == 0) {
    log_msg("Could not auto-detect groups for ", gse_id,
            ". Inspect pData columns and set manually.")
    fwrite(pdata, file.path(geo_results_dir, paste0(gse_id, "_pData.csv")))
    return(invisible(NULL))
  }

  group <- factor(ifelse(mut_idx, "KEAP1-KO", "Control"), levels = c("Control", "KEAP1-KO"))
  design <- model.matrix(~ group)
  fit <- lmFit(expr_sym[, mut_idx | wt_idx], design)
  fit <- eBayes(fit)
  dea <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
  dea$gene_symbol <- rownames(dea)
  dea <- dea |> relocate(gene_symbol)

  fwrite(dea, file.path(geo_results_dir, paste0(gse_id, "_dea.csv")))

  pred_dea <- dea |> filter(gene_symbol %in% PREDATORY_GENES)
  fwrite(pred_dea, file.path(geo_results_dir, paste0(gse_id, "_predatory_dea.csv")))
  log_msg(gse_id, " predatory matrix DEA:")
  print(pred_dea[, c("gene_symbol", "logFC", "adj.P.Val")])

  # Correlation in KO samples
  ko_samples <- colnames(expr_sym)[mut_idx]
  if (length(ko_samples) >= 3) {
    sub <- extract_predatory_expr(expr_sym[, ko_samples, drop = FALSE])
    ct  <- cor_test_matrix(sub)
    ct$p_adj <- adjust_cor_pvalues(ct$p)
    fwrite(
      as.data.table(ct$r, keep.rownames = "gene"),
      file.path(geo_results_dir, paste0(gse_id, "_correlation_KO.csv"))
    )
  }

  # Boxplot
  plot_samples <- colnames(expr_sym)[mut_idx | wt_idx]
  plot_df <- data.frame(
    sample = plot_samples,
    group  = group[mut_idx | wt_idx]
  )
  for (g in PREDATORY_GENES) {
    if (g %in% rownames(expr_sym)) plot_df[[g]] <- expr_sym[g, plot_samples]
  }
  plot_df <- tidyr::pivot_longer(plot_df, cols = all_of(intersect(PREDATORY_GENES, colnames(plot_df))),
                                   names_to = "gene", values_to = "expr")

  p <- ggplot(plot_df, aes(x = group, y = expr, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.12, size = 1, alpha = 0.6) +
    facet_wrap(~ gene, scales = "free_y") +
    scale_fill_manual(values = c("Control" = "#3C5488", "KEAP1-KO" = "#E64B35")) +
    labs(title = paste0(gse_id, ": Predatory Matrix Genes"), x = NULL, y = "Expression") +
    theme_bw()

  ggsave(file.path(PATHS$figures, paste0(gse_id, "_predatory_boxplot.pdf")), p,
         width = 9, height = 6, dpi = FIG_DPI)

  invisible(dea)
}

# ---------------------------------------------------------------------------
# 5B. GSE68465 — LUAD bulk (co-expression structure validation)
# ---------------------------------------------------------------------------
run_gse68465 <- function() {
  gse_id <- "GSE68465"
  log_msg("Processing ", gse_id, " (co-expression only)...")

  cache_rds <- file.path(PATHS$geo, paste0(gse_id, "_expr.rds"))
  if (!file.exists(cache_rds)) {
    gse <- getGEO(gse_id, GSEMatrix = TRUE, getGPL = FALSE)[[1]]
    expr <- exprs(gse)
    fdata <- fData(gse)
    symbol_col <- grep("Symbol|GENE_SYMBOL", colnames(fdata), value = TRUE)[1]
    symbols <- if (!is.na(symbol_col)) fdata[[symbol_col]] else rownames(expr)
    expr_sym <- collapse_by_symbol(expr, symbols)
    saveRDS(expr_sym, cache_rds)
  } else {
    expr_sym <- readRDS(cache_rds)
  }

  sub <- extract_predatory_expr(expr_sym)
  ct  <- cor_test_matrix(sub)
  ct$p_adj <- adjust_cor_pvalues(ct$p)

  fwrite(
    as.data.table(ct$r, keep.rownames = "gene"),
    file.path(geo_results_dir, paste0(gse_id, "_correlation.csv"))
  )

  # Quick heatmap via pheatmap (lighter weight)
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    pdf(file.path(PATHS$figures, paste0(gse_id, "_correlation_heatmap.pdf")),
        width = 6, height = 5)
    pheatmap::pheatmap(
      ct$r,
      display_numbers = TRUE,
      number_format   = "%.2f",
      color           = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
      breaks          = seq(-1, 1, length.out = 101),
      main            = paste0(gse_id, ": Predatory Matrix Co-expression"),
      fontsize_row    = 11,
      fontsize_col    = 11
    )
    dev.off()
  }

  log_msg(gse_id, " mean pairwise rho: ",
          round(mean(ct$r[upper.tri(ct$r)]), 3))
}

# Run validations
tryCatch(run_gse142694(), error = function(e) {
  log_msg("GSE142694 failed: ", conditionMessage(e))
})
tryCatch(run_gse68465(), error = function(e) {
  log_msg("GSE68465 failed: ", conditionMessage(e))
})

log_msg("Step 5 complete.")
