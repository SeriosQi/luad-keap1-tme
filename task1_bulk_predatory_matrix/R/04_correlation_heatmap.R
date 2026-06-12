# =============================================================================
# Step 4: Predatory matrix co-expression network & correlation heatmap
# =============================================================================
# Usage: Rscript R/04_correlation_heatmap.R
# Output: Publication-quality ComplexHeatmap (PDF/PNG)
# =============================================================================

source("config.R")
source("R/utils.R")

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(data.table)
  library(dplyr)
})

log_msg("=== Step 4: Predatory matrix correlation analysis ===")

tpm       <- readRDS(file.path(PATHS$processed, "tcga_luad_tpm_matrix.rds"))
status_df <- fread(file.path(PATHS$processed, "keap1_status.csv"))

# Log2-transform TPM for correlation stability
log_expr <- log2(tpm + 1)

# Split by KEAP1 status
mut_samples <- status_df$sample_id[status_df$keap1_status == "KEAP1-MUT"]
wt_samples  <- status_df$sample_id[status_df$keap1_status == "KEAP1-WT"]
all_samples <- status_df$sample_id

groups <- list(
  "All LUAD"   = all_samples,
  "KEAP1-MUT"  = mut_samples,
  "KEAP1-WT"   = wt_samples
)

cor_results <- list()

for (grp_name in names(groups)) {
  samp <- intersect(groups[[grp_name]], colnames(log_expr))
  sub  <- extract_predatory_expr(log_expr[, samp, drop = FALSE])
  ct   <- cor_test_matrix(sub, method = COR_METHOD)
  ct$p_adj <- adjust_cor_pvalues(ct$p)
  cor_results[[grp_name]] <- ct

  # Save numeric tables
  fwrite(
    as.data.table(ct$r, keep.rownames = "gene") |> setnames("gene", "Gene"),
    file.path(PATHS$results, paste0("correlation_r_", gsub(" ", "_", grp_name), ".csv"))
  )
  fwrite(
    as.data.table(ct$p_adj, keep.rownames = "gene") |> setnames("gene", "Gene"),
    file.path(PATHS$results, paste0("correlation_padj_", gsub(" ", "_", grp_name), ".csv"))
  )
}

# --- Significance annotation function for heatmap cells ---
sig_mark <- function(r_mat, p_mat, alpha = 0.05) {
  marks <- matrix("", nrow = nrow(r_mat), ncol = ncol(r_mat),
                  dimnames = dimnames(r_mat))
  idx <- upper.tri(r_mat) & p_mat < alpha
  marks[idx] <- ifelse(abs(r_mat[idx]) >= 0.5, "**",
                       ifelse(abs(r_mat[idx]) >= 0.3, "*", ""))
  marks
}

# --- Draw multi-panel correlation heatmap ---
draw_correlation_heatmap <- function(cor_list, out_prefix) {
  col_fun <- colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B"))

  ht_list <- NULL
  for (i in seq_along(cor_list)) {
    grp  <- names(cor_list)[i]
    ct   <- cor_list[[grp]]
    r    <- ct$r
    p    <- ct$p_adj
    n_samp <- length(groups[[grp]])
    marks <- sig_mark(r, p)

    cell_fun <- function(j, i, x, y, w, h, fill) {
      if (i > j) {
        grid.text(sprintf("%.2f", r[i, j]), x, y, gp = gpar(fontsize = 9))
        if (marks[i, j] != "") {
          grid.text(marks[i, j], x, y + unit(2, "mm"),
                    gp = gpar(fontsize = 8, col = "black"))
        }
      }
    }

    ht <- Heatmap(
      r,
      name              = paste0("rho (", grp, ")"),
      col               = col_fun,
      cluster_rows      = TRUE,
      cluster_columns   = TRUE,
      show_row_names    = TRUE,
      show_column_names = TRUE,
      row_names_gp      = gpar(fontsize = 11, fontface = "italic"),
      column_names_gp   = gpar(fontsize = 11, fontface = "italic"),
      cell_fun          = cell_fun,
      column_title      = paste0(grp, "\n(Spearman rho, n=", n_samp, " samples)"),
      column_title_gp   = gpar(fontsize = 12, fontface = "bold"),
      width             = unit(4, "cm"),
      height            = unit(4, "cm"),
      heatmap_legend_param = list(
        title = "Correlation",
        at    = c(-1, -0.5, 0, 0.5, 1),
        legend_height = unit(3, "cm")
      )
    )
    ht_list <- ht_list + ht
  }

  pdf(paste0(out_prefix, ".pdf"), width = 14, height = 5)
  draw(ht_list, merge_legend = TRUE,
       heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()

  png(paste0(out_prefix, ".png"), width = 14, height = 5, units = "in", res = FIG_DPI)
  draw(ht_list, merge_legend = TRUE,
       heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
}

draw_correlation_heatmap(
  cor_results,
  file.path(PATHS$figures, "heatmap_predatory_correlation_by_group")
)

# --- Enhanced single-panel heatmap: KEAP1-MUT group (primary hypothesis) ---
if (length(mut_samples) >= 3) {
  ct_mut <- cor_results[["KEAP1-MUT"]]
  r_mut  <- ct_mut$r
  p_mut  <- ct_mut$p_adj

  # Combine r and significance stars into display matrix
  display_mat <- matrix(
    sprintf("%.2f%s", r_mut, ifelse(upper.tri(r_mut) & p_mut < 0.05, "*", "")),
    nrow = nrow(r_mut), dimnames = dimnames(r_mut)
  )

  col_fun <- colorRamp2(c(-1, 0, 1), c("#2166AC", "#F7F7F7", "#B2182B"))

  ht_mut <- Heatmap(
    r_mut,
    name = "Spearman rho",
    col  = col_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    cell_fun = function(j, i, x, y, w, h, fill) {
      grid.text(sprintf("%.2f", r_mut[i, j]), x, y, gp = gpar(fontsize = 11))
      if (i != j && p_mut[i, j] < 0.05) {
        grid.text("*", x, y - unit(2.5, "mm"), gp = gpar(fontsize = 10))
      }
    },
    row_names_gp    = gpar(fontsize = 12, fontface = "italic"),
    column_names_gp = gpar(fontsize = 12, fontface = "italic"),
    column_title = paste0(
      "Predatory Matrix Co-expression\nKEAP1-MUT LUAD (n=", length(mut_samples), ")"
    ),
    column_title_gp = gpar(fontsize = 13, fontface = "bold"),
    heatmap_legend_param = list(
      title = expression(Spearman ~ rho),
      at = c(-1, -0.5, 0, 0.5, 1)
    )
  )

  pdf(file.path(PATHS$figures, "heatmap_predatory_correlation_KEAP1_MUT.pdf"),
      width = 7, height = 6)
  draw(ht_mut)
  dev.off()

  png(file.path(PATHS$figures, "heatmap_predatory_correlation_KEAP1_MUT.png"),
      width = 7, height = 6, units = "in", res = FIG_DPI)
  draw(ht_mut)
  dev.off()
}

# --- Compare mean pairwise correlation across groups (summary statistic) ---
mean_rho <- sapply(cor_results, function(ct) {
  r <- ct$r
  mean(r[upper.tri(r)], na.rm = TRUE)
})
mean_rho_df <- data.frame(group = names(mean_rho), mean_pairwise_rho = mean_rho)
fwrite(mean_rho_df, file.path(PATHS$results, "mean_pairwise_correlation_summary.csv"))

log_msg("Mean pairwise Spearman rho by group:")
print(mean_rho_df)

# --- Fisher z-test: correlation difference MUT vs WT (SLC7A11-GGT1 key pair) ---
if (length(mut_samples) >= 5 && length(wt_samples) >= 5) {
  key_pairs <- list(
    c("SLC7A11", "GGT1"),
    c("GGT1", "SLC1A5"),
    c("SLC7A11", "SLC1A5")
  )
  pair_tests <- lapply(key_pairs, function(pair) {
    g1 <- pair[1]; g2 <- pair[2]
    x_mut <- log_expr[g1, mut_samples]; y_mut <- log_expr[g2, mut_samples]
    x_wt  <- log_expr[g1, wt_samples];  y_wt  <- log_expr[g2, wt_samples]
    if (!all(c(g1, g2) %in% rownames(log_expr))) return(NULL)

    r_mut <- cor(x_mut, y_mut, method = "spearman")
    r_wt  <- cor(x_wt,  y_wt,  method = "spearman")
    # Fisher z transformation difference
    z_diff <- atanh(r_mut) - atanh(r_wt)
    se <- sqrt(1 / (length(mut_samples) - 3) + 1 / (length(wt_samples) - 3))
    z_stat <- z_diff / se
    p_val  <- 2 * pnorm(-abs(z_stat))
    data.frame(gene1 = g1, gene2 = g2, rho_mut = r_mut, rho_wt = r_wt,
               z_diff = z_diff, p_value = p_val)
  })
  pair_tests <- bind_rows(pair_tests)
  fwrite(pair_tests, file.path(PATHS$results, "key_pair_correlation_mut_vs_wt.csv"))
  log_msg("Key pair correlation comparison (MUT vs WT):")
  print(pair_tests)
}

log_msg("Step 4 complete.")
