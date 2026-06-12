# =============================================================================
# Step 6: H358 WT vs KEAP1-KO — predatory matrix validation
# =============================================================================
# Place lab data in data/raw/h358/ (see templates/h358_README.txt)
# Usage: Rscript R/06_h358_predatory_matrix.R
# =============================================================================

source("config.R")
source("R/utils.R")

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(ggpubr)
  library(tidyr)
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
})

log_msg("=== Step 6: H358 KEAP1-KO predatory matrix validation ===")

h358_results <- file.path(PATHS$results, "h358")
h358_figures <- file.path(PATHS$figures, "h358")
dir.create(h358_results, recursive = TRUE, showWarnings = FALSE)
dir.create(h358_figures, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(H358$expr_file) || !file.exists(H358$meta_file)) {
  log_msg("H358 data not found. Expected:")
  log_msg("  ", H358$expr_file)
  log_msg("  ", H358$meta_file)
  log_msg("Copy templates/h358_sample_metadata.csv and your expression matrix, then re-run.")
  log_msg("Reference public dataset: ", H358$reference, " (H358 KEAP1 KO clones).")
  quit(save = "no", status = 0)
}

meta <- fread(H358$meta_file)
expr <- fread(H358$expr_file, data.table = FALSE, check.names = FALSE)
gene_col <- colnames(expr)[1]
genes <- expr[[gene_col]]
mat   <- as.matrix(expr[, -1, drop = FALSE])
rownames(mat) <- genes

if ("cell_line" %in% colnames(meta)) {
  meta <- meta[cell_line == H358$cell_line | grepl("358", cell_line, ignore.case = TRUE)]
}
meta <- meta[meta$group %in% c(H358$wt_label, H358$ko_label), ]
common <- intersect(meta$sample_id, colnames(mat))
if (length(common) < 4) stop("Too few H358 samples matched between metadata and expression matrix.")

mat  <- mat[, common, drop = FALSE]
meta <- meta[match(common, meta$sample_id), ]
group <- factor(meta$group, levels = c(H358$wt_label, H358$ko_label))

if (max(mat, na.rm = TRUE) > 50) mat <- log2(mat + 1)

# --- DEA ---
design <- model.matrix(~ group)
fit    <- lmFit(mat, design)
fit    <- eBayes(fit)
dea    <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
dea$gene_symbol <- rownames(dea)
pred_dea <- dea[dea$gene_symbol %in% PREDATORY_GENES, ]
pred_dea <- pred_dea[order(pred_dea$adj.P.Val), ]
fwrite(dea, file.path(h358_results, "dea_h358_all_genes.csv"))
fwrite(pred_dea, file.path(h358_results, "dea_h358_predatory_matrix.csv"))
log_msg("H358 predatory matrix DEA:")
print(pred_dea[, c("gene_symbol", "logFC", "adj.P.Val")])

# --- Correlation heatmap (KO samples if >=3, else all) ---
ko_samp <- meta$sample_id[meta$group == H358$ko_label]
use_samp <- if (length(ko_samp) >= 3) ko_samp else common
sub <- extract_predatory_expr(mat[, use_samp, drop = FALSE])
ct  <- cor_test_matrix(sub)
ct$p_adj <- adjust_cor_pvalues(ct$p)

fwrite(as.data.table(ct$r, keep.rownames = "Gene"),
       file.path(h358_results, "correlation_r_h358.csv"))

col_fun <- colorRamp2(c(-1, 0, 1), c("#2166AC", "#F7F7F7", "#B2182B"))
row_ha  <- rowAnnotation(
  Module = predatory_module_annotation(rownames(sub)),
  col = list(Module = MODULE_COLORS)
)
ht <- Heatmap(
  ct$r, name = "Spearman rho", col = col_fun,
  left_annotation = row_ha,
  cluster_rows = TRUE, cluster_columns = TRUE,
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.text(sprintf("%.2f", ct$r[i, j]), x, y, gp = gpar(fontsize = 10))
    if (i != j && ct$p_adj[i, j] < 0.05) grid.text("*", x, y - unit(2, "mm"))
  },
  row_names_gp = gpar(fontsize = 11, fontface = "italic"),
  column_title = paste0("H358 Predatory Matrix (n=", length(use_samp), " samples)")
)
pdf(file.path(h358_figures, "heatmap_h358_predatory_correlation.pdf"), width = 9, height = 8)
draw(ht, merge_legend = TRUE)
dev.off()
png(file.path(h358_figures, "heatmap_h358_predatory_correlation.png"),
    width = 9, height = 8, units = "in", res = FIG_DPI)
draw(ht, merge_legend = TRUE)
dev.off()

# --- Boxplots ---
plot_df <- meta
for (g in PREDATORY_GENES) {
  if (g %in% rownames(mat)) plot_df[[g]] <- as.numeric(mat[g, plot_df$sample_id])
}
plot_df <- plot_df |>
  pivot_longer(cols = all_of(intersect(PREDATORY_GENES, colnames(plot_df))),
               names_to = "gene", values_to = "expr")

p <- ggplot(plot_df, aes(x = group, y = expr, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.12, size = 1.2) +
  facet_wrap(~ gene, scales = "free_y", ncol = 3) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  scale_fill_manual(values = c("WT" = "#3C5488", "KEAP1_KO" = "#E64B35")) +
  labs(title = "H358: Predatory Matrix Genes", x = NULL, y = "log2(expression + 1)") +
  theme_bw()
ggsave(file.path(h358_figures, "boxplot_h358_predatory_genes.pdf"), p, width = 12, height = 9)
ggsave(file.path(h358_figures, "boxplot_h358_predatory_genes.png"), p, width = 12, height = 9, dpi = FIG_DPI)

log_msg("Step 6 complete. Outputs: ", h358_results, " and ", h358_figures)
