# =============================================================================
# Step 3: Differential expression — KEAP1-MUT vs KEAP1-WT
# =============================================================================
# Usage: Rscript R/03_differential_expression.R
# Method: DESeq2 on raw counts (recommended for bulk RNA-seq)
# =============================================================================

source("config.R")
source("R/utils.R")

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
  library(tidyr)
  library(data.table)
})

log_msg("=== Step 3: Differential expression analysis ===")

se         <- readRDS(file.path(PATHS$processed, "tcga_luad_counts_se.rds"))
status_df  <- fread(file.path(PATHS$processed, "keap1_status.csv"))
dea_rds    <- file.path(PATHS$results, "dea_keap1_mut_vs_wt.rds")
dea_csv    <- file.path(PATHS$results, "dea_keap1_mut_vs_wt.csv")

# Align counts with status
sample_ids <- normalize_sample_id(colnames(se))
idx <- match(status_df$sample_id, sample_ids)
keep <- !is.na(idx)
se_sub <- se[, idx[keep]]
col_status <- status_df$keap1_status[keep]
col_status <- factor(col_status, levels = c("KEAP1-WT", "KEAP1-MUT"))

log_msg("DESeq2 input: ", ncol(se_sub), " samples (",
        sum(col_status == "KEAP1-MUT"), " MUT / ",
        sum(col_status == "KEAP1-WT"), " WT)")

if (!file.exists(dea_rds)) {
  counts <- assay(se_sub, "unstranded")
  if (is.null(counts)) counts <- assay(se_sub)

  gene_ids <- rowData(se_sub)$gene_id
  symbols  <- map_ensembl_to_symbol(gene_ids)

  # Pre-filter low counts
  keep_genes <- rowSums(counts >= 10) >= ceiling(MIN_SAMPLE_FRACTION * ncol(counts))
  counts <- counts[keep_genes, ]
  symbols  <- symbols[keep_genes]

  dds <- DESeqDataSetFromMatrix(
    countData = round(counts),
    colData   = data.frame(keap1_status = col_status, row.names = colnames(counts)),
    design    = ~ keap1_status
  )
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("keap1_status", "KEAP1-MUT", "KEAP1-WT"))
  res_df <- as.data.frame(res) |>
    tibble::rownames_to_column("ensembl_id") |>
    mutate(gene_symbol = symbols[match(ensembl_id, gene_ids[keep_genes])]) |>
    relocate(gene_symbol, .before = baseMean)

  saveRDS(res_df, dea_rds)
  fwrite(res_df, dea_csv)
  log_msg("DESeq2 complete. Saved: ", dea_csv)
} else {
  res_df <- readRDS(dea_rds)
  log_msg("Loaded cached DEA results.")
}

# --- Focus: predatory matrix genes ---
pred_dea <- res_df |>
  filter(gene_symbol %in% PREDATORY_GENES) |>
  arrange(padj)

fwrite(pred_dea, file.path(PATHS$results, "dea_predatory_matrix_genes.csv"))
log_msg("Predatory matrix DEA results:")
print(pred_dea[, c("gene_symbol", "log2FoldChange", "padj")])

# --- Volcano plot (predatory genes highlighted) ---
res_df$significant <- !is.na(res_df$padj) & res_df$padj < DEA_PADJ &
  abs(res_df$log2FoldChange) >= DEA_LOG2FC
res_df$highlight <- ifelse(res_df$gene_symbol %in% PREDATORY_GENES, "Predatory Matrix", "Other")
res_df$highlight[!(res_df$gene_symbol %in% PREDATORY_GENES) & res_df$significant] <- "Other DEG"

p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = highlight), alpha = 0.5, size = 1.2) +
  scale_color_manual(values = c(
    "Predatory Matrix" = "#E64B35",
    "Other DEG"        = "#4DBBD5",
    "Other"            = "grey70"
  )) +
  geom_vline(xintercept = c(-DEA_LOG2FC, DEA_LOG2FC), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(DEA_PADJ), linetype = "dashed", color = "grey40") +
  ggrepel::geom_text_repel(
    data = subset(res_df, gene_symbol %in% PREDATORY_GENES),
    aes(label = gene_symbol),
    size = 3.5, max.overlaps = 20, color = "#E64B35"
  ) +
  labs(
    title = "TCGA-LUAD: KEAP1-MUT vs KEAP1-WT",
    subtitle = "Predatory matrix genes highlighted",
    x = expression(log[2] ~ "Fold Change (MUT/WT)"),
    y = expression(-log[10] ~ "adjusted" ~ italic(P)),
    color = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(
  file.path(PATHS$figures, "volcano_predatory_matrix.pdf"),
  p_volcano, width = 8, height = 7, dpi = FIG_DPI
)
ggsave(
  file.path(PATHS$figures, "volcano_predatory_matrix.png"),
  p_volcano, width = 8, height = 7, dpi = FIG_DPI
)

# --- Boxplots: individual predatory genes ---
tpm       <- readRDS(file.path(PATHS$processed, "tcga_luad_tpm_matrix.rds"))
plot_df <- as.data.frame(status_df)
for (g in PREDATORY_GENES) {
  if (g %in% rownames(tpm)) {
    plot_df[[g]] <- as.numeric(tpm[g, plot_df$sample_id])
  }
}
plot_genes <- intersect(PREDATORY_GENES, colnames(plot_df))
plot_df <- plot_df |>
  pivot_longer(cols = all_of(plot_genes), names_to = "gene", values_to = "tpm") |>
  mutate(log2_expr = log2(tpm + 1))

p_box <- ggplot(plot_df, aes(x = keap1_status, y = log2_expr, fill = keap1_status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.4) +
  facet_wrap(~ gene, scales = "free_y", ncol = 3) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y.npc = 0.95) +
  scale_fill_manual(values = c("KEAP1-WT" = "#3C5488", "KEAP1-MUT" = "#E64B35")) +
  labs(
    title = "Predatory Matrix Gene Expression",
    subtitle = "TCGA-LUAD: KEAP1-MUT vs KEAP1-WT",
    x = NULL, y = expression(log[2](TPM + 1)), fill = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", strip.background = element_rect(fill = "grey95"))

ggsave(
  file.path(PATHS$figures, "boxplot_predatory_genes.pdf"),
  p_box, width = 10, height = 7, dpi = FIG_DPI
)

# --- Module score comparison ---
module_scores <- compute_module_score(tpm)
score_df <- status_df |>
  mutate(predatory_score = module_scores[sample_id])

p_score <- ggplot(score_df, aes(x = keap1_status, y = predatory_score, fill = keap1_status)) +
  geom_violin(trim = FALSE, alpha = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  scale_fill_manual(values = c("KEAP1-WT" = "#3C5488", "KEAP1-MUT" = "#E64B35")) +
  labs(
    title = "Predatory Matrix Module Score",
    subtitle = "Mean z-score of SLC7A11, GGT1, SLC1A5, ABCC1/2/3",
    x = NULL, y = "Module Score", fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

ggsave(
  file.path(PATHS$figures, "module_score_keap1_groups.pdf"),
  p_score, width = 5, height = 5, dpi = FIG_DPI
)

fwrite(score_df, file.path(PATHS$results, "predatory_module_scores.csv"))
log_msg("Step 3 complete.")
