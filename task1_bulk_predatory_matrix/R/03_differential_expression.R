# =============================================================================
# Step 3: Differential expression — KEAP1-MUT vs KEAP1-WT
# Scheme A: limma on log2(FPKM-UQ + 1)  |  GDC: DESeq2 on raw counts
# =============================================================================

source("config.R")
source("R/utils.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
  library(tidyr)
  library(data.table)
})

log_msg("=== Step 3: Differential expression analysis ===")

`%||%` <- function(a, b) if (!is.null(a)) a else b

status_df <- fread(file.path(PATHS$processed, "keap1_status.csv"))
dea_rds   <- file.path(PATHS$results, "dea_keap1_mut_vs_wt.rds")
dea_csv   <- file.path(PATHS$results, "dea_keap1_mut_vs_wt.csv")
src_flag  <- file.path(PATHS$processed, "data_source.txt")
data_src  <- if (file.exists(src_flag)) readLines(src_flag, n = 1) else "gdc"

run_limma_dea <- function(expr_mat, status_df) {
  suppressPackageStartupMessages(library(limma))

  common <- intersect(status_df$sample_id, colnames(expr_mat))
  expr   <- expr_mat[, common, drop = FALSE]
  st     <- status_df[match(common, status_df$sample_id), ]
  group  <- factor(st$keap1_status, levels = c("KEAP1-WT", "KEAP1-MUT"))

  log_msg("limma input: ", ncol(expr), " samples (",
          sum(group == "KEAP1-MUT"), " MUT / ",
          sum(group == "KEAP1-WT"), " WT)")

  # FPKM-UQ: log2 transform if not already log-scaled
  if (max(expr, na.rm = TRUE) > 50) {
    log_expr <- log2(expr + 1)
  } else {
    log_expr <- expr
    log_msg("Expression appears log-scaled; using as-is.")
  }

  # Clean row names for limma
  valid <- !is.na(rownames(log_expr)) & rownames(log_expr) != ""
  log_expr <- log_expr[valid, , drop = FALSE]
  log_expr <- log_expr[!duplicated(rownames(log_expr)), , drop = FALSE]

  keep_genes <- rowSums(log_expr > 0) >= ceiling(MIN_SAMPLE_FRACTION * ncol(log_expr))
  log_expr   <- log_expr[keep_genes, , drop = FALSE]

  design <- model.matrix(~ group)
  fit    <- lmFit(log_expr, design)
  fit    <- eBayes(fit)
  tt     <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  tt$gene_symbol <- rownames(tt)
  tt <- tt |>
    tibble::rownames_to_column("ensembl_id") |>
    mutate(
      log2FoldChange = logFC,
      padj = adj.P.Val,
      baseMean = AveExpr
    ) |>
    relocate(gene_symbol, .before = baseMean)

  list(res_df = tt, log_expr = log_expr, method = "limma_log2_FPKM-UQ")
}

run_deseq2_dea <- function(se, status_df) {
  suppressPackageStartupMessages(library(DESeq2))

  sample_ids <- normalize_sample_id(colnames(se))
  idx <- match(status_df$sample_id, sample_ids)
  keep <- !is.na(idx)
  se_sub <- se[, idx[keep]]
  col_status <- factor(status_df$keap1_status[keep], levels = c("KEAP1-WT", "KEAP1-MUT"))

  log_msg("DESeq2 input: ", ncol(se_sub), " samples (",
          sum(col_status == "KEAP1-MUT"), " MUT / ",
          sum(col_status == "KEAP1-WT"), " WT)")

  counts <- assay(se_sub, "unstranded")
  if (is.null(counts)) counts <- assay(se_sub)
  gene_ids <- rowData(se_sub)$gene_id
  symbols  <- map_ensembl_to_symbol(gene_ids)

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

  list(res_df = res_df, method = "DESeq2_GDC_counts")
}

if (!file.exists(dea_rds)) {
  if (data_src == "xena_gdc_hub") {
    expr_mat <- readRDS(file.path(PATHS$processed, "tcga_luad_tpm_matrix.rds"))
    out <- run_limma_dea(expr_mat, status_df)
    res_df <- out$res_df
    dea_method <- out$method
  } else {
    se <- readRDS(file.path(PATHS$processed, "tcga_luad_counts_se.rds"))
    out <- run_deseq2_dea(se, status_df)
    res_df <- out$res_df
    dea_method <- out$method
  }
  attr(res_df, "method") <- dea_method
  saveRDS(res_df, dea_rds)
  fwrite(res_df, dea_csv)
  log_msg("DEA complete (", dea_method, "). Saved: ", dea_csv)
} else {
  res_df <- readRDS(dea_rds)
  dea_method <- attr(res_df, "method") %||% "cached"
  log_msg("Loaded cached DEA results (", dea_method, ").")
}

# --- Predatory matrix genes ---
pred_dea <- res_df |>
  filter(gene_symbol %in% PREDATORY_GENES) |>
  arrange(padj)

fwrite(pred_dea, file.path(PATHS$results, "dea_predatory_matrix_genes.csv"))
log_msg("Predatory matrix DEA results:")
print(pred_dea[, c("gene_symbol", "log2FoldChange", "padj")])

# --- Volcano plot ---
res_df$significant <- !is.na(res_df$padj) & res_df$padj < DEA_PADJ &
  abs(res_df$log2FoldChange) >= DEA_LOG2FC
res_df$highlight <- ifelse(res_df$gene_symbol %in% PREDATORY_GENES, "Predatory Matrix", "Other")
res_df$highlight[!(res_df$gene_symbol %in% PREDATORY_GENES) & res_df$significant] <- "Other DEG"

dea_label <- if (grepl("limma", dea_method)) "limma, log2(FPKM-UQ+1)" else "DESeq2, GDC counts"

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
    subtitle = paste0("Predatory matrix highlighted | ", dea_label),
    x = expression(log[2] ~ "Fold Change (MUT/WT)"),
    y = expression(-log[10] ~ "adjusted" ~ italic(P)),
    color = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(PATHS$figures, "volcano_predatory_matrix.pdf"),
       p_volcano, width = 8, height = 7, dpi = FIG_DPI)
ggsave(file.path(PATHS$figures, "volcano_predatory_matrix.png"),
       p_volcano, width = 8, height = 7, dpi = FIG_DPI)

# --- Boxplots ---
tpm <- readRDS(file.path(PATHS$processed, "tcga_luad_tpm_matrix.rds"))
plot_df <- as.data.frame(status_df)
for (g in PREDATORY_GENES) {
  if (g %in% rownames(tpm)) {
    plot_df[[g]] <- as.numeric(tpm[g, plot_df$sample_id])
  }
}
plot_genes <- intersect(PREDATORY_GENES, colnames(plot_df))
plot_df <- plot_df |>
  pivot_longer(cols = all_of(plot_genes), names_to = "gene", values_to = "fpkm") |>
  mutate(log2_expr = if (max(fpkm, na.rm = TRUE) > 50) log2(fpkm + 1) else fpkm)

p_box <- ggplot(plot_df, aes(x = keap1_status, y = log2_expr, fill = keap1_status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.4) +
  facet_wrap(~ gene, scales = "free_y", ncol = 3) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y.npc = 0.95) +
  scale_fill_manual(values = c("KEAP1-WT" = "#3C5488", "KEAP1-MUT" = "#E64B35")) +
  labs(
    title = "Predatory Matrix Gene Expression",
    subtitle = "TCGA-LUAD: KEAP1-MUT vs KEAP1-WT (FPKM-UQ)",
    x = NULL, y = expression(log[2](FPKM-UQ + 1)), fill = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", strip.background = element_rect(fill = "grey95"))

ggsave(file.path(PATHS$figures, "boxplot_predatory_genes.pdf"),
       p_box, width = 12, height = 9, dpi = FIG_DPI)
ggsave(file.path(PATHS$figures, "boxplot_predatory_genes.png"),
       p_box, width = 12, height = 9, dpi = FIG_DPI)

# --- Module score ---
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
    subtitle = "Mean z-score: pump + efflux + scissors + uptake (10 genes)",
    x = NULL, y = "Module Score", fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

ggsave(file.path(PATHS$figures, "module_score_keap1_groups.pdf"),
       p_score, width = 5, height = 5, dpi = FIG_DPI)

fwrite(score_df, file.path(PATHS$results, "predatory_module_scores.csv"))
log_msg("Step 3 complete.")
