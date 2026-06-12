# =============================================================================
# LUAD KEAP1 TME — Task 1: Bulk RNA-seq "Predatory Matrix" Analysis
# Configuration
# =============================================================================

PROJECT_ROOT <- "/home/xiruisi9394/luad_keap1_tme/task1_bulk_predatory_matrix"

PATHS <- list(
  raw       = file.path(PROJECT_ROOT, "data", "raw"),
  processed = file.path(PROJECT_ROOT, "data", "processed"),
  results   = file.path(PROJECT_ROOT, "results"),
  figures   = file.path(PROJECT_ROOT, "results", "figures"),
  geo       = file.path(PROJECT_ROOT, "data", "geo"),
  h358      = file.path(PROJECT_ROOT, "data", "raw", "h358")
)

for (p in PATHS) dir.create(p, recursive = TRUE, showWarnings = FALSE)

# --- Predatory matrix gene set (expanded Task 1 / 战役 1) ---
PREDATORY_GENES <- c(
  "SLC7A11",                    # Cystine pump (主力抽水机)
  "ABCC1", "ABCC2", "ABCC3",   # GSH efflux (外排泄露泵)
  "GGT1", "GGT5", "DPEP1",      # GSH scissors (废品回收剪刀)
  "SLC1A4", "SLC1A5"           # Cysteine uptake — ASCT1 / ASCT2
)

PREDATORY_GENE_MODULE <- c(
  "SLC7A11" = "Cystine pump",
  "ABCC1"   = "GSH efflux", "ABCC2" = "GSH efflux", "ABCC3" = "GSH efflux",
  "GGT1"    = "GSH scissors", "GGT5" = "GSH scissors", "DPEP1" = "GSH scissors",
  "SLC1A4"  = "Cysteine uptake", "SLC1A5" = "Cysteine uptake"
)

# Key co-expression pairs for mechanism testing
KEY_GENE_PAIRS <- list(
  c("SLC7A11", "GGT1"),   # pump + shredder
  c("SLC7A11", "GGT5"),
  c("GGT1", "GGT5"),
  c("GGT1", "DPEP1"),
  c("GGT5", "DPEP1"),
  c("SLC7A11", "SLC1A5"),
  c("SLC7A11", "SLC1A4"),
  c("GGT1", "SLC1A5"),
  c("SLC7A11", "ABCC1"),
  c("ABCC1", "ABCC2")
)

# NRF2 pathway anchor genes (optional validation panel)
NRF2_TARGETS <- c("NFE2L2", "HMOX1", "NQO1", "GCLC", "GCLM", "TXNRD1")

# Data source: "xena" (Scheme A, default) or "gdc" (TCGAbiolinks STAR counts)
DATA_SOURCE <- "xena"

# GDC Xena Hub — Scheme A (2-file fast download)
XENA_BASE <- "https://gdc.xenahubs.net/download"
XENA_FILES <- list(
  expression = list(
    url  = paste0(XENA_BASE, "/TCGA-LUAD.star_fpkm-uq.tsv.gz"),
    dest = "TCGA-LUAD.star_fpkm-uq.tsv.gz"
  ),
  mutation = list(
    url  = paste0(XENA_BASE, "/TCGA-LUAD.somaticmutation_wxs.tsv.gz"),
    dest = "TCGA-LUAD.somaticmutation_wxs.tsv.gz"
  )
)

# H358 WT vs KEAP1-KO — place your lab matrices here (see templates/h358/)
H358 <- list(
  expr_file  = file.path(PATHS$h358, "expression_matrix.csv"),
  meta_file  = file.path(PATHS$h358, "sample_metadata.csv"),
  cell_line  = "H358",
  wt_label   = "WT",
  ko_label   = "KEAP1_KO",
  # Public reference cohort (ArrayExpress E-MTAB-9724, H358 KEAP1 KO clones)
  reference  = "E-MTAB-9724"
)

# TCGA settings
TCGA_PROJECT <- "TCGA-LUAD"
MIN_TUMOR_SAMPLES_MUT <- 5
MIN_TPM <- 1
MIN_SAMPLE_FRACTION <- 0.10

# DEA parameters
DEA_PADJ <- 0.05
DEA_LOG2FC <- log2(1.25)

# Correlation
COR_METHOD <- "spearman"
COR_PADJ_METHOD <- "BH"

# GEO validation cohorts
GEO_DATASETS <- list(
  GSE142694 = list(
    title = "KEAP1 knockout H1299 RNA-seq",
    group_col_hint = NULL,
    mut_label = "KEAP1_KO",
    wt_label  = "Control"
  ),
  GSE68465 = list(
    title = "LUAD bulk RNA-seq (co-expression validation)",
    use_for = "correlation_only"
  )
)

# Figure theme
FIG_DPI <- 300
FIG_WIDTH <- 10
FIG_HEIGHT <- 8

# Heatmap module colors
MODULE_COLORS <- c(
  "Cystine pump"     = "#E64B35",
  "GSH efflux"       = "#4DBBD5",
  "GSH scissors"     = "#00A087",
  "Cysteine uptake"  = "#3C5488"
)
