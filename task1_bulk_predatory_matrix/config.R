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
  geo       = file.path(PROJECT_ROOT, "data", "geo")
)

for (p in PATHS) dir.create(p, recursive = TRUE, showWarnings = FALSE)

# --- Predatory matrix gene set (dual metabolic monopoly hypothesis) ---
PREDATORY_GENES <- c(
  "SLC7A11",  # cystine pump
  "GGT1",     # GSH shredder
  "SLC1A5",   # ASCT2, cysteine re-uptake
  "ABCC1", "ABCC2", "ABCC3"  # GSH efflux transporters
)

# NRF2 pathway anchor genes (optional validation panel)
NRF2_TARGETS <- c("NFE2L2", "HMOX1", "NQO1", "GCLC", "GCLM", "TXNRD1")

# TCGA settings
TCGA_PROJECT <- "TCGA-LUAD"
MIN_TUMOR_SAMPLES_MUT <- 5   # minimum KEAP1-MUT samples for stable stats
MIN_TPM <- 1                 # expression filter threshold
MIN_SAMPLE_FRACTION <- 0.10  # gene must be expressed in >= 10% samples

# DEA parameters
DEA_PADJ <- 0.05
DEA_LOG2FC <- log2(1.25)     # modest effect size for metabolic genes

# Correlation
COR_METHOD <- "spearman"
COR_PADJ_METHOD <- "BH"

# GEO validation cohorts (KEAP1 / NRF2 / LUAD related)
GEO_DATASETS <- list(
  # KEAP1 CRISPR KO vs control in H1299 (lung adenocarcinoma cell line)
  GSE142694 = list(
    title = "KEAP1 knockout H1299 RNA-seq",
    group_col_hint = NULL,  # auto-detect from pData
    mut_label = "KEAP1_KO",
    wt_label  = "Control"
  ),
  # Independent LUAD bulk cohort (prognostic, no mutation — used for co-expression only)
  GSE68465 = list(
    title = "LUAD bulk RNA-seq (co-expression validation)",
    use_for = "correlation_only"
  )
)

# Figure theme
FIG_DPI <- 300
FIG_WIDTH <- 8
FIG_HEIGHT <- 7
