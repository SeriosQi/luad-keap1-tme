#!/usr/bin/env Rscript
# =============================================================================
# Master pipeline — Task 1: Bulk RNA-seq Predatory Matrix Analysis
# =============================================================================
# Usage:
#   cd /home/xiruisi9394/luad_keap1_tme/task1_bulk_predatory_matrix
#   Rscript run_task1_pipeline.R [--skip-download] [--skip-geo] [--gdc]
# Default download: Scheme A (GDC Xena Hub, 2 files)
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
SKIP_DOWNLOAD <- "--skip-download" %in% args
SKIP_GEO      <- "--skip-geo" %in% args
USE_GDC       <- "--gdc" %in% args

PROJECT_ROOT <- "/home/xiruisi9394/luad_keap1_tme/task1_bulk_predatory_matrix"
setwd(PROJECT_ROOT)

source("config.R")
source("R/utils.R")

log_msg("========================================")
log_msg("LUAD KEAP1 TME — Task 1 Pipeline Start")
log_msg("========================================")

download_step <- if (USE_GDC) "R/01_download_tcga.R" else "R/01_download_tcga_xena.R"

steps <- c(
  download_step,
  "R/02_keap1_mutation_status.R",
  "R/03_differential_expression.R",
  "R/04_correlation_heatmap.R",
  "R/05_geo_validation.R",
  "R/06_h358_predatory_matrix.R"
)

if (SKIP_DOWNLOAD) {
  steps <- steps[-1]
  log_msg("Skipping TCGA download (using cached data)")
}
if (SKIP_GEO) {
  steps <- steps[steps != "R/05_geo_validation.R"]
  log_msg("Skipping GEO validation")
}

for (step in steps) {
  log_msg("Running: ", step)
  tryCatch(
    source(step, local = new.env()),
    error = function(e) {
      log_msg("ERROR in ", step, ": ", conditionMessage(e))
      stop(e)
    }
  )
}

save_session_info(file.path(PATHS$results, "session_info.txt"))

log_msg("========================================")
log_msg("Pipeline complete. Key outputs:")
log_msg("  DEA table       : results/dea_keap1_mut_vs_wt.csv")
log_msg("  Predatory DEA   : results/dea_predatory_matrix_genes.csv")
log_msg("  Correlation     : results/correlation_r_KEAP1-MUT.csv")
log_msg("  Heatmap (MUT)   : results/figures/heatmap_predatory_correlation_KEAP1_MUT.pdf")
log_msg("  Module score    : results/figures/module_score_keap1_groups.pdf")
log_msg("  Volcano plot    : results/figures/volcano_predatory_matrix.pdf")
log_msg("========================================")
