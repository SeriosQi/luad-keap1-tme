# =============================================================================
# Install R dependencies for Task 1
# Run once: Rscript install_dependencies.R
# =============================================================================

cran_pkgs <- c(
  "data.table", "dplyr", "tibble", "tidyr", "ggplot2", "ggpubr",
  "ggrepel", "pheatmap", "circlize", "RColorBrewer"
)

bioc_pkgs <- c(
  "TCGAbiolinks", "SummarizedExperiment", "DESeq2", "limma", "edgeR",
  "ComplexHeatmap", "maftools", "GEOquery", "org.Hs.eg.db", "AnnotationDbi"
)

install_cran <- function(pkgs) {
  missing <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(missing) > 0) {
    message("Installing CRAN: ", paste(missing, collapse = ", "))
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
}

install_bioc <- function(pkgs) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  missing <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(missing) > 0) {
    message("Installing Bioconductor: ", paste(missing, collapse = ", "))
    BiocManager::install(missing, ask = FALSE, update = FALSE)
  }
}

install_cran(cran_pkgs)
install_bioc(bioc_pkgs)

message("All dependencies ready.")
