# =============================================================================
# Step 2: Classify TCGA-LUAD samples into KEAP1-MUT vs KEAP1-WT
# =============================================================================
# Usage: Rscript R/02_keap1_mutation_status.R
# =============================================================================

source("config.R")
source("R/utils.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

log_msg("=== Step 2: KEAP1 mutation classification ===")

maf       <- readRDS(file.path(PATHS$processed, "tcga_luad_maf.rds"))
tpm       <- readRDS(file.path(PATHS$processed, "tcga_luad_tpm_matrix.rds"))
out_file  <- file.path(PATHS$processed, "keap1_status.csv")

# Non-silent variant classes (standard TCGA MAF annotation)
nonsilent_classes <- c(
  "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
  "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
  "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation"
)

maf_dt <- as.data.table(maf)

# Identify KEAP1 mutant patients
keap1_mut <- maf_dt[
  Hugo_Symbol == "KEAP1" &
    Variant_Classification %in% nonsilent_classes &
    !is.na(Tumor_Sample_Barcode)
]
keap1_mut_patients <- unique(normalize_patient_id(keap1_mut$Tumor_Sample_Barcode))

log_msg("KEAP1 non-silent mutations found in ", length(keap1_mut_patients), " patients")

# Map TPM columns to patient IDs
sample_ids  <- colnames(tpm)
patient_ids <- normalize_patient_id(sample_ids)

status_df <- data.frame(
  sample_id  = sample_ids,
  patient_id = patient_ids,
  keap1_status = ifelse(patient_ids %in% keap1_mut_patients, "KEAP1-MUT", "KEAP1-WT"),
  stringsAsFactors = FALSE
)

# Optional: flag concurrent NFE2L2 (NRF2) mutations
nrf2_mut <- maf_dt[
  Hugo_Symbol == "NFE2L2" &
    Variant_Classification %in% nonsilent_classes
]
nrf2_patients <- unique(normalize_patient_id(nrf2_mut$Tumor_Sample_Barcode))
status_df$nfe2l2_status <- ifelse(status_df$patient_id %in% nrf2_patients, "NFE2L2-MUT", "NFE2L2-WT")
status_df$nrf2_keap1_co <- status_df$keap1_status == "KEAP1-MUT" & status_df$nfe2l2_status == "NFE2L2-MUT"

n_mut <- sum(status_df$keap1_status == "KEAP1-MUT")
n_wt  <- sum(status_df$keap1_status == "KEAP1-WT")
log_msg("Tumor samples — KEAP1-MUT: ", n_mut, " | KEAP1-WT: ", n_wt)

if (n_mut < MIN_TUMOR_SAMPLES_MUT) {
  warning(
    "KEAP1-MUT sample count (", n_mut, ") is below recommended minimum (",
    MIN_TUMOR_SAMPLES_MUT, "). Consider expanding to pan-lung or using GEO validation."
  )
}

# Save mutation detail table for supplementary
keap1_detail <- keap1_mut[, .(
  patient_id = normalize_patient_id(Tumor_Sample_Barcode),
  Variant_Classification,
  HGVSp_Short,
  HGVSc,
  Tumor_Sample_Barcode
)]
fwrite(keap1_detail, file.path(PATHS$results, "keap1_mutation_details.csv"))
fwrite(status_df, out_file)

log_msg("Saved sample status: ", out_file)
log_msg("Step 2 complete.")
