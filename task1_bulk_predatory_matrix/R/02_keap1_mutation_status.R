# =============================================================================
# Step 2: Classify TCGA-LUAD samples into KEAP1-MUT vs KEAP1-WT
# Supports: Xena somaticmutation (Scheme A) or GDC MAF (Scheme B/GDC)
# =============================================================================

source("config.R")
source("R/utils.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

log_msg("=== Step 2: KEAP1 mutation classification ===")

tpm      <- readRDS(file.path(PATHS$processed, "tcga_luad_tpm_matrix.rds"))
out_file <- file.path(PATHS$processed, "keap1_status.csv")
src_flag <- file.path(PATHS$processed, "data_source.txt")
data_src <- if (file.exists(src_flag)) readLines(src_flag, n = 1) else "gdc"

nonsilent_pattern <- paste0(
  "Missense_Mutation|Nonsense_Mutation|Frame_Shift|In_Frame|Splice_Site|",
  "Translation_Start_Site|Nonstop_Mutation|Fusion|Deletion|Insertion"
)

get_keap1_mut_patients_xena <- function(mut_dt) {
  nonsilent_effects <- paste0(
    "missense|stop_gained|frameshift|splice|inframe|start_lost|stop_lost"
  )

  # Long format (GDC Xena somaticmutation_wxs): sample + gene + effect columns
  if ("gene" %in% colnames(mut_dt) && "sample" %in% colnames(mut_dt)) {
    keap1_rows <- mut_dt[
      gene == "KEAP1" &
        grepl(nonsilent_effects, effect, ignore.case = TRUE) &
        !grepl("synonymous", effect, ignore.case = TRUE)
    ]
    mut_samples <- unique(keap1_rows$sample)
    mut_patients <- unique(normalize_patient_id(mut_samples))
    keap1_detail <- as.data.frame(keap1_rows[, .(
      patient_id = normalize_patient_id(sample),
      mutation_annotation = effect,
      sample = sample
    )])
    nrf2_patients <- unique(normalize_patient_id(
      mut_dt[gene == "NFE2L2" & grepl(nonsilent_effects, effect, ignore.case = TRUE), sample]
    ))
    return(list(
      mut_patients = mut_patients,
      keap1_detail = keap1_detail,
      nrf2_patients = nrf2_patients
    ))
  }

  # Wide format fallback: samples x genes
  gene_col <- colnames(mut_dt)[1]
  if ("KEAP1" %in% colnames(mut_dt)) {
    samp_ids <- mut_dt[[gene_col]]
    keap1_vals <- mut_dt[["KEAP1"]]
    mut_samples <- samp_ids[
      !is.na(keap1_vals) & keap1_vals != "" &
        grepl(nonsilent_pattern, keap1_vals, ignore.case = TRUE)
    ]
    mut_patients <- unique(normalize_patient_id(mut_samples))
    keap1_detail <- data.frame(
      patient_id = normalize_patient_id(mut_samples),
      mutation_annotation = keap1_vals[
        !is.na(keap1_vals) & keap1_vals != "" &
          grepl(nonsilent_pattern, keap1_vals, ignore.case = TRUE)
      ],
      stringsAsFactors = FALSE
    )
    nrf2_patients <- character(0)
    if ("NFE2L2" %in% colnames(mut_dt)) {
      nrf2_vals <- mut_dt[["NFE2L2"]]
      nrf2_patients <- unique(normalize_patient_id(samp_ids[
        !is.na(nrf2_vals) & grepl(nonsilent_pattern, nrf2_vals, ignore.case = TRUE)
      ]))
    }
    return(list(
      mut_patients = mut_patients,
      keap1_detail = keap1_detail,
      nrf2_patients = nrf2_patients
    ))
  }

  stop("KEAP1 not found in Xena mutation table.")
}

if (data_src == "xena_gdc_hub") {
  mut_dt <- readRDS(file.path(PATHS$processed, "xena_mutation_dt.rds"))
  parsed <- get_keap1_mut_patients_xena(mut_dt)
  keap1_mut_patients <- parsed$mut_patients
  keap1_detail <- parsed$keap1_detail
  nrf2_patients <- parsed$nrf2_patients
  log_msg("Data source: GDC Xena Hub somaticmutation_wxs")
} else {
  maf <- readRDS(file.path(PATHS$processed, "tcga_luad_maf.rds"))
  nonsilent_classes <- c(
    "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
    "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
    "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation"
  )
  maf_dt <- as.data.table(maf)
  keap1_mut <- maf_dt[
    Hugo_Symbol == "KEAP1" &
      Variant_Classification %in% nonsilent_classes &
      !is.na(Tumor_Sample_Barcode)
  ]
  keap1_mut_patients <- unique(normalize_patient_id(keap1_mut$Tumor_Sample_Barcode))
  keap1_detail <- keap1_mut[, .(
    patient_id = normalize_patient_id(Tumor_Sample_Barcode),
    Variant_Classification,
    HGVSp_Short,
    HGVSc,
    Tumor_Sample_Barcode
  )]
  nrf2_mut <- maf_dt[
    Hugo_Symbol == "NFE2L2" &
      Variant_Classification %in% nonsilent_classes
  ]
  nrf2_patients <- unique(normalize_patient_id(nrf2_mut$Tumor_Sample_Barcode))
  log_msg("Data source: GDC MAF")
}

log_msg("KEAP1 non-silent mutations found in ", length(keap1_mut_patients), " patients")

sample_ids  <- colnames(tpm)
patient_ids <- normalize_patient_id(sample_ids)

status_df <- data.frame(
  sample_id    = sample_ids,
  patient_id   = patient_ids,
  keap1_status = ifelse(patient_ids %in% keap1_mut_patients, "KEAP1-MUT", "KEAP1-WT"),
  stringsAsFactors = FALSE
)
status_df$nfe2l2_status <- ifelse(patient_ids %in% nrf2_patients, "NFE2L2-MUT", "NFE2L2-WT")
status_df$nrf2_keap1_co <- status_df$keap1_status == "KEAP1-MUT" & status_df$nfe2l2_status == "NFE2L2-MUT"

n_mut <- sum(status_df$keap1_status == "KEAP1-MUT")
n_wt  <- sum(status_df$keap1_status == "KEAP1-WT")
log_msg("Tumor samples — KEAP1-MUT: ", n_mut, " | KEAP1-WT: ", n_wt)

if (n_mut < MIN_TUMOR_SAMPLES_MUT) {
  warning(
    "KEAP1-MUT sample count (", n_mut, ") is below recommended minimum (",
    MIN_TUMOR_SAMPLES_MUT, ")."
  )
}

fwrite(as.data.table(keap1_detail), file.path(PATHS$results, "keap1_mutation_details.csv"))
fwrite(status_df, out_file)

log_msg("Saved sample status: ", out_file)
log_msg("Step 2 complete.")
