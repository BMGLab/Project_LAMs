# Load necessary libraries
library(dplyr)
library(readr)  
library(TCGAbiolinks)
library(SummarizedExperiment)

query.exp <- GDCquery(
  project = c("TCGA-LUAD", "TCGA-LUSC"),
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  access = "open"
)

GDCdownload(query.exp, directory = "data/")

data <- GDCprepare(query.exp, directory = "data/")

tpm_data <- assay(data, "tpm_unstrand")
saveRDS(tpm_data, "data/TCGA-Luad-Lusc-tpm_data.rds")

clinical_data <- colData(data)
saveRDS(clinical_data, "data/TCGA-Luad-Lusc-clinical_data.rds")



# Convert to data frame if needed
clean_clinical_data <- function(df) {
  clinical_columns <- c("primary_diagnosis", "vital_status", 
                        "days_to_last_follow_up", "ajcc_pathologic_stage")

  # Select available columns safely
  existing_columns <- clinical_columns[clinical_columns %in% colnames(df)]
  df <- df %>% dplyr::select(all_of(existing_columns))

  # Handle missing survival columns
  has_os <- "overall_survival" %in% colnames(df) && !all(is.na(df$overall_survival))
  has_death <- "days_to_death" %in% colnames(df) && !all(is.na(df$days_to_death))

  # Fill survival time
  if (has_os) {
    df$overall_survival <- as.numeric(df$overall_survival)
  } else if (has_death) {
    df$overall_survival <- as.numeric(df$days_to_death)
  } else {
    warning("No usable survival time found in 'overall_survival' or 'days_to_death'. Filling with NA.")
    df$overall_survival <- NA_real_
  }

  # Handle vital_status
  if (!"vital_status" %in% colnames(df)) {
    warning("'vital_status' column not found. Filling with NA.")
    df$vital_status <- NA
  }

  # Handle follow-up
  if (!"days_to_last_follow_up" %in% colnames(df)) {
    warning("'days_to_last_follow_up' column not found. Filling with NA.")
    df$days_to_last_follow_up <- NA_real_
  }

  # Final calculations
  df <- df %>%
    mutate(
      deceased = vital_status != "Alive",
      overall_survival = ifelse(vital_status == "Alive",
                                days_to_last_follow_up,
                                overall_survival),
      vital_status = as.factor(vital_status)
    )

  return(as.data.frame(df))
}



# Clean and save each type
only_primary_lusc_clean <- clean_clinical_data(only_primary_lusc)
saveRDS(only_primary_lusc_clean, "data/lusc_clinic_for_bulk.rds")

only_primary_luad_clean <- clean_clinical_data(only_primary_luad)
saveRDS(only_primary_luad_clean, "data/luad_clinic_for_bulk.rds")

# Prepare bulk TPM data
prepare_bulk_data <- function(clinical_df, tpm_data, output_path) {
  common_samples <- intersect(rownames(clinical_df), colnames(tpm_data))
  bulk_data <- tpm_data[, common_samples]
  rownames(bulk_data) <- sub("\\..*", "", rownames(bulk_data))
  saveRDS(bulk_data, output_path)
}

# Apply for each disease type
prepare_bulk_data(only_primary_luad_clean, tpm_data, "data/luad_bulk.rds")
prepare_bulk_data(only_primary_lusc_clean, tpm_data, "data/lusc_bulk.rds")
