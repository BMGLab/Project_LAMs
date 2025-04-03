# Load necessary libraries
library(dplyr)
library(readr)
library(SummarizedExperiment)


if (!file.exists("data/TCGA-Luad-Lusc-clinical_data.rds") | !file.exists("data/TCGA-Luad-Lusc-tpm_data.rds") ) {
  print("need to download from TCGA!")
  library(TCGAbiolinks)

  query.exp <- GDCquery(
    project = c("TCGA-LUAD", "TCGA-LUSC"),
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    access = "open")

  GDCdownload(query.exp, directory = "data/")

  data <- GDCprepare(query.exp, directory = "data/")

  tpm_data <- assay(data, "tpm_unstrand")
  saveRDS(tpm_data, "data/TCGA-Luad-Lusc-tpm_data.rds")
  clinical_data <- colData(data)
  saveRDS(clinical_data, "data/TCGA-Luad-Lusc-clinical_data.rds")
}else{
  tpm_data <- readRDS("data/TCGA-Luad-Lusc-tpm_data.rds")
  clinical_data <- readRDS("data/TCGA-Luad-Lusc-clinical_data.rds")
}


# Vector of diagnoses
pr_diag <- c("Squamous cell carcinoma, NOS", "Adenocarcinoma, NOS")

# Define a general function
process_clinical_data <- function(df, diagnosis, output_path) {
  processed <- df %>%
    as.data.frame() %>%
    mutate(
      Tumor.stage = ifelse(paper_Tumor.stage == "[Not Available]", NA, paper_Tumor.stage),
      Tumor.stage = ifelse(is.na(ajcc_pathologic_stage), NA, ajcc_pathologic_stage)
    ) %>%
    filter(!is.na(Tumor.stage)) %>%
    filter(sample_type == "Primary Tumor") %>%
    filter(primary_diagnosis == diagnosis) %>%
    mutate(
      overall_survival = ifelse(
        vital_status == "Alive",
        as.numeric(days_to_last_follow_up),
        as.numeric(days_to_death)
      ),
      deceased = vital_status != "Alive",
      vital_status = as.factor(vital_status)
    ) %>%
    select(primary_diagnosis, vital_status, days_to_last_follow_up,
           overall_survival, Tumor.stage, deceased)

  saveRDS(processed, output_path)
  return(processed)
}

# Prepare bulk TPM data
prepare_bulk_data <- function(clinical_df, tpm_data, output_path) {
  common_samples <- intersect(rownames(clinical_df), colnames(tpm_data))
  bulk_data <- tpm_data[, common_samples]
  rownames(bulk_data) <- sub("\\..*", "", rownames(bulk_data))
  saveRDS(bulk_data, output_path)
}

# Run for each diagnosis
for (diag in pr_diag) {
  label <- if (grepl("Squamous", diag)) "sq" else "ad"
  cln_outfile <- paste0("data/clinic_data_", label, ".rds")
  exp_outfile <- paste0("data/bulk_data_", label, ".rds")
  cln_df <- process_clinical_data(clinical_data, diagnosis = diag, output_path = cln_outfile)
  prepare_bulk_data(cln_df, tpm_data = tpm_data, exp_outfile)
}
