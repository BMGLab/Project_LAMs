# Load required libraries
library(survminer)
library(survival)


# Define the cell types
macs_7_types <- c("LA_TAMs", "Angio_TAMs", "Inflam_TAMs", "Prolif_TAMs", "Reg_TAMs", "RTM_TAMs", "IFN_TAMs")

# Function to perform survival analysis and save plots
perform_survival_analysis <- function(deconv_file, clinic_file, output_pdf) {
  
  # Load data
  deconv_data <- readRDS(deconv_file)
  clinic_data <- readRDS(clinic_file)
  
  # Find common samples
  common_samples <- intersect(rownames(clinic_data), rownames(deconv_data))
  deconv_data <- deconv_data[common_samples, ]
  clinic_data <- clinic_data[common_samples, ]
  surv.data <- cbind(deconv_data, clinic_data)
  
  # Open a PDF to save plots
  pdf(output_pdf, width = 8, height = 6)
  
  for (cell_type in macs_7_types) {
    cat("\nProcessing:", cell_type, "\n")
    
    # Try to compute cutpoint
    res.cut <- tryCatch({
      surv_cutpoint(surv.data, time = "overall_survival", event = "deceased", variables = cell_type)
    }, error = function(e) NULL)
    
    # Skip if no valid cutpoint
    if (is.null(res.cut)) {
      cat("Skipping:", cell_type, "due to cutpoint computation error.\n")
      next
    }
    
    cutpoint <- res.cut[["cutpoint"]][["cutpoint"]]
    surv.data$MacStatus <- ifelse(surv.data[[cell_type]] > cutpoint, "High", "Low")
    
    # Fit survival model
    fit_MacStatus <- survfit(Surv(overall_survival, deceased) ~ MacStatus, data = surv.data)
    
    # Generate and save the survival plot
    plot <- ggsurvplot(fit_MacStatus, 
                       surv.data, 
                       pval = TRUE,    
                       risk.table = TRUE,
                       xlab = "Time (days)",
                       title = paste("Survival Analysis for", cell_type))
    
    print(plot)  # Print the plot to save it in the PDF
  }
  
  dev.off()  # Close the PDF device
}

# Run the function for LUSC dataset
perform_survival_analysis("data/macs_deconved_bulk_sq.rds",
                          "data/clinic_data_sq.rds",
                          "figures/LUSC_Survival_Plots.pdf")

# Run the function for LUAD dataset
perform_survival_analysis("data/macs_deconved_bulk_ad.rds",
                          "data/clinic_data_ad.rds",
                          "figures/LUAD_Survival_Plots.pdf")

