# Load necessary libraries
library(HPAanalyze)
library(dplyr)
library(ggplot2)

# Create output directory if it doesn't exist
if (!dir.exists("figures")) dir.create("figures")

# Read input CSV
deg_data <- read.csv("results/macs_PB-DEGs_repeated.csv", stringsAsFactors = FALSE)

# Filter for SubType = LA_TAMs
la_genes <- deg_data %>%
  filter(SubType == "LA_TAMs") %>%
  distinct(gene_name) %>%
  pull(gene_name)

# Generate and save cancer plot
pdf("figures/histopathology_lung_cancer_LA_TAMs.pdf", width = 10, height = 3)
hpaVis(targetGene = la_genes,
       visType = "Cancer",
       targetCancer = "lung cancer")
dev.off()

# Generate and save normal lung tissue plot
pdf("figures/histopathology_lung_tissue_LA_TAMs.pdf", width = 10, height = 6)
hpaVisTissue(targetGene = la_genes,
             targetTissue = "lung")
dev.off()
