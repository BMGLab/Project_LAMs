library(Seurat)
library(devtools)
library(HieRFIT)
library(tidyverse)
library(SingleCellExperiment)
library(zellkonverter)

luca_mod <- readRDS("data/luca_hierfit_training_data.rds")

library(dplyr)
library(ggplot2)
library(purrr)
library(forcats)

meta_df <- luca_mod@meta.data

clinical_features <- c("tumor_stage", "uicc_stage", "origin", "disease")

# Remove NAs for tumor_stage and uicc_stage
meta_df <- meta_df %>%
  filter(!is.na(tumor_stage), !is.na(uicc_stage))

# Count cells per sample per clinical category
count_df <- meta_df %>%
  select(sample, all_of(clinical_features)) %>%
  group_by(sample, across(all_of(clinical_features))) %>%
  summarise(cell_count = n(), .groups = "drop")

# Function to plot distribution for a given feature
plot_feature <- function(feature_name) {
  ggplot(count_df, aes_string(x = feature_name, y = "cell_count")) +
    geom_boxplot(outlier.alpha = 0.4, fill = "skyblue") +
    geom_jitter(width = 0.2, alpha = 0.5, color = "darkblue") +
    theme_bw() +
    labs(
      title = paste("Distribution of cell counts by", feature_name),
      x = feature_name,
      y = "Cell count"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

# Loop over features and save plots
walk(clinical_features, function(f) {
  p <- plot_feature(f)
  print(p)  # show in RStudio
  ggsave(filename = paste0("figures/cell_count_by_", f, ".pdf"),
         plot = p, width = 6, height = 4)
})
