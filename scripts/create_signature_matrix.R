
library(omnideconv)
library(Seurat)
library(yaml)
# Load the config
config <- yaml::read_yaml("config.yaml")
# Access credentials
username <- config$cibersortx$username
token <- config$cibersortx$token
# Set credentials
omnideconv::set_cibersortx_credentials(username, token)

#load re-annotated macs data
macs <- readRDS("data/luca_query_reannotated.rds")
table(macs@meta.data$Projection_CellType)

macs_7_types <- c("LA_TAMs", "Angio_TAMs", "Inflam_TAMs", "Prolif_TAMs", "Reg_TAMs", "RTM_TAMs", "IFN_TAMs")

Idents(macs) <- "Projection_CellType"
macs <- subset(macs, idents = macs_7_types)
table(macs@meta.data$Projection_CellType)

process_disease <- function(macs, disease_name, output_filename, downsample_n = 200,
                            method = "cibersortx", g_min = 50, g_max = 100) {
  
  # Subset based on disease and primary tumor origin
  macs_sub <- subset(macs, subset = disease == disease_name)
  macs_sub <- subset(macs_sub, subset = origin == "tumor_primary")
  
  # Print cell type counts
  print(table(macs_sub@meta.data$Projection_CellType))
  
  # Downsample
  macs_sub <- subset(x = macs_sub, downsample = downsample_n)
  
  # Extract data and annotations
  norm_cnt <- macs_sub@assays[["RNA"]]@data
  ct <- macs_sub@meta.data$Projection_CellType
  
  # Build signature matrix
  sig_mtx <- omnideconv::build_model(
    single_cell_object = norm_cnt,
    cell_type_annotations = ct,
    method = method,
    g_min = g_min,
    g_max = g_max,
    verbose=TRUE
  )

  # Save to file
  saveRDS(sig_mtx, output_filename)
  
  return(sig_mtx)
}

sig_mtx_ad <- process_disease(
  macs = macs,
  disease_name = "lung adenocarcinoma",
  output_filename = "data/lung_adenocarcinoma_cibersortx_sig_mtrx.rds"
)

sig_mtx_sq <- process_disease(
  macs = macs,
  disease_name = "squamous cell lung carcinoma",
  output_filename = "data/squamous_cell_lung_carcinoma_cibersortx_sig_mtrx.rds"
)
