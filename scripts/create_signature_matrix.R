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

createSignature <- function(seu, disease_name, output_filename,
                            downsample_n = 200,
                            method = "cibersortx",
                            g_min = 50,
                            g_max = 100,
                            celltype_col = "Projection_CellType") {
  
  # Subset based on disease and primary tumor origin
  seu_sub <- subset(seu, subset = disease == disease_name)
  seu_sub <- subset(seu_sub, subset = origin == "tumor_primary")
  
  # Downsample
  seu_sub <- subset(x = seu_sub, downsample = downsample_n)
  
  # Extract data and dynamic cell type annotations
  norm_cnt <- seu_sub@assays[["RNA"]]@data
  
  if (!celltype_col %in% colnames(seu_sub@meta.data)) {
    stop(paste("Cell type column", celltype_col, "not found in meta.data"))
  }
  
  ct <- seu_sub@meta.data[[celltype_col]]
  
  # Build signature matrix
  sig_mtx <- omnideconv::build_model(
    single_cell_object = norm_cnt,
    cell_type_annotations = ct,
    method = method,
    g_min = g_min,
    g_max = g_max,
    verbose = TRUE
  )
  
  # Save to file
  saveRDS(sig_mtx, output_filename)
  
  return(sig_mtx)
}

#load re-annotated macs data
macs <- readRDS("data/luca_query_reannotated.rds")
table(macs@meta.data$Projection_CellType)

macs_7_types <- c("LA_TAMs", "Angio_TAMs", "Inflam_TAMs", "Prolif_TAMs", "Reg_TAMs", "RTM_TAMs", "IFN_TAMs")

Idents(macs) <- "Projection_CellType"
macs <- subset(macs, idents = macs_7_types)
print(table(macs@meta.data$Projection_CellType))

print("Creating a signature matrix for lung adenocarcinoma macrophages...")
sig_mtx_ad <- createSignature(
  seu = macs,
  disease_name = "lung adenocarcinoma",
  output_filename = "data/lung_adenocarcinoma_cibersortx_macs_sig_mtrx.rds",
  celltype_col = "Projection_CellType"
)

print("Creating a signature matrix for squamous cell lung carcinoma macrophages...")
sig_mtx_sq <- createSignature(
  seu = macs,
  disease_name = "squamous cell lung carcinoma",
  output_filename = "data/squamous_cell_lung_carcinoma_macs_cibersortx_sig_mtrx.rds",
  celltype_col = "Projection_CellType"
)
