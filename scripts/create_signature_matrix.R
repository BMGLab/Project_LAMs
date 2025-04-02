
library(omnideconv)
library(Seurat)
set_cibersortx_credentials("180323041@ogr.cbu.edu.tr", "b823e05f9e499673fdff7001ea895fa8")#Will be masked later on.

#load re-annotated macs data
macs <- readRDS("data/luca_query_reannotated.rds")
table(macs@meta.data$Projection_CellType)

macs_7_types <- c("LA_TAMs", "Angio_TAMs", "Inflam_TAMs", "Prolif_TAMs", "Reg_TAMs", "RTM_TAMs", "IFN_TAMs")

Idents(macs) <- "Projection_CellType"
macs <- subset(macs, idents = macs_7_types)
table(macs@meta.data$Projection_CellType)
library(Seurat)
library(omnideconv)

process_disease <- function(macs, disease_name, output_filename, downsample_n = 200,
                            method = "cibersortx", g_min = 50, g_max = 100) {
  
  set_cibersortx_credentials("180323041@ogr.cbu.edu.tr", "b823e05f9e499673fdff7001ea895fa8")#Will be masked later on.
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






celltype <- readRDS("~/Documents/luca_celltype_ann.rds")
matrx <- readRDS("~/Documents/luca_data_matrix.rds")
celltype <- as.vector(celltype)
class(celltype)

omnideconv::build_model_cibersortx()


omnideconv::set_cibersortx_credentials("180323041@ogr.cbu.edu.tr","b823e05f9e499673fdff7001ea895fa8")
signature_matrix_adeno <- omnideconv::build_model(single_cell_object = matrx, cell_type_annotations = celltype,method = "cibersortx")
saveRDS(signature_matrix_adeno, file = "~/Documents/luca-7Macs.cbrsrtx.signature.mtx.rds")

luadbulk <- readRDS("~/Documents/Adeno_bulk.rds")
luscbulk <- readRDS("~/Documents/Squamous_bulk.rds")

deconvolution_adeno <- omnideconv::deconvolute(luadbulk, signature_matrix_adeno, method = "cibersortx")

saveRDS(deconvolution_adeno, file = "~/Documents/TCGA-LUAD-7Macs.deconvoluted.rds")


deconvolution_squamoz <- omnideconv::deconvolute(luscbulk, signature_matrix_adeno, method = "cibersortx")

saveRDS(deconvolution_squamoz, file = "~/Documents/TCGA-LUSC-7Macs.deconvoluted.rds")