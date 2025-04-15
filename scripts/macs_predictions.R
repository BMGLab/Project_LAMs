library(Seurat)
library(devtools)
library(HieRFIT)
library(DiagrammeR)
library(alluvial)
library(ggalluvial)
library(tidyverse)
library(SingleCellExperiment)
library(zellkonverter)

#load whole LuCa dataset
luca <- readRDS("data/local.rds")

#Subsetting the only macrophages
luca_macs <- subset(luca, ann_fine %in% c("Macrophage alveolar", "Macrophage"))
saveRDS(luca_macs, "data/luca_macs_subset.rds")

#Labels with max AUCell scores
subtypes <- read.csv("data/subtype_matrix.csv")

colnames(subtypes)[colnames(subtypes) == "X"] <- "Barcodes"

luca_m <- subset(luca_macs, cells = subtypes$Barcodes)

luca_m@meta.data <- cbind(luca_m@meta.data, subtypes)


#Barcodes of top250 cells from each subtypes
Angio_TAMs_barcodes <- read_csv("data/Angio_TAMs_barcodes.csv")
IFN_TAMs_barcodes <- read_csv("data/IFN_TAMs_barcodes.csv")
Inflam_TAMs_barcodes <- read_csv("data/Inflam_TAMs_barcodes.csv")
LA_TAMs_barcodes <- read_csv("data/LA_TAMs_barcodes.csv")
Prolif_TAMs_barcodes <- read_csv("data/Prolif_TAMs_barcodes.csv")
Reg_TAMs_barcodes <- read_csv("data/Reg_TAMs_barcodes.csv")
RTM_TAMs_barcodes <- read_csv("data/RTM_TAMs_barcodes.csv")

barcodes <- list(Angio_TAMs_barcodes, IFN_TAMs_barcodes, Inflam_TAMs_barcodes, LA_TAMs_barcodes, 
                 Prolif_TAMs_barcodes, Reg_TAMs_barcodes, RTM_TAMs_barcodes)

barcodes_df <- do.call(rbind, barcodes)

#### Training and Query Dataset for HieRFIT Model

#Training Dataset;
luca_mod <- subset(luca_m, cells = barcodes_df$Barcodes)

#Query Dataset;
luca_new <- subset(luca_m, cells = setdiff(Cells(luca_m), barcodes_df$Barcodes))


#Saving the objects
saveRDS(luca_mod, file = "data/luca_hierfit_training_data.rds")
saveRDS(luca_new, file = "data/luca_hierfit_query_data.rds")

#Create a hierfit object
refmod <- CreateHieR(RefData = luca_mod[["RNA"]]@data,
                     ClassLabels = luca_mod@meta.data$Subtype,
                     species = "hsapiens")

#Saving the model
SaveHieRMod(refMod = refmod, filePrefix = "data/Luca_Subtype_HierMod")

### Reannotation of Query Data with HieRFIT Model;
#Project the cell class labels on the new dataset:
hierObj <- HieRFIT(Query = luca_new[["RNA"]]@data, refMod = refmod)


p.pred <- PlotBarStats(HieRobj = hierObj)

ggsave(p.pred, filename = "figures/macs.predictions.barplot.pdf", width = 8, height = 4, device = "pdf")


luca_new@meta.data <- cbind(luca_new@meta.data, hierObj@Evaluation$Projection)
colnames(luca_new@meta.data)[colnames(luca_new@meta.data) == "hierObj@Evaluation$Projection"] <- "Projection_CellType"

saveRDS(luca_new, file = "data/luca_query_reannotated.rds")

# Load your Seurat object
luca_new <- readRDS("data/luca_query_reannotated.rds")
sce <- as.SingleCellExperiment(luca_new)
writeH5AD(sce, file = "/data/luca_query_reannotated.h5ad")
