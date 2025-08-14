library(Seurat)
library(devtools)
library(HieRFIT)
#library(DiagrammeR)
#library(alluvial)
#library(ggalluvial)
library(tidyverse)
library(SingleCellExperiment)
library(zellkonverter)
#Load whole LuCa dataset
luca_macs <- readRDS("data/luca_macs_subset.rds")

#subset the macrophages for COPD and SCLC
luca_new <- subset(luca_macs, disease %in% c("chronic obstructive pulmonary disease", "SCLC"))
#number of donors
luca_new@meta.data[,c("donor_id", "disease")] %>% unique %>% select(disease) %>% table
#number of samples
luca_new@meta.data[,c("sample", "disease")] %>% unique %>% select(disease) %>% table

#Load hierfit model

refmod <- readRDS("data/Luca_Subtype_HierMod.RDS")

### Reannotation of Query Data with HieRFIT Model;
#Project the cell class labels on the new dataset:
hierObj <- HieRFIT(Query = luca_new[["RNA"]]@data, refMod = refmod)


p.pred <- PlotBarStats(HieRobj = hierObj)

ggsave(p.pred, filename = "figures/COPD.macs.predictions.barplot.pdf", width = 8, height = 4, device = "pdf")


luca_new@meta.data <- cbind(luca_new@meta.data, hierObj@Evaluation$Projection)
colnames(luca_new@meta.data)[colnames(luca_new@meta.data) == "hierObj@Evaluation$Projection"] <- "Projection_CellType"


saveRDS(luca_new, file = "data/COPD.luca_query_reannotated.rds")

# Load your Seurat object
luca_new <- readRDS("data/COPD.luca_query_reannotated.rds")
sce <- as.SingleCellExperiment(luca_new)
writeH5AD(sce, file = "/data/COPD.luca_query_reannotated.h5ad")


#composition of macrophages in COPD
luca_new <- readRDS("data/COPD.luca_query_reannotated.rds")

#Compostion of macrophages in COPD gender comparison
p.ea <- luca_new@meta.data %>% 
  select(sample, Projection_CellType, development_stage, sex, ever_smoker, disease) %>% 
  group_by(sample, Projection_CellType) %>% 
  mutate(Cnt=n()) %>% 
  unique() %>% 
  group_by(sample) %>% 
  mutate(pct=100*Cnt/sum(Cnt), Less=ifelse(sum(Cnt)<50, "low", "ok")) %>% 
  filter(Less == "ok") %>% 
  filter(!Projection_CellType %in% c('Int.Node.3', 'Int.Node.4', 'Int.Node.5')) %>%
  filter(sex %in% c('female', 'male')) %>%
  ggplot(aes(Projection_CellType, pct, fill = sex )) + geom_boxplot() + theme_bw()

ggsave(p.ea, filename = "figures/COPD.pct.sex.pdf", width = 8, height = 4, device = "pdf")

#composition of macrophages in COPD ever smoker comparison
p.ea <- luca_new@meta.data %>% 
  select(sample, Projection_CellType, development_stage, sex, ever_smoker, disease) %>% 
  group_by(sample, Projection_CellType) %>% 
  mutate(Cnt=n()) %>% 
  unique() %>% 
  group_by(sample) %>% 
  mutate(pct=100*Cnt/sum(Cnt), Less=ifelse(sum(Cnt)<50, "low", "ok")) %>% 
  filter(Less == "ok") %>% 
  filter(!Projection_CellType %in% c('Int.Node.3', 'Int.Node.4', 'Int.Node.5')) %>%
  filter(ever_smoker %in% c('yes', 'no')) %>%
  ggplot(aes(Projection_CellType, pct, fill = ever_smoker )) + geom_boxplot() + theme_bw()


ggsave(p.ea, filename = "figures/COPD.pct.smoking.pdf", width = 8, height = 4, device = "pdf")
