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
sclc <- readRDS("data/20250813_sclc.rds")

#subset the macrophages for SCLC
sclc_macs <- subset(sclc, cell_type %in% c("Macrophage_monocyte"))
#number of donors
sclc_macs@meta.data[,c("patient", "condition")] %>% unique %>% select(patient) %>% table
#number of samples
#luca_new@meta.data[,c("sample", "disease")] %>% unique %>% select(disease) %>% table

#normalize the data
sclc_macs <- NormalizeData(sclc_macs, scale.factor = 10000)


#Load hierfit model

refmod <- readRDS("data/Luca_Subtype_HierMod.RDS")

### Reannotation of Query Data with HieRFIT Model;
#Project the cell class labels on the new dataset:
hierObj <- HieRFIT(Query = sclc_macs[["RNA"]]$data, refMod = refmod)

p.pred <- PlotBarStats(HieRobj = hierObj)

ggsave(p.pred, filename = "figures/sclc.macs.predictions.barplot.pdf", width = 8, height = 4, device = "pdf")


sclc_macs@meta.data <- cbind(sclc_macs@meta.data, hierObj@Evaluation$Projection)
colnames(sclc_macs@meta.data)[colnames(sclc_macs@meta.data) == "hierObj@Evaluation$Projection"] <- "Projection_CellType"


saveRDS(sclc_macs, file = "data/sclc.luca_query_reannotated.rds")

# Load your Seurat object
luca_new <- readRDS("data/sclc.luca_query_reannotated.rds")
sce <- as.SingleCellExperiment(luca_new)
writeH5AD(sce, file = "/data/sclc.luca_query_reannotated.h5ad")


#composition of macrophages in COPD
luca_new <- readRDS("data/sclc.luca_query_reannotated.rds")

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
