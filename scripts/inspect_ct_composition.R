library("HieRFIT")
library("Seurat")
library(tidyverse)
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
library("ggdensity")

#Plot the tree topology
refmod <- readRDS("data/Luca_Subtype_HierMod.RDS")
plt <- PlotTopoNodeAcc(refMod = refmod)

# Export to SVG and then save as PDF
plt_svg <- export_svg(plt)
charToRaw(plt_svg) %>% rsvg_pdf("figures/topo.tree.pdf")

#Read macrophages reannotated
macs <- readRDS("data/luca_query_reannotated.rds")


macs@meta.data$tumor_stage <- droplevels(macs@meta.data$tumor_stage)
macs@meta.data$tumor_stage <- factor(macs@meta.data$tumor_stage,
       levels = c("non-cancer", "early", "advanced"), ordered = TRUE)

table(macs@meta.data$Projection_CellType)

#Plot mac composition early vs advanced
p.ea <- macs@meta.data %>% 
  select(sample, Projection_CellType, tumor_stage, uicc_stage, origin, disease) %>% 
  group_by(sample, Projection_CellType) %>% 
  mutate(Cnt=n()) %>% 
  unique() %>% 
  group_by(sample) %>% 
  mutate(pct=100*Cnt/sum(Cnt), Less=ifelse(sum(Cnt)<50, "low", "ok")) %>% 
  filter(Less == "ok") %>% 
  filter(!Projection_CellType %in% c('Int.Node.3', 'Int.Node.4', 'Int.Node.5')) %>%
  filter(tumor_stage %in% c('early', 'advanced', 'non-cancer')) %>%
  ggplot(aes(Projection_CellType, pct, fill = tumor_stage )) + geom_boxplot() + theme_bw()


ggsave(p.ea, filename = "figures/pct.early-advanced.pdf", width = 8, height = 4, device = "pdf")

#Plot mac composition across uicc stages
macs@meta.data %>% 
  select(sample, Projection_CellType, tumor_stage, uicc_stage, origin, disease) %>% 
  group_by(sample, Projection_CellType) %>% 
  mutate(Cnt=n()) %>% 
  unique() %>% 
  group_by(sample) %>% 
  mutate(pct=100*Cnt/sum(Cnt), Less=ifelse(sum(Cnt)<50, "low", "ok")) %>% 
  filter(Less == "ok") %>% 
  filter(!Projection_CellType %in% c('Int.Node.3', 'Int.Node.4', 'Int.Node.5')) %>%
  filter(tumor_stage %in% c('early', 'advanced', 'non-cancer')) %>%
  mutate(So = if_else(uicc_stage == "IV", paste(uicc_stage, origin), uicc_stage)) %>% 
  filter(So != "IV normal_adjacent") -> So.data
  
So.data$So <- factor(So.data$So,
                        levels = c("non-cancer", "I", "II", "III", "III or IV", "IV tumor_primary", "IV tumor_metastasis"), ordered = TRUE)

p.stg <- ggplot(So.data, aes(Projection_CellType, pct, fill = So )) + geom_boxplot() + theme_bw()
ggsave(p.stg, filename = "figures/pct.stages.pdf", width = 10, height = 4, device = "pdf")

##LA vs RTM antagonist
macs@meta.data %>% 
  select(sample, Projection_CellType, tumor_stage, uicc_stage, origin, disease) %>% 
  group_by(sample, Projection_CellType) %>% 
  mutate(Cnt=n()) %>% 
  unique() %>% 
  group_by(sample) %>% 
  mutate(pct=100*Cnt/sum(Cnt), Less=ifelse(sum(Cnt)<50, "low", "ok")) %>% 
  filter(Less == "ok") %>% 
  filter(Projection_CellType %in% c("LA_TAMs", "RTM_TAMs")) %>% 
  select(Projection_CellType, pct, uicc_stage) %>% 
  reshape2::dcast( formula = uicc_stage + sample ~ Projection_CellType, value.var = "pct") %>% 
  drop_na() -> df3 
df3$uicc_stage <- factor(df3$uicc_stage, levels= c("non-cancer", "I", "II", "III", "III or IV", "IV"), ordered = TRUE)

p.dens <- ggplot(df3, aes(LA_TAMs, RTM_TAMs, color=uicc_stage))+
  geom_hdr_lines(xlim = c(0, 100), ylim = c(0, 100)) +
  geom_point()+
  facet_wrap(vars(uicc_stage)) + theme_bw()

ggsave(p.dens, filename = "figures/pct.density-LAvsRTM.pdf", width = 10, height = 4, device = "pdf")

