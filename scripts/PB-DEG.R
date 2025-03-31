library("HieRFIT")
library("Seurat")
library(tidyverse)
library(dplyr)

library(org.Hs.eg.db)
library(AnnotationDbi)

#load re-annotated macs data
macs <- readRDS("data/luca_query_reannotated.rds")

macs@meta.data$tumor_stage <- droplevels(macs@meta.data$tumor_stage)
macs@meta.data$tumor_stage <- factor(macs@meta.data$tumor_stage,
                                     levels = c("non-cancer", "early", "advanced"), ordered = TRUE)

#Pseudobulk
pseudobulk <- AggregateExpression(macs, assays = "RNA", return.seurat = T, group.by = c("sample", "Projection_CellType"))
pseudobulk

Idents(pseudobulk) <- "Projection_CellType"


bulk.all.7macs.de <- FindAllMarkers(object = pseudobulk,
                             test.use = "MAST", latent.vars = c("sample"))

saveRDS(bulk.all.7macs.de, "data/bulk.all.7macs.de.rds")
#bulk.all.7macs.de <- readRDS("data/bulk.all.7macs.de.rds")

bulk.all.7macs.de %>% 
    filter(pct.1 > 0.5) %>% 
    arrange(desc(avg_log2FC)) %>% 
    filter(avg_log2FC > 1.0) %>% 
    filter(p_val_adj < 0.05) -> All.up.markers.ps


All.up.markers.ps$GeneSymbol <- mapIds(org.Hs.eg.db, keys = All.up.markers.ps$gene, column = "SYMBOL", keytype = "ENSEMBL")
All.up.markers.ps

saveRDS(All.up.markers.ps, "data/PS.up.7macs.deg.rds")


### Explore the PB-based DEG results

pseudobulk@meta.data %>% 
  filter(Projection_CellType %in% c("LA-TAMs", "Angio-TAMs", "Inflam-TAMs", "Prolif-TAMs", "Reg-TAMs", "RTM-TAMs", "IFN-TAMs") ) %>%
  rownames() -> cells.7macs.class

macs_7_types <- c("LA-TAMs", "Angio-TAMs", "Inflam-TAMs", "Prolif-TAMs", "Reg-TAMs", "RTM-TAMs", "IFN-TAMs")


pseudobulk@meta.data$SampleID <- rownames(pseudobulk@meta.data)
custom_colors <- colorRampPalette(c("navy", "white", "red"))(50)


pseudobulk@assays$RNA$scale.data[unique(All.up.markers.ps$gene), cells.7macs.class] %>% 
  reshape2::melt() %>% 
  dplyr::rename(gene = Var1, SampleID = Var2 ) %>% 
  left_join(., unique(All.up.markers.ps[, c("GeneSymbol", "gene")]), by=c("gene")) %>% 
  left_join(., pseudobulk@meta.data, by=c("SampleID")) %>% 
  group_by(GeneSymbol, Projection_CellType) %>% 
  summarize(mean_value = mean(value, na.rm = TRUE) ) %>% 
  dplyr::select(mean_value) %>% pull -> all.data

#print the heatmaps to separate pdf files

for(ct in macs_7_types){

print(ct)

pseudobulk@assays$RNA$scale.data[unique(All.up.markers.ps$gene), cells.7macs.class] %>% 
  reshape2::melt() %>% 
  dplyr::rename(gene = Var1, SampleID = Var2 ) %>% 
  left_join(., unique(All.up.markers.ps[which(All.up.markers.ps$cluster == ct), c("GeneSymbol", "gene")]), by=c("gene")) %>% 
  left_join(., pseudobulk@meta.data, by=c("SampleID")) %>%
  group_by(GeneSymbol, Projection_CellType) %>%
  summarize(mean_value = mean(value, na.rm = TRUE) ) %>%
  filter(!is.na(GeneSymbol)) %>%
  ggplot(aes(x = Projection_CellType, y = factor(GeneSymbol, levels = rev(unique(GeneSymbol))), fill = mean_value)) +
  geom_tile() +
  scale_fill_gradientn(colors = custom_colors, na.value = "grey50") + # Use the custom color palette. and range limits= fixed_range,
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotate x-axis labels
  labs(
    x = "Projection Cell Type",
    y = "Gene Symbol",
    fill = "Mean Value",
    title = "Heatmap of Mean Values"
  ) -> pb.all.ups

h.size <- ceiling(length(All.up.markers.ps[which(All.up.markers.ps$cluster == ct),"GeneSymbol"])/5)

ggsave(plot = pb.all.ups, filename = paste("data/", ct,".genes.pb.results.heatmap.pdf", sep = ""),width = 4, height = h.size, limitsize = FALSE)

}

