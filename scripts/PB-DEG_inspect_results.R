library("HieRFIT")
library("Seurat")
library(tidyverse)
library(dplyr)

library(org.Hs.eg.db)
library(AnnotationDbi)
library("zellkonverter")

#load re-annotated macs data
macs <- readRDS("data/luca_query_reannotated.rds")

macs@meta.data$tumor_stage <- droplevels(macs@meta.data$tumor_stage)
macs@meta.data$tumor_stage <- factor(macs@meta.data$tumor_stage,
                                     levels = c("non-cancer", "early", "advanced"), ordered = TRUE)

#load Pseudobulk h5ad object.
pseudobulk <- zellkonverter::readH5AD("data/ps_adata_macs.h5ad")
pseudobulk
#convert to seurat object
ps_macs_seu <- as.Seurat(pseudobulk, counts = "X", data = NULL)
ps_macs_seu@meta.data |> head()
Idents(ps_macs_seu) <- "Projection_CellType"
ps_macs_seu <- ScaleData(ps_macs_seu)
#load DEGs
ps_macs_deg <- read.csv("results/macs_PB-DEGs.csv")

table(ps_macs_deg$gene_name)
# bulk.all.7macs.de %>% 
#     filter(pct.1 > 0.5) %>% 
#     arrange(desc(avg_log2FC)) %>% 
#     filter(avg_log2FC > 1.0) %>% 
#     filter(p_val_adj < 0.05) -> All.up.markers.ps


#All.up.markers.ps$GeneSymbol <- mapIds(org.Hs.eg.db, keys = All.up.markers.ps$gene, column = "SYMBOL", keytype = "ENSEMBL")
#All.up.markers.ps

#saveRDS(All.up.markers.ps, "data/PS.up.7macs.deg.rds")


### Explore the PB-based DEG results

ps_macs_seu@meta.data %>% 
  filter(Projection_CellType %in% c("LA_TAMs", "Angio_TAMs", "Inflam_TAMs", "Prolif_TAMs", "Reg_TAMs", "RTM_TAMs", "IFN_TAMs") ) %>%
  rownames() -> cells.7macs.class

macs_7_types <- c("LA_TAMs", "Angio_TAMs", "Inflam_TAMs", "Prolif_TAMs", "Reg_TAMs", "RTM_TAMs", "IFN_TAMs")


ps_macs_seu@meta.data$SampleID <- rownames(ps_macs_seu@meta.data)
custom_colors <- colorRampPalette(c("navy", "white", "red"))(50)

#print the heatmaps to separate pdf files
for (ct in macs_7_types){

print(ct)

# Step 1: Compute mean_value for each gene and Projection_CellType
all_means <- ps_macs_seu@assays$originalexp$scale.data[unique(ps_macs_deg$gene_symbols), cells.7macs.class] %>% 
  reshape2::melt() %>% 
  dplyr::rename(gene_symbols = Var1, SampleID = Var2 ) %>% 
  left_join(., unique(ps_macs_deg[which(ps_macs_deg$SubType == ct), c("gene_name", "gene_symbols")]), by = "gene_symbols") %>% 
  left_join(., ps_macs_seu@meta.data, by = "SampleID") %>%
  group_by(gene_name, gene_symbols, Projection_CellType) %>%
  summarize(mean_value = mean(value, na.rm = TRUE), .groups = "drop")

# Step 2: Get max Projection_CellType and its value for each gene
max_celltype <- all_means %>%
  group_by(gene_symbols) %>%
  slice_max(order_by = mean_value, n = 1, with_ties = FALSE) %>%
  dplyr::select(gene_symbols, max_celltype = Projection_CellType, max_value = mean_value)

# Step 3: Join, filter, and plot heatmap
all_means %>%
  left_join(max_celltype, by = "gene_symbols") %>%
  mutate(gene_name = ifelse(ct == max_celltype, gene_name, NA)) %>%
  filter(!is.na(gene_name) & max_value >= 0.25) %>%
  ggplot(aes(x = Projection_CellType, y = factor(gene_name, levels = rev(unique(gene_name))), fill = mean_value)) +
  geom_tile() +
  scale_fill_gradientn(colors = custom_colors, na.value = "grey50") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Projection Cell Type",
    y = "Gene Symbol",
    fill = "Mean Value",
    title = "Heatmap of Mean Values"
  ) -> pb.all.ups



plotted_genes <- all_means %>%
  left_join(max_celltype, by = "gene_symbols") %>%
  mutate(gene_name = ifelse(ct == max_celltype, gene_name, NA)) %>%
  filter(!is.na(gene_name) & max_value >= 0.15) %>%
  distinct(gene_name) %>%
  pull(gene_name)

if(length(plotted_genes) > 5){
  h.size <- ceiling(length(plotted_genes) / 5) # Adjust divisor to control vertical size
}else{
  h.size <- 3
}

ggsave(
  plot = pb.all.ups,
  filename = paste0("figures/", ct, ".genes.pb.results.heatmap.pdf"),
  width = 4,
  height = h.size,
  limitsize = FALSE
)
}
