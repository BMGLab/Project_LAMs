import warnings
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import seaborn as sns
import os
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")
warnings.simplefilter("ignore")

# Load AnnData
adata = sc.read_h5ad('data/luca_query_reannotated.h5ad')
os.makedirs("figures", exist_ok=True)

# --- Step 0: Filter out unwanted cell types ---
exclude_types = ["Int.Node.3", "Int.Node.4", "Int.Node.5"]
adata = adata[~adata.obs["Projection_CellType"].isin(exclude_types)].copy()
print(f"Filtered out {exclude_types}, remaining cells: {adata.n_obs}")

# --- Step 1: Load and parse the GTF file to get gene_id ↔ gene_name mapping ---
print("Parsing GTF for gene name mapping...")
with open('data/Homo_sapiens.GRCh38.104.gtf') as f:
    gtf = [line for line in f if not line.startswith('#') and 'gene_id "' in line and 'gene_name "' in line]

gtf_mapping = dict(map(lambda x: (
    x.split('gene_id "')[1].split('"')[0],
    x.split('gene_name "')[1].split('"')[0]
), gtf))

# Create mapping DataFrame
gtf_df = pd.DataFrame.from_dict(gtf_mapping, orient='index', columns=['gene_symbol'])
gtf_df.index.name = 'ensembl_id'
gtf_df = gtf_df.reset_index()

# --- Step 2: Load marker genes and assign category ---
print("Loading Martinez marker genes...")
m1_m2_table = pd.read_csv('data/Martinez_etal_M1M2_markers.txt', sep='\t')
m1_m2_table['M1:M2Ratio'] = pd.to_numeric(m1_m2_table['M1:M2Ratio'], errors='coerce')
m1_m2_table['category'] = m1_m2_table['M1:M2Ratio'].apply(lambda x: 'M1' if x > 0 else 'M2')

m1_genes_raw = m1_m2_table[m1_m2_table['category'] == 'M1']['GeneSymbol'].dropna().unique().tolist()
m2_genes_raw = m1_m2_table[m1_m2_table['category'] == 'M2']['GeneSymbol'].dropna().unique().tolist()

# --- Step 3: Map gene symbols to Ensembl IDs using GTF ---
adata_genes = adata.var_names.astype(str)
symbol_to_ens = gtf_df[gtf_df['ensembl_id'].isin(adata_genes)]

m1_genes = symbol_to_ens[symbol_to_ens['gene_symbol'].isin(m1_genes_raw)]['ensembl_id'].tolist()
m2_genes = symbol_to_ens[symbol_to_ens['gene_symbol'].isin(m2_genes_raw)]['ensembl_id'].tolist()

# Report missing genes
missing_m1 = set(m1_genes_raw) - set(symbol_to_ens[symbol_to_ens['ensembl_id'].isin(m1_genes)]['gene_symbol'])
missing_m2 = set(m2_genes_raw) - set(symbol_to_ens[symbol_to_ens['ensembl_id'].isin(m2_genes)]['gene_symbol'])

print(f"✅ M1 genes matched: {len(m1_genes)} | ❌ Missing M1: {len(missing_m1)} → {missing_m1}")
print(f"✅ M2 genes matched: {len(m2_genes)} | ❌ Missing M2: {len(missing_m2)} → {missing_m2}")

if len(m1_genes) == 0:
    raise ValueError("No M1 genes found in adata.var_names after mapping.")
if len(m2_genes) == 0:
    raise ValueError("No M2 genes found in adata.var_names after mapping.")

# --- Step 4: Compute expression scores ---
adata.obs['M1_score'] = adata[:, m1_genes].X.mean(axis=1)
adata.obs['M2_score'] = adata[:, m2_genes].X.mean(axis=1)

# --- Step 5: FACS-style scatter plot ---
plt.figure(figsize=(6, 6))
scatter = sns.scatterplot(
    x=adata.obs['M1_score'].A1 if hasattr(adata.obs['M1_score'], 'A1') else adata.obs['M1_score'],
    y=adata.obs['M2_score'].A1 if hasattr(adata.obs['M2_score'], 'A1') else adata.obs['M2_score'],
    hue=adata.obs['Projection_CellType'],
    alpha=0.5,
    s=2  # marker size in plot
)

# Customize legend marker size
handles, labels = scatter.get_legend_handles_labels()
scatter.legend(
    handles=handles,
    labels=labels,
    title="Macrophage Subtype",
    bbox_to_anchor=(1.05, 1),
    loc='upper left',
    borderaxespad=0.,
    markerscale=.5  # 🔥 Increase this value to enlarge legend dot size
)
plt.xlabel("Mean M1 Gene Expression")
plt.ylabel("Mean M2 Gene Expression")
plt.title("FACS-style Plot: M1 vs M2 Expression")
plt.legend(title="Macrophage type", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True)
plt.tight_layout()
plt.savefig("figures/facs_M1_vs_M2_expression_by_tumor_stage.pdf")
plt.close()

# --- Unified heatmap: column-scaled (z-score) M1 and M2 scores per cell type ---

# Compute raw mean scores
heatmap_joint_df = (
    adata.obs
    .groupby("Projection_CellType")[["M1_score", "M2_score"]]
    .mean()
    .sort_index()
)

# Scale by column (z-score)
scaler = StandardScaler()
scaled_values = scaler.fit_transform(heatmap_joint_df.values)
heatmap_scaled_df = pd.DataFrame(
    scaled_values,
    index=heatmap_joint_df.index,
    columns=heatmap_joint_df.columns
)

# Transpose to reverse rows and columns
heatmap_scaled_df_T = heatmap_scaled_df.T

# Plot the transposed (columns <-> rows) scaled heatmap
plt.figure(figsize=(max(6, 0.5 * heatmap_scaled_df_T.shape[1]), 3))
sns.heatmap(
    heatmap_scaled_df_T,
    annot=True,
    cmap="vlag",
    center=0,
    fmt=".2f",
    cbar_kws={'orientation': 'vertical'}
)
plt.title("Scaled M1 and M2 Expression by Macrophage Subtype")
plt.xlabel("Macrophage Subtype")
plt.ylabel("Gene Signature")
plt.xticks(rotation=90)  # rotate column labels (macrophage subtypes)
plt.yticks(rotation=0)   # vertical row labels (M1/M2)
plt.tight_layout()
plt.savefig("figures/heatmap_zscore_M1_M2_expression_by_celltype_transposed.pdf")
plt.close()

# Save transposed heatmap data
heatmap_scaled_df_T.to_csv("figures/heatmap_zscore_M1_M2_expression_data_transposed.csv")
