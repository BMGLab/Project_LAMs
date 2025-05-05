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
print("Loading known marker genes...")
markergenes = pd.read_csv('data/known.markers.txt', sep='\t')

# Get marker gene symbols
marker_gene_symbols = markergenes['GeneSymbol'].dropna().unique().tolist()

# Subset mapping table to only include genes present in adata
available_genes = gtf_df[gtf_df['gene_symbol'].isin(marker_gene_symbols)]
available_genes = available_genes[available_genes['ensembl_id'].isin(adata.var_names)]

# Report which genes are missing
mapped_symbols = available_genes['gene_symbol'].tolist()
missing_genes = set(marker_gene_symbols) - set(mapped_symbols)
print(f"✅ Found {len(mapped_symbols)} marker genes in dataset")
if missing_genes:
    print(f"❌ Missing {len(missing_genes)} genes: {missing_genes}")


# Prepare expression data for plotting
exp_data = []

for _, row in available_genes.iterrows():
    ensembl_id = row['ensembl_id']
    gene_symbol = row['gene_symbol']
    
    expr = adata[:, ensembl_id].X
    expr_array = expr.toarray().flatten() if hasattr(expr, "toarray") else expr.flatten()

    df = pd.DataFrame({
        'Expression': expr_array,
        'CellType': adata.obs['Projection_CellType'].values,
        'Gene': gene_symbol
    })
    exp_data.append(df)

exp_df = pd.concat(exp_data, axis=0)

# --- Step 3: Create and save clustered heatmap of marker gene expression ---

# Average expression per gene per cell type
print("Creating heatmap of average gene expression per Projection_CellType...")

# Pivot exp_df for heatmap
heatmap_df = exp_df.groupby(['CellType', 'Gene'])['Expression'].mean().unstack().fillna(0)

# Standardize expression values per gene (i.e., columnwise scaling)
scaler = StandardScaler()
heatmap_scaled = pd.DataFrame(
    scaler.fit_transform(heatmap_df),
    index=heatmap_df.index,
    columns=heatmap_df.columns
)

# Create clustermap with hierarchical clustering on both axes
g = sns.clustermap(
    heatmap_scaled,
    cmap="vlag",
    figsize=(16, 8),
    linewidths=0.5,
    metric='euclidean',   # Distance metric for clustering
    method='average',     # Linkage method for clustering
    cbar_kws={"label": "Z-scored Expression"},
    row_cluster=True,
    col_cluster=True
)

# Optionally add a title to the cluster map figure
g.fig.suptitle("Clustered Marker Gene Expression (Z-scored) by Projection_CellType", y=1.05)

# Save the clustered heatmap as a PDF in the figures folder
plt.savefig("figures/heatmap_cluster_marker_expression.pdf", bbox_inches='tight')
plt.close()

