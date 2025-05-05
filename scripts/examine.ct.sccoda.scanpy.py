import scanpy as sc
import pandas as pd
import random
import numpy as np
import matplotlib as plt
import gc 
import ctypes
import scvi
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import gc 
import ctypes
import sccoda
import pertpy as pt

# Setup
import importlib
import warnings
warnings.filterwarnings("ignore")
import pickle as pkl
import matplotlib.pyplot as plt

from sccoda.util import comp_ana as mod
from sccoda.util import cell_composition_data as dat
from sccoda.util import data_visualization as viz

import sccoda.datasets as scd

adata = sc.read_h5ad('data/luca_query_reannotated.h5ad')

df = adata.obs[['sample', 'Projection_CellType']]

cell_counts = df.groupby(['sample', 'Projection_CellType']).size().reset_index(name='count')

count_matrix = cell_counts.pivot(index='sample', columns='Projection_CellType', values='count').fillna(0)

count_matrix = count_matrix.drop(columns=['Int.Node.3', 'Int.Node.4', 'Int.Node.5'])

print(count_matrix)

row_sums = count_matrix.sum(axis=1)

rows_below_50 = count_matrix.index[row_sums < 50]

df_cleaned = count_matrix.drop(index=rows_below_50)

print(df_cleaned)

df_cleaned['sample'] = df_cleaned.index

data_all = dat.from_pandas(df_cleaned, covariate_columns=["sample"])
data_all.obs.index.name = None

# 1. Build a mapping from sample to tumor_stage
sample_to_stage = adata.obs[['sample', 'tumor_stage']].drop_duplicates().set_index('sample')['tumor_stage']

# 2. Apply the mapping to data_all.obs
data_all.obs['tumor_stage'] = data_all.obs['sample'].map(sample_to_stage)

data_cancer = data_all[data_all.obs["tumor_stage"].isin(["early", "advanced"])]
print(data_cancer.obs)


# Generate the plot
viz.boxplots(data_cancer, feature_name="tumor_stage")

# Save as PDF in the 'figures/' folder
plt.savefig("figures/tumor_stage_boxplot_sccoda.pdf", format="pdf", bbox_inches="tight")

model_cancer = mod.CompositionalAnalysis(data_cancer, formula="tumor_stage", reference_cell_type="Prolif_TAMs")

sim_results = model_cancer.sample_hmc()

sim_results.summary()
print(sim_results.credible_effects())
sim_results.set_fdr(est_fdr=0.4)
sim_results.summary()

# saving
path = "data/sccoda.hmc.results.pkl"
sim_results.save(path)

#sccoda_model.plot_effects_barplot(sccoda_data, "coda", "condition")
#plt.show()