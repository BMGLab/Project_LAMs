import warnings
import pandas as pd
warnings.filterwarnings("ignore")
warnings.simplefilter("ignore")

import matplotlib
import matplotlib.pyplot as plt
import os
import numpy as np
import pertpy as pt
import scanpy as sc
import seaborn as sns


adata = sc.read_h5ad('data/luca_query_reannotated.h5ad')

# Select relevant columns
df = adata.obs[['donor_id', 'disease']]

# Group by disease and count unique donors
unique_donors_per_disease = df.groupby('disease')['donor_id'].nunique()

print(unique_donors_per_disease)