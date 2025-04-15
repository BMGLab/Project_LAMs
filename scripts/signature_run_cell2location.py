import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import os
import cell2location
import scvi
import torch
import matplotlib as mpl
from matplotlib import rcParams
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.backends.backend_pdf as mpdf
import warnings
warnings.filterwarnings('ignore')
data_type = 'float32'
from visium_qc_and_visualisation import read_and_qc
from config import lung_config

# Load the data
o_sig_matrix = pd.read_csv("/home/biolab/Projects/LAMs_wd/data/lung_adenocarcinoma_cibersortx_macs_sig_mtrx.csv", index_col=0)

# Run Cell2location
def run_cell2loc(ifolder, environment):
    seed = 42
    torch.manual_seed(seed)
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

    outDir = f"/result/cell2location/spatial_model_[SC_Y_sig_matrix_adeno]_[LUAD_patients]/"
    regDir = f"/result/cell2location/regression_model/"
 
    inf_aver = None
    # Export estimated expression signature in each cluster
    
    inf_aver = o_sig_matrix.copy()

    # Scale up by average sample scaling factor. This corrects for sequencing depth
    # inf_aver = inf_aver * adata_sc_model.uns["mod"]["post_sample_means"]["sample_scaling"].mean() # Need to check if this is still necessary

    # Pick up Visium
    sample_names = lung_config[environment]["spatial"]
    print(sample_names)
    #inDir = "/lustre/scratch118/opentargets/opentargets/OTAR2060/nelson/data/lung/visium/samples/{}/".format(environment.lower())
    inDir = "/data/E-MTAB-13530/"
    slides = [read_and_qc(sample_name, path=inDir) for sample_name in sample_names]
    adata_vis = slides[0].concatenate(slides[1:], batch_key="sample", uns_merge="unique", batch_categories=sample_names, index_unique=None)
    # Mitochondria-encoded (MT) genes should be removed for spatial mapping
    # Find mitochondria-encoded (MT) genes
    adata_vis.var["MT_gene"] = [gene.startswith("MT-") for gene in adata_vis.var["SYMBOL"]]
    # Remove MT genes for spatial mapping (keeping their counts in the object)
    adata_vis.obsm["MT"] = adata_vis[:, adata_vis.var["MT_gene"].values].X.toarray()
    adata_vis = adata_vis[:, ~adata_vis.var["MT_gene"].values]
    # Reset the index for compatibility with the scRNA-seq counts
    adata_vis.var.set_index("SYMBOL", inplace=True)

    # Find shared genes and subset both anndata and reference signatures
    adata_vis.var_names_make_unique()
    
    inf_aver_index = np.array(inf_aver.index, dtype=str)
    adata_vis_var_names = np.array(adata_vis.var_names, dtype=str)
    intersect = np.intersect1d(adata_vis_var_names, inf_aver_index)
    #intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
    #adata_vis = adata_vis[:, intersect] # This doesn't work -- think there's some incompatibility between ScanPy and AnnData versions. Use the line below instead
    adata_vis = adata_vis[:, adata_vis.var.index.isin(intersect)].copy()
    inf_aver = inf_aver.loc[intersect, :]
    # Reindex for cell2location
    inf_aver = inf_aver.reindex(adata_vis.var.index)

    # Save on memory
#     del adata_sc_model
#     del mod

    # Training cell2location
    cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="sample")
    mod = cell2location.models.Cell2location(adata_vis,
                                             cell_state_df=inf_aver,
                                             detection_alpha=200, # Controls the normalisation of the within-experiment variation of the RNA detection (default value for now)
                                             N_cells_per_location=8
                                             )
    '''
                                             **model_kwargs={
                                                # Prior on the number of cells, cell types and co-located groups. Hyperparameter inputs go here
                                                "cell_number_prior": {
                                                    # - N - the expected number of cells per location:
                                                    "cells_per_spot": 8,
                                                    # - A - the expected number of cell types per location:
                                                    "factors_per_spot": 9,
                                                    # - Y - the expected number of co-located cell type groups per location
                                                    "combs_per_spot": 5
                                                },
                                                # Prior beliefs on the sensitivity of spatial technology:
                                                "gene_level_prior": {
                                                    # Prior on the mean
                                                    "mean": 0.5,
                                                    # Prior on standard deviation
                                                    # A good choice of this value should be at least 2 times lower that the mean
                                                    "sd": 0.15
                                                }
                                             }
                                             )
    '''

    # Try 2500 batches initially, rather than training on the entire data-set
    mod.train(batch_size=2500, train_size=1, max_epochs=3000)
    mod.save(f"{outDir}{environment}", overwrite=True)
    adata_vis = mod.export_posterior(adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': 2500})
    adata_vis.write(f"{outDir}{environment}/sp_our_sig_matrix_all_cts_{environment}.h5ad")

    train_performance = mod.history["elbo_train"]
    train_performance.plot(logy=True)
    plt.xlabel("Training epochs")
    plt.ylabel("-ELBO loss")
    plt.title(f"Lung ({environment})")
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(f"{outDir}{environment}/ELBO_Cell2location_{environment}.png")
    plt.close()

run_cell2loc("dataset", "Background")
