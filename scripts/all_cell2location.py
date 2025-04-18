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
sys.path.append("/home/biolab/Projects/LAMs_wd/scripts")
warnings.filterwarnings('ignore')
data_type = 'float32'
from visium_qc_and_visualisation import read_and_qc
from config import lung_config

o_sig_matrix = pd.read_csv("/home/biolab/Projects/LAMs_wd/data/signature_matrix_sq_24_cell_type_renamed.csv", index_col=0)

def run_cell2loc(ifolder, environment):
    seed = 42
    torch.manual_seed(seed)
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

    outDir = f"/home/biolab/Projects/LAMs_wd/results/cell2location/spatial_model_SC_Y_sig_matrix_adeno_LUSC_patients/{environment}/"
    regDir = f"/home/biolab/Projects/LAMs_wd/data/E-MTAB-13530/{ifolder}/cell2location/regression_model/"
    os.makedirs(outDir, exist_ok=True)

    inf_aver = o_sig_matrix.copy()

    if environment == "Tumour":
        sample_names = lung_config[environment]["spatial"]["LUSC"]
    else:
        sample_names = lung_config[environment]["spatial"]

    print(sample_names)
    inDir = "/home/biolab/Projects/LAMs_wd/data/E-TAMB-13530/"
    slides = [read_and_qc(sample_name, path=inDir) for sample_name in sample_names]
    adata_vis = slides[0].concatenate(slides[1:], batch_key="sample", uns_merge="unique", batch_categories=sample_names, index_unique=None)

    adata_vis.var["MT_gene"] = [gene.startswith("MT-") for gene in adata_vis.var["SYMBOL"]]
    adata_vis.obsm["MT"] = adata_vis[:, adata_vis.var["MT_gene"].values].X.toarray()
    adata_vis = adata_vis[:, ~adata_vis.var["MT_gene"].values]
    adata_vis.var.set_index("SYMBOL", inplace=True)
    adata_vis.var_names_make_unique()

    inf_aver_index = np.array(inf_aver.index, dtype=str)
    adata_vis_var_names = np.array(adata_vis.var_names, dtype=str)
    intersect = np.intersect1d(adata_vis_var_names, inf_aver_index)

    adata_vis = adata_vis[:, adata_vis.var.index.isin(intersect)].copy()
    inf_aver = inf_aver.loc[intersect, :]
    inf_aver = inf_aver.reindex(adata_vis.var.index)

    cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="sample")
    mod = cell2location.models.Cell2location(
        adata_vis,
        cell_state_df=inf_aver,
        detection_alpha=200,
        N_cells_per_location=8
    )

    mod.train(batch_size=2500, train_size=1, max_epochs=3000)
    mod.save(outDir, overwrite=True)
    adata_vis = mod.export_posterior(adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': 2500})
    adata_vis.write(f"{outDir}/sp_our_sig_matrix_all_cts_{environment}.h5ad")

    train_performance = mod.history["elbo_train"]
    train_performance.plot(logy=True)
    plt.xlabel("Training epochs")
    plt.ylabel("-ELBO loss")
    plt.title(f"Lung ({environment})")
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(f"{outDir}/ELBO_Cell2location_{environment}.png")
    plt.close()

for environment in ["Tumour", "Background", "Healthy"]:
    run_cell2loc("dataset", environment)