import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import os
import cell2location
import scvi
import torch
import matplotlib.pyplot as plt
import warnings
from visium_qc_and_visualisation import read_and_qc
from config import lung_config

warnings.filterwarnings('ignore')
data_type = 'float32'

def run_cell2loc(ifolder, environment, o_sig_matrix, tag):
    seed = 42
    torch.manual_seed(seed)
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

    if environment == "Tumour":
        if tag == "sq":
            tumour_type = "LUSC"
        elif tag == "adeno":
            tumour_type = "LUAD"
        else:
            raise ValueError(f"Unknown tag for Tumour environment: {tag}")
        outDir = f"/home/biolab/Projects/LAMs_wd/results/cell2location_task1/spatial_model_Tumour_{tumour_type}/"
        sample_names = lung_config[environment]["spatial"][tumour_type]
    else:
        outDir = f"/home/biolab/Projects/LAMs_wd/results/cell2location_task1/spatial_model_{environment}/"
        sample_names = lung_config[environment]["spatial"][environment]

    regDir = f"/home/biolab/Projects/LAMs_wd/data/E-MTAB-13530/{ifolder}/cell2location/regression_model/"
    os.makedirs(outDir, exist_ok=True)

    inf_aver = o_sig_matrix.copy()

    print(f"Running environment: {environment}, samples: {sample_names}")

    inDir = f"/home/biolab/Projects/LAMs_wd/data/E-MTAB-13530/{ifolder}/"
    slides = [read_and_qc(sample_name, path=inDir) for sample_name in sample_names]

    # Clean .obsm if necessary
    for ad in slides:
        for k in ad.obsm.keys():
            if isinstance(ad.obsm[k], pd.DataFrame):
                ad.obsm[k] = ad.obsm[k].to_numpy()

    adata_vis = slides[0].concatenate(slides[1:], batch_key="sample", uns_merge="unique",
                                      batch_categories=sample_names, index_unique=None)

    adata_vis.var["MT_gene"] = [gene.startswith("MT-") for gene in adata_vis.var["SYMBOL"]]
    adata_vis.obsm["MT"] = adata_vis[:, adata_vis.var["MT_gene"].values].X.toarray()
    adata_vis = adata_vis[:, ~adata_vis.var["MT_gene"].values]
    adata_vis.var.set_index("SYMBOL", inplace=True)
    adata_vis.var_names_make_unique()

    intersect = np.intersect1d(np.array(adata_vis.var_names, dtype=str),
                               np.array(inf_aver.index, dtype=str))
    adata_vis = adata_vis[:, adata_vis.var.index.isin(intersect)].copy()
    inf_aver = inf_aver.loc[intersect, :].reindex(adata_vis.var.index)

    cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="sample")
    mod = cell2location.models.Cell2location(adata_vis,
                                             cell_state_df=inf_aver,
                                             detection_alpha=200,
                                             N_cells_per_location=8)

    mod.train(batch_size=2500, train_size=1, max_epochs=3000)
    mod.save(outDir, overwrite=True)
    adata_vis = mod.export_posterior(adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': 2500})
    adata_vis.write(f"{outDir}/sp_our_sig_matrix_all_cts_{environment.replace(' ', '_')}.h5ad")

    train_performance = mod.history["elbo_train"]
    train_performance.plot(logy=True)
    plt.xlabel("Training epochs")
    plt.ylabel("-ELBO loss")
    plt.title(f"Lung ({environment})")
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(f"{outDir}/ELBO_Cell2location_{environment.replace(' ', '_')}.png")
    plt.close()

if __name__ == "__main__":
    o_sig_matrix_sq = pd.read_csv("/home/biolab/Projects/LAMs_wd/data/signature_matrix_sq_24_cell_type_renamed.csv", index_col=0)
    o_sig_matrix_adeno = pd.read_csv("/home/biolab/Projects/LAMs_wd/data/signature_adeno_matrix_24_cell_type_renamed.csv", index_col=0)

    environments = ["Healthy", "Tumour", "Background"]
    for env in environments:
        run_cell2loc("dataset", env, o_sig_matrix_sq, tag="sq")
        run_cell2loc("dataset", env, o_sig_matrix_adeno, tag="adeno")