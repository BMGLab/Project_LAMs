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

if not os.path.exists("data/sccoda_data.h5mu"):
    adata = sc.read_h5ad('data/luca_query_reannotated.h5ad')
    print(adata.obs)

    sccoda_model = pt.tl.Sccoda()
    sccoda_data = sccoda_model.load(
        adata,
        type="cell_level",
        generate_sample_level=True,
        cell_type_identifier="Projection_CellType",
        sample_identifier="sample",
        covariate_obs=["tumor_stage"]
    )

    print(sccoda_data)

    os.makedirs("figures", exist_ok=True)

    sccoda_model.plot_boxplots(
        sccoda_data,
        modality_key="coda",
        feature_name="tumor_stage",
        figsize=(12, 5),
        add_dots=False,
        args_swarmplot={},
    )
    plt.savefig("figures/tumor_stage_boxplot_sccoda-1.pdf", format="pdf", bbox_inches="tight")
    plt.close()

    sccoda_model.plot_stacked_barplot(
        sccoda_data,
        modality_key="coda",
        feature_name="tumor_stage",
        figsize=(4, 2)
    )
    plt.savefig("figures/tumor_stage_boxplot_sccoda-2.pdf", format="pdf", bbox_inches="tight")
    plt.close()

    sccoda_data_coda = sccoda_data.mod["coda"]
    print(sccoda_data_coda.obs["tumor_stage"].isna().sum())
    sccoda_data = sccoda_data_coda[~sccoda_data_coda.obs["tumor_stage"].isna()].copy()

    sccoda_data = sccoda_model.prepare(
        sccoda_data,
        modality_key="coda",
        formula="tumor_stage",
        reference_cell_type="automatic"
    )
    sccoda_model.run_nuts(sccoda_data, modality_key="coda", rng_key=1234)
    sccoda_model.set_fdr(sccoda_data, 0.2)
    sccoda_model.credible_effects(sccoda_data, modality_key="coda")

    sccoda_data.write("data/sccoda_data.h5ad")
else:
    sccoda_data = sc.read("data/sccoda_data.h5ad")

# Plot the effects barplot
sccoda_model = pt.tl.Sccoda()
sccoda_model.plot_effects_barplot(sccoda_data, "coda", "tumor_stage")
plt.savefig("figures/tumor_stage_boxplot_sccoda-3.pdf", format="pdf", bbox_inches="tight")
plt.close()
