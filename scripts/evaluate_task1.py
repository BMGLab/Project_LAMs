import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import os

ids = {'spaceranger110_count_38385_OTAR_LNGsp10206160_GRCh38-2020-A': 'P17_T2',
       'spaceranger110_count_36209_OTAR_LNGsp9476038_GRCh38-2020-A': 'P10_T1',
       'spaceranger110_count_39586_OTAR_LNGsp10391238_GRCh38-2020-A': 'P24_T2',
       'spaceranger110_count_40612_OTAR_LNGsp10782313_GRCh38-2020-A': 'P19_B2',
       'spaceranger110_count_39586_OTAR_LNGsp10391237_GRCh38-2020-A': 'P24_T1',
       'spaceranger110_count_36209_OTAR_LNGsp9476039_GRCh38-2020-A': 'P10_T2',
       'spaceranger110_count_38262_OTAR_LNGsp10206166_GRCh38-2020-A': 'P16_T2',
       'spaceranger110_count_38385_OTAR_LNGsp10206159_GRCh38-2020-A': 'P17_T1',
       'spaceranger110_count_36209_OTAR_LNGsp9476040_GRCh38-2020-A': 'P10_T3',
       'spaceranger110_count_38385_OTAR_LNGsp10206158_GRCh38-2020-A': 'P15_T2',
       'spaceranger110_count_38385_OTAR_LNGsp10206157_GRCh38-2020-A': 'P15_T1',
       'spaceranger110_count_38384_OTAR_LNGsp10206164_GRCh38-2020-A': 'D1_2',
       'spaceranger110_count_39586_OTAR_LNGsp10391236_GRCh38-2020-A': 'P25_T2',
       'spaceranger110_count_36209_OTAR_LNGsp9476041_GRCh38-2020-A': 'P10_T4',
       'spaceranger110_count_38384_OTAR_LNGsp10206161_GRCh38-2020-A': 'D2_1',
       'spaceranger110_count_40951_OTAR_LNGsp10922366_GRCh38-2020-A': 'P17_B2',
       'spaceranger110_count_40951_OTAR_LNGsp10922365_GRCh38-2020-A': 'P17_B1',
       'spaceranger110_count_36210_OTAR_LNGsp9476045_GRCh38-2020-A': 'P11_T4',
       'spaceranger110_count_39586_OTAR_LNGsp10391235_GRCh38-2020-A': 'P25_T1',
       'spaceranger110_count_40751_OTAR_LNGsp10782317_GRCh38-2020-A': 'P10_B2',
       'spaceranger110_count_40612_OTAR_LNGsp10782311_GRCh38-2020-A': 'P15_B2',
       'spaceranger110_count_40952_OTAR_LNGsp10922370_GRCh38-2020-A': 'P25_B2',
       'spaceranger110_count_40612_OTAR_LNGsp10782312_GRCh38-2020-A': 'P19_B1',
       'spaceranger110_count_40951_OTAR_LNGsp10922368_GRCh38-2020-A': 'P16_B2',
       'spaceranger110_count_36210_OTAR_LNGsp9476044_GRCh38-2020-A': 'P11_T3',
       'spaceranger110_count_40751_OTAR_LNGsp10782314_GRCh38-2020-A': 'P11_B1',
       'spaceranger110_count_40751_OTAR_LNGsp10782316_GRCh38-2020-A': 'P10_B1',
       'spaceranger110_count_40952_OTAR_LNGsp10922371_GRCh38-2020-A': 'P24_B1',
       'spaceranger110_count_40612_OTAR_LNGsp10782310_GRCh38-2020-A': 'P15_B1',
       'spaceranger110_count_38262_OTAR_LNGsp10206168_GRCh38-2020-A': 'P19_T2',
       'spaceranger110_count_38262_OTAR_LNGsp10206167_GRCh38-2020-A': 'P19_T1',
       'spaceranger110_count_40951_OTAR_LNGsp10922367_GRCh38-2020-A': 'P16_B1',
       'spaceranger110_count_38384_OTAR_LNGsp10206163_GRCh38-2020-A': 'D1_1',
       'spaceranger110_count_38384_OTAR_LNGsp10206162_GRCh38-2020-A': 'D2_2',
       'spaceranger110_count_40952_OTAR_LNGsp10922369_GRCh38-2020-A': 'P25_B1',
       'spaceranger110_count_38262_OTAR_LNGsp10206165_GRCh38-2020-A': 'P16_T1',
       'spaceranger110_count_36210_OTAR_LNGsp9476042_GRCh38-2020-A': 'P11_T1',
       'spaceranger110_count_36210_OTAR_LNGsp9476043_GRCh38-2020-A': 'P11_T2',
       'spaceranger110_count_40751_OTAR_LNGsp10782315_GRCh38-2020-A': 'P11_B2',
       'spaceranger110_count_40952_OTAR_LNGsp10922372_GRCh38-2020-A': 'P24_B2'}

env = "cell2location_task1/spatial_model_SC_Y_sig_matrix_adeno_LUSC_patients"  # Healthy kısmı kaldırıldı
folder = "/home/biolab/Projects/LAMs_wd/results/"
cell2location_res_all_cts = sc.read_h5ad(f"{folder}{env}/sp_our_sig_matrix_all_cts_Healthy.h5ad")
patients = list(cell2location_res_all_cts.uns["spatial"].keys())

#dict to sort the output for each patient
slices = {}

for patient in patients:
    slices[patient] = {}
    #get specific spots for that patient
    spots = [s for s in cell2location_res_all_cts.obs_names if patient in s]
    
    slice_res = cell2location_res_all_cts[spots].copy()
    #slice_res.obs_names = [s.replace(f"{patient}_","") for s in slice_res.obs_names]
    
    #retrieve the cell abundance matrix
    weights = slice_res.obsm["q05_cell_abundance_w_sf"]
    weights.columns = [c.replace("q05cell_abundance_w_sf_", "") for c in weights]

    slices[patient]["w"] = weights.copy()
    #add a column that contains the cell type id with the highest weight
    slices[patient]["w"]['max_w'] = slices[patient]["w"].idxmax(axis=1)
    
    #add two columns that consist of the spot's spatial data
    slices[patient]["w"]["x"] = slice_res.obsm["spatial"][:,0]
    slices[patient]["w"]["y"] = slice_res.obsm["spatial"][:,1]

macro_spots = {}
spots = set()
for patient in patients:
    macro_spots[patient] = slices[patient]["w"][slices[patient]["w"]["max_w"].isin(["Macrophage_alveolar", "Macrophage"])]
    print(ids[patient], ":", macro_spots[patient].shape[0])
    spots = spots.union(set(macro_spots[patient].index))

len(spots)
env_cleaned = env.strip("/").replace("/", "_")  # Replace slashes with underscores

# Now write to a file with a valid name
with open(f"luad_macrophage_spots_{env_cleaned}.txt", "w") as o:
    for spot in spots:
        o.write(f"{spot}\n")

if not os.path.exists(f"output/{env}/deconv_1/"):
    os.makedirs(f"output/{env}/deconv_1/")

for patient in patients:
    slices[patient]["w"].to_csv(f"output/{env}/deconv_1/{ids[patient]}_deconv_1_cell_abundance.xlsx")
