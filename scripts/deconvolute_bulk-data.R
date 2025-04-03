library(omnideconv)
library(yaml)
# Load the config
config <- yaml::read_yaml("config.yaml")
# Access credentials
username <- config$cibersortx$username
token <- config$cibersortx$token
# Set credentials
omnideconv::set_cibersortx_credentials(username, token)

print("Loading bulk expression data of TCGA for LUAD and LUSC...")
bulk_ad <- readRDS("data/bulk_data_ad.rds")
bulk_sq <- readRDS("data/bulk_data_sq.rds")

print("Deconvolution of bulk LUAD with adenocarcinoma specific macrophage signature...")
sig_mtx_ad <- readRDS("data/lung_adenocarcinoma_cibersortx_macs_sig_mtrx.rds")
deconved_bulk_ad <- omnideconv::deconvolute(bulk_ad, sig_mtx_ad, method = "cibersortx")
saveRDS(deconved_bulk_ad, "data/macs_deconved_bulk_ad.rds")
print("done.")

print("Deconvolution of bulk LUSC with squamous cell lc specific macrophage signature...")
sig_mtx_sq <- readRDS("data/squamous_cell_lung_carcinoma_macs_cibersortx_sig_mtrx.rds")
deconved_bulk_sq <- omnideconv::deconvolute(bulk_sq, sig_mtx_sq, method = "cibersortx")
saveRDS(deconved_bulk_sq, "data/macs_deconved_bulk_sq.rds")
print("done.")
