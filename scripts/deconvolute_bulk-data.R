library(omnideconv)
library(yaml)
# Load the config
config <- yaml::read_yaml("config.yaml")
# Access credentials
username <- config$cibersortx$username
token <- config$cibersortx$token
# Set credentials
omnideconv::set_cibersortx_credentials(username, token)


signature_adeno <- readRDS("data/lung_adenocarcinoma_cibersortx_sig_mtrx.rds")
adeno_bulk <- readRDS("data/luad_bulk.rds")
deconvolution_adeno <- omnideconv::deconvolute(adeno_bulk, signature_adeno, method = "cibersortx")
saveRDS(deconvolution_adeno, "data/14_02_2025_deconvolution_adeno.rds")


sq_bulk <- readRDS("/home/biolab/OzlemTuna/lusc_bulk.rds")
signature_sq <- readRDS("data/squamous_cell_lung_carcinoma_cibersortx_sig_mtrx.rds")
deconvolution_sq <- omnideconv::deconvolute(sq_bulk, signature_sq, method = "cibersortx")
saveRDS(deconvolution_sq, "/home/biolab/OzlemTuna/14_02_2025_deconvolution_sq.rds")
