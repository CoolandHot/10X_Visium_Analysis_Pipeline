# https://yma-lab.github.io/CARD/documentation/04_CARD_Example.html

library(Seurat)
library(SeuratObject)
library(CARD)

# spatialTrans_data_prefix <- "/parallel_scratch/hh01116/Visium_HD_Expr1/"
# singleCell_data_prefix <- "/parallel_scratch/hh01116/raw_data/sc_refer_data/"
spatialTrans_data_prefix <- "/app/data/"
singleCell_data_prefix <- "/app/sc_refer_data/"
output.dir <- paste0(spatialTrans_data_prefix, "deconvoluted/")
if (!dir.exists(output.dir)) {
  dir.create(output.dir, recursive = TRUE)
}

###################################################################
# load single cell data as cell-type profile reference
# then compute and deconvolute, output as RDS
###################################################################
sc_meta <- readRDS(paste0(singleCell_data_prefix, "meta_data.rds"))
sc_count <- readRDS(paste0(singleCell_data_prefix, "raw_counts_sparse.rds"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  cat("No arguments provided.\n")
  quit(status = 1)
} else {
  batch <- args[1]
}

for(batch in c("SAL", "KD", "RAD", "KD_RAD")){
  st_location <- readRDS(paste0(spatialTrans_data_prefix, "rds_data/",  batch, "_st_location.rds"))
  st_count <- readRDS(paste0(spatialTrans_data_prefix, "rds_data/", batch, "_st_count.rds"))

  CARD_obj <- CARD::createCARDObject(
      sc_count = sc_count,
      sc_meta = sc_meta,
      spatial_count = st_count,
      spatial_location = st_location,
      ct.varname = "cell_type",
      ct.select = unique(sc_meta$cell_type),
      sample.varname = "sample",
      minCountGene = 500,
      minCountSpot = 10
  ) |> CARD::CARD_deconvolution()

  saveRDS(CARD_obj, paste0(output.dir, batch, "_CARD_minGene500_minSpot10.Rds"))
}


# print(CARD_obj@Proportion_CARD[1:2, ])
