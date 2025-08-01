library(Seurat)

# spatialTrans_data_prefix <- "/parallel_scratch/hh01116/Visium_HD_Expr1/"
# singleCell_data_prefix <- "/parallel_scratch/hh01116/raw_data/sc_refer_data/"
spatialTrans_data_prefix <- "/app/data/"


for(batch in c("SAL", "KD", "RAD", "KD_RAD")){
    obj <- readRDS(paste0(spatialTrans_data_prefix, "rds_data/",  batch, "_raw.rds"))
  GetTissueCoordinates(obj) |> 
  saveRDS(paste0(spatialTrans_data_prefix, "rds_data/",  batch, "_st_location.rds"))
    obj[['Spatial.008um']]$counts |>
  saveRDS(st_count, paste0(spatialTrans_data_prefix, "rds_data/", batch, "_st_count.rds"))
}