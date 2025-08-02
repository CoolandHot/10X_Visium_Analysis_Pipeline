source("util_headers.r")

temp_dir <- paste0(project_dir, "rds_data/temp_", output.file.prefix, "/")
image_name_remove_regex <- ifelse(VisiumHD, "_slice\\.008um$", "_slice$")

ann_output_file <- paste0(project_dir, "rds_data/", output.file.prefix, "_merged.h5ad")

merged_obj <- readRDS(paste0(rds_data_dir, output.file.prefix, "_merged.rds"))


# Create temporary directory for intermediate files
if (!dir.exists(temp_dir)) {
  dir.create(temp_dir, recursive = TRUE)
}

# Remove unwanted metadata columns
merged_obj$seurat_cluster <- NULL
merged_obj$seurat_clusters <- NULL
merged_obj$predict_GBM_scSeq_cellType_1 <- NULL
merged_obj$predict_GBM_scSeq_cellType_2 <- NULL
merged_obj$manual_cell_type <- NULL
merged_obj$clusters_ext <- NULL
merged_obj$nFeature_Spatial.016um <- NULL
merged_obj$nCount_Spatial.016um <- NULL
merged_obj$projected_seurat_cluster.score <- NULL
merged_obj$projected_seurat_cluster <- NULL
merged_obj$sketch_seurat_cluster <- NULL


# --- Keep only the 'Spatial.008um' assay and remove all others ---
merged_obj@assays <- merged_obj@assays[ifelse(VisiumHD, "sketch", "Spatial")]
merged_obj@active.assay <- ifelse(VisiumHD, "sketch", "Spatial")

# Remove all dimensional reductions
merged_obj@reductions <- list()

# --- Keep only spatial images (slices) ending with '.008um' ---
if (VisiumHD) {
  keep_imgs <- grep("\\.008um$", names(merged_obj@images), value = TRUE)
  merged_obj@images <- merged_obj@images[keep_imgs]
} else {
  keep_imgs <- names(merged_obj@images)
}


# Subset columns to only those present in .008um slices
keep_barcodes <- unlist(lapply(keep_imgs, function(img) {
  rownames(Seurat::GetTissueCoordinates(merged_obj, image = img, scale = NULL))
}))
# Intersect with barcodes present in the assay
keep_barcodes <- intersect(keep_barcodes, colnames(merged_obj@assays[[1]]))
merged_obj <- subset(merged_obj, cells = keep_barcodes)

# ========== EXTRACT AND SAVE ALL NECESSARY INFORMATION ==========

# 1. Extract and save spatial coordinates
image_names <- Seurat::Images(merged_obj)
coords_list <- lapply(image_names, function(img) {
  Seurat::GetTissueCoordinates(merged_obj, image = img, scale = NULL)
})
coords_all <- dplyr::bind_rows(coords_list)
coords_all <- coords_all[match(colnames(merged_obj), rownames(coords_all)), c("y", "x")]
coords_mat <- as.matrix(coords_all)
rownames(coords_mat) <- colnames(merged_obj)
saveRDS(coords_mat, paste0(temp_dir, "spatial_coords.rds"))


# 2. Extract and save library IDs
library_ids <- unlist(lapply(image_names, function(img) {
  coords <- Seurat::GetTissueCoordinates(merged_obj, image = img, scale = NULL)
  rep(img, nrow(coords))
}))
library_ids <- library_ids[match(colnames(merged_obj), unlist(lapply(image_names, function(img) {
  rownames(Seurat::GetTissueCoordinates(merged_obj, image = img, scale = NULL))
})))]
library_ids <- gsub(image_name_remove_regex, "", library_ids)
saveRDS(library_ids, paste0(temp_dir, "library_ids.rds"))


# 3. Extract and save high-resolution images
# config$batch_file_names <- paste0("/vol/research/brainTumorST/gbm_data/raw_data/", config$batch_file_names)
hires_images <- lapply(config$batch_file_names, function(batch_file_name) {
  Read10X_Image(
    paste0(batch_file_name, ifelse(VisiumHD, "/outs/binned_outputs/square_008um/spatial", "/spatial")),
    image.name = "tissue_hires_image.png",
    assay = "Spatial",
    slice = batch_file_name,
    filter.matrix = TRUE,
    image.type = "VisiumV2"
  )
})
names(hires_images) <- config$batch_names


# 4. Extract and save spatial metadata for .uns['spatial']
spatial_uns <- list()
for (img in image_names) {
  coords <- Seurat::GetTissueCoordinates(merged_obj, image = img, scale = NULL)
  tissue_positions <- as.matrix(coords[, c("x", "y")])
  rownames(tissue_positions) <- rownames(coords)

  scalefactors <- merged_obj@images[[img]]@scale.factors
  scalefactors[["tissue_hires_scalef"]] <- scalefactors$hires
  scalefactors$hires <- NULL
  scalefactors[["tissue_lowres_scalef"]] <- scalefactors$lowres
  scalefactors$lowres <- NULL
  scalefactors[["spot_diameter_fullres"]] <- scalefactors$spot
  scalefactors$spot <- NULL
  scalefactors[["fiducial_diameter_fullres"]] <- scalefactors$fiducial
  scalefactors$fiducial <- NULL

  lowres_img <- merged_obj@images[[img]]@image


  # Build spatial_uns here
  images_list <- list()
  if (!is.null(lowres_img)) {
    images_list[["lowres"]] <- as.array(lowres_img)
    # Get corresponding hires image
    img_key <- gsub(image_name_remove_regex, "", img)
    if (img_key %in% names(hires_images)) {
      images_list[["hires"]] <- as.array(hires_images[[img_key]]@image)
    }
  }

  spatial_uns[[img]] <- list(
    images = images_list,
    scalefactors = scalefactors,
    metadata = list(
      chemistry_description = "Visium",
      software_version = "Seurat"
    )
  )
}

# Apply consistent naming transformation
names(spatial_uns) <- gsub(image_name_remove_regex, "", names(spatial_uns))

saveRDS(spatial_uns, paste0(temp_dir, "spatial_uns.rds"))

# 5. Extract and save counts matrix in correct format for later use
counts <- Seurat::GetAssayData(merged_obj, layer = "counts")
saveRDS(list(
  genes = rownames(counts),
  cells = colnames(counts)
), paste0(temp_dir, "counts_info.rds"))

# counts <- counts[colnames(ann), ]
# counts <- counts[, rownames(ann)]
# if (!is.matrix(counts)) counts <- as.matrix(counts)
# ann$X <- t(counts)

# ========== CONVERT TO ANNDATA AND CLEAR MEMORY ==========

# Convert to AnnData
# scCustomize::as.anndata(x = merged_obj, file_path = "/scratch/hh01116/download", file_name = "DIPG_B7_expr_six_merged_merged.h5ad")
scCustomize::as.anndata(x = merged_obj, file_path = paste0(project_dir, "rds_data/"), file_name = paste0(output.file.prefix, "_merged.h5ad"))

# Clear merged_obj from memory
rm(merged_obj)
gc() # Force garbage collection

# ========== LOAD ANNDATA AND ADD BACK INFORMATION ==========

# Skip pip upgrade
Sys.setenv(RETICULATE_AUTOCONFIGURE = "FALSE")
ann <- anndata::read_h5ad(ann_output_file)

# Load saved information
coords_mat <- readRDS(paste0(temp_dir, "spatial_coords.rds"))
library_ids <- readRDS(paste0(temp_dir, "library_ids.rds"))
spatial_uns <- readRDS(paste0(temp_dir, "spatial_uns.rds"))
counts_info <- readRDS(paste0(temp_dir, "counts_info.rds"))

# 2. Add spatial coordinates to AnnData .obsm['spatial']
ann$obsm[["spatial"]] <- coords_mat

# 3. Add library_id to AnnData obs
ann$obs[["library_id"]] <- library_ids

# 4. Add .uns['spatial'] for Scanpy compatibility
ann$uns[["spatial"]] <- spatial_uns

# Save the updated AnnData object
ann$write_h5ad(ann_output_file)

# Clean up temporary files
unlink(temp_dir, recursive = TRUE)

cat("Successfully created AnnData file with spatial information:", ann_output_file, "\n")
