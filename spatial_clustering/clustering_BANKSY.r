# ================================================
# Identifying spatially-defined tissue domains using BANKSY
# https://benchmarkst-reproducibility.readthedocs.io/en/latest/BANKSY.html
# ================================================

source("util_headers.r")
source(paste0(project_dir, "cluster_labelling/cluster_labelling_utils.r"))

merged_obj <- readRDS(paste0(rds_data_dir, output.file.prefix, "_merged.rds"))

merged_obj_diet <- DietSeurat(merged_obj,
    assays = ifelse(VisiumHD, "sketch", "Spatial"),
    dimreducs = NULL,
    graphs = NULL
)

# extract spatial coordinates
coords <- dplyr::bind_rows(lapply(
    Images(merged_obj_diet),
    function(x) GetTissueCoordinates(merged_obj_diet, image = x, scale = "hires")
))

if (!is.null(coords) && "x" %in% colnames(coords) && "y" %in% colnames(coords)) {
    coords_aligned <- data.frame(row.names = Cells(merged_obj_diet))
    coords_aligned$spatial_x <- coords[rownames(coords_aligned), "x"]
    coords_aligned$spatial_y <- coords[rownames(coords_aligned), "y"]

    if (any(is.na(coords_aligned$spatial_x)) || any(is.na(coords_aligned$spatial_y))) {
        warning("NA values present in spatial_x or spatial_y after alignment. Check image coordinates and cell names.")
    }
    merged_obj_diet$spatial_x <- coords_aligned$spatial_x
    merged_obj_diet$spatial_y <- coords_aligned$spatial_y


    # Add batch information if available in original data
    if ("batch" %in% colnames(merged_obj@meta.data)) {
        if (!("batch" %in% colnames(merged_obj_diet@meta.data))) {
            merged_obj_diet$batch <- merged_obj@meta.data$batch
        } else {
            cat("Batch information already exists in the merged object, skipping addition.\n")
        }
    } else {
        cat("No batch information found in the original merged object.\n")
    }

    # Assign images after subsetting to common_cells
    if (exists("VisiumHD") && VisiumHD) {
        images_filtered <- merged_obj@images[grepl("008um", names(merged_obj@images))]
        names(images_filtered) <- gsub("_slice\\.008um$", "", names(images_filtered))
    } else {
        images_filtered <- merged_obj@images
    }
    # Filter each image to only include common_cells
    for (img_name in names(images_filtered)) {
        img <- images_filtered[[img_name]]
        img_cells <- intersect(Cells(img), Cells(merged_obj_diet))
        images_filtered[[img_name]] <- subset(img, cells = img_cells)
    }
    merged_obj_diet@images <- images_filtered
} else {
    stop("Failed to extract 'x' and 'y' from spatial coordinates. Check image data in merged_obj_diet.")
}

rm(coords, coords_aligned, merged_obj)

# ==========================================
# spatial domain clustering using BANKSY
# https://satijalab.org/seurat/articles/visiumhd_analysis_vignette#identifying-spatially-defined-tissue-domains
# Building Aggregates with a Neighborhood Kernel and Spatial Yardstick (BANKSY), Singhal et al.,
# augments a spot's expression pattern with both
# the mean and the gradient of gene expression levels
# in a spot's broader neighborhood.
# ==========================================

# remotes::install_github('satijalab/seurat-wrappers')
# BiocManager::install("Banksy")

library(Banksy)
merged_obj_diet <- FindVariableFeatures(merged_obj_diet, assay = "sketch", selection.method = "vst")
# `lambda` : Influence of the neighborhood. cell-type clustering (low λ, emphasizing cell's own transcriptome)
# and domain segmentation (high λ, emphasizing the neighborhood context)
# `k_geom` : Local neighborhood size. Larger values will yield larger domains
merged_obj_diet <- SeuratWrappers::RunBanksy(merged_obj_diet,
    lambda = 0.35, verbose = TRUE,
    group = "batch",
    assay = ifelse(VisiumHD, "sketch", "Spatial"), slot = "data", features = "variable",
    k_geom = 18,
    dimx = "spatial_x",
    dimy = "spatial_y"
)
DefaultAssay(merged_obj_diet) <- "BANKSY"
merged_obj_diet <- RunPCA(merged_obj_diet,
    assay = "BANKSY",
    reduction.name = "banksy_pca",
    features = rownames(merged_obj_diet),
    npcs = 30
) |>
    FindNeighbors(reduction = "banksy_pca", dims = 1:30) |>
    FindClusters(cluster.name = "banksy_cluster", resolution = 0.3)

# Export BANKSY cluster assignments to CSV
banksy_results <- data.frame(
    Cell = Cells(merged_obj_diet),
    banksy_cluster = merged_obj_diet$banksy_cluster,
    stringsAsFactors = FALSE
)
write.csv(banksy_results,
    file = paste0(output_dirs$clustering, "banksy_cluster.csv"),
    row.names = FALSE
)






# Read the CSV and set 'Cell' as rownames
banksy_clusters_df <- read.csv(paste0(output_dirs$clustering, "banksy_cluster.csv"), stringsAsFactors = FALSE, check.names = FALSE)
rownames(banksy_clusters_df) <- banksy_clusters_df$Cell
banksy_clusters_df$Cell <- NULL

merged_obj_diet$banksy_cluster <- banksy_clusters_df
Idents(merged_obj_diet) <- "banksy_cluster" # Set idents on the full object for plotting

cluster_labelling_env$plot_spatial_clusters_with_cell_percentage(merged_obj_diet,
    group_by = "banksy_cluster", file_name = "AllCells_banksy_clusters", plot_title = "All cells",
    output_dir = output_dirs$clustering_spatialAware, color_map = color_map, total_cells_per_batch = table(merged_obj_diet$batch), batch_names = batch_names
)

merged_obj_diet <- FindVariableFeatures(merged_obj_diet, assay = "BANKSY", selection.method = "vst")
# Run UMAP for visualization
merged_obj_diet <- RunUMAP(merged_obj_diet,
    reduction = "banksy_pca", dims = 1:30,
    reduction.name = "banksy_umap", n.components = 3
)

cluster_labelling_env$plot_umap(merged_obj_diet, "banksy_umap", "banksy_cluster", "AllCells_banksy_clusters", "All cells", output_dirs$clustering_spatialAware, color_map)

cat("=====================================\n")
cat("Done inferring spatial domain clustering using BANKSY.\n")
cat("=====================================\n")
quit("no")
