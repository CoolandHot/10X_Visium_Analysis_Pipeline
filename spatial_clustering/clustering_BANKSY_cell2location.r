# ================================================
# BANKSY clustering on cell2location results
# Using cell abundances from cell2location as features for spatial clustering
# ================================================

source("util_headers.r")
source(paste0("./cluster_labelling/cluster_labelling_utils.r"))

# Read cell2location results
cell2loc_data <- read.csv(paste0(output_dirs$cell_type_cell2loc, "cell_abundances_and_clusters.csv"),
    stringsAsFactors = FALSE, row.names = 1
)

# Check for NA or NaN values in cell2loc_data
na_count <- sum(is.na(cell2loc_data))
null_count <- sum(is.null(cell2loc_data))
nan_count <- sum(is.nan(as.matrix(cell2loc_data)))
if (na_count > 0 || nan_count > 0 || null_count > 0) {
    warning(sprintf("cell2loc_data contains %d NA, %d NaN, and %d NULL values.", na_count, nan_count, null_count))
}

if (na_count > 0 || nan_count > 0 || null_count > 0) {
    warning("Removing rows with NA or NaN values from cell2loc_data.")
    cell2loc_data <- cell2loc_data[complete.cases(cell2loc_data) & !apply(is.nan(as.matrix(cell2loc_data)), 1, any), , drop = FALSE]
}

# Read original merged object to extract spatial coordinates
merged_obj_original <- readRDS(paste0(rds_data_dir, output.file.prefix, "_merged.rds"))

# extract the `batch` column if it exists
if ("batch" %in% colnames(cell2loc_data)) {
    batch_info <- cell2loc_data[, "batch", drop = FALSE]
} else {
    batch_info <- NULL
}
# remove the `batch` column from abundance data if it exists
if (!is.null(batch_info)) {
    cell2loc_data <- cell2loc_data[, !colnames(cell2loc_data) %in% c("batch", "region_cluster"), drop = FALSE]
}



# Create new Seurat object from cell abundance data
merged_obj_diet <- CreateSeuratObject(
    counts = t(as.matrix(cell2loc_data)),
    project = "cell2location_BANKSY",
    assay = "CellAbundance"
)
merged_obj_diet$batch <- batch_info


# Extract spatial coordinates from original object
coords <- dplyr::bind_rows(lapply(
    Images(merged_obj_original),
    function(x) GetTissueCoordinates(merged_obj_original, image = x, scale = "hires")
))

# Align coordinates with cell2location data
if (!is.null(coords) && "x" %in% colnames(coords) && "y" %in% colnames(coords)) {
    # Match cells between datasets
    common_cells <- intersect(rownames(cell2loc_data), rownames(coords))

    if (length(common_cells) == 0) {
        stop("No common cells found between cell2location data and spatial coordinates")
    }

    # Subset objects to common cells
    merged_obj_diet <- subset(merged_obj_diet, cells = common_cells)
    coords_subset <- coords[common_cells, ]

    # Add spatial coordinates to metadata
    merged_obj_diet$spatial_x <- coords_subset$x
    merged_obj_diet$spatial_y <- coords_subset$y

    # Add batch information if available in original data
    if ("batch" %in% colnames(merged_obj_original@meta.data)) {
        batch_info_origin <- merged_obj_original@meta.data[common_cells, "batch", drop = FALSE]
        if (!("batch" %in% colnames(merged_obj_diet@meta.data))) {
            merged_obj_diet$batch <- batch_info_origin$batch
        } else {
            cat("Batch information already exists in the merged object, skipping addition.\n")
        }
    }

    # Assign images after subsetting to common_cells
    if (exists("VisiumHD") && VisiumHD) {
        images_filtered <- merged_obj_original@images[grepl("008um", names(merged_obj_original@images))]
        names(images_filtered) <- gsub("_slice\\.008um$", "", names(images_filtered))
    } else {
        images_filtered <- merged_obj_original@images
    }
    # Filter each image to only include common_cells
    for (img_name in names(images_filtered)) {
        img <- images_filtered[[img_name]]
        img_cells <- intersect(Cells(img), common_cells)
        images_filtered[[img_name]] <- subset(img, cells = img_cells)
    }
    merged_obj_diet@images <- images_filtered

    cat("Successfully aligned", length(common_cells), "cells with spatial coordinates\n")
} else {
    stop("Failed to extract spatial coordinates from original object")
}

# Find variable features (cell types with highest variation)
merged_obj_diet <- FindVariableFeatures(merged_obj_diet, assay = "CellAbundance", selection.method = "vst")

# ==========================================
# BANKSY clustering on cell abundance data
# ==========================================

library(Banksy)

# Set default assay
DefaultAssay(merged_obj_diet) <- "CellAbundance"

# Run BANKSY with cell abundance features
merged_obj_diet <- SeuratWrappers::RunBanksy(merged_obj_diet,
    lambda = 0.35, verbose = TRUE,
    group = "batch",
    assay = "CellAbundance", slot = "counts", features = "variable",
    k_geom = 18,
    dimx = "spatial_x",
    dimy = "spatial_y"
)

# Run PCA on BANKSY features
DefaultAssay(merged_obj_diet) <- "BANKSY"
merged_obj_diet <- RunPCA(merged_obj_diet,
    assay = "BANKSY",
    reduction.name = "banksy_pca",
    features = rownames(merged_obj_diet),
    npcs = 30
) |>
    FindNeighbors(reduction = "banksy_pca", dims = 1:30) |>
    FindClusters(cluster.name = "banksy_cluster_cell2loc", resolution = 0.3)

# Run UMAP for visualization
merged_obj_diet <- RunUMAP(merged_obj_diet,
    reduction = "banksy_pca", dims = 1:30,
    reduction.name = "banksy_umap", n.components = 3
)

# Export BANKSY cluster assignments
banksy_results <- data.frame(
    Cell = Cells(merged_obj_diet),
    banksy_cluster_cell2loc = merged_obj_diet$banksy_cluster_cell2loc,
    stringsAsFactors = FALSE
)

write.csv(banksy_results,
    file = paste0(output_dirs$clustering, "banksy_cluster_cell2location.csv"),
    row.names = FALSE
)

# Create output directory for cell2location BANKSY results
output_dir_cell2loc <- output_dirs$clustering_spatialAware
dir.create(output_dir_cell2loc, recursive = TRUE, showWarnings = FALSE)

# # Generate visualizations
# banksy_results <- read.csv(paste0(output_dirs$clustering, "banksy_cluster_cell2location.csv"),
#     stringsAsFactors = FALSE, row.names = 1
# )
# # Ensure banksy_cluster_cell2loc is a vector, not a data frame
# merged_obj_diet$banksy_cluster_cell2loc <- banksy_results$banksy_cluster_cell2loc
Idents(merged_obj_diet) <- "banksy_cluster_cell2loc"

# Spatial plot
if ("batch" %in% colnames(merged_obj_diet@meta.data)) {
    batch_names <- unique(merged_obj_diet$batch)
    total_cells_per_batch <- table(merged_obj_diet$batch)
} else {
    batch_names <- "All"
    total_cells_per_batch <- length(Cells(merged_obj_diet))
    names(total_cells_per_batch) <- "All"
}

# Ensure color_map is defined and has enough colors for the number of clusters
num_clusters <- length(unique(na.omit(merged_obj_diet$banksy_cluster_cell2loc)))
if (num_clusters == 0) {
    stop("No clusters found in banksy_cluster_cell2loc. Cannot generate color map.")
}
if (!exists("color_map") || length(color_map) < num_clusters) {
    if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
        install.packages("RColorBrewer")
    }
    library(RColorBrewer)
    # Use Set3 or another palette with enough colors
    palette_name <- "Set3"
    max_colors <- brewer.pal.info[palette_name, "maxcolors"]
    if (num_clusters <= max_colors) {
        color_map <- brewer.pal(num_clusters, palette_name)
    } else {
        # If more clusters than palette supports, interpolate colors
        color_map <- colorRampPalette(brewer.pal(max_colors, palette_name))(num_clusters)
    }
}


# Plot spatial clusters
cluster_labelling_env$plot_spatial_clusters_with_cell_percentage(merged_obj_diet,
    group_by = "banksy_cluster_cell2loc",
    file_name = "Cell2Location_banksy_clusters",
    plot_title = "BANKSY clusters (Cell2Location)",
    output_dir = output_dir_cell2loc,
    color_map = color_map,
    total_cells_per_batch = total_cells_per_batch,
    batch_names = batch_names
)

# Plot UMAP
cluster_labelling_env$plot_umap(
    merged_obj_diet, "banksy_umap", "banksy_cluster_cell2loc",
    "Cell2Location_banksy_clusters", "BANKSY clusters (Cell2Location)",
    output_dir_cell2loc, color_map
)

cat("=====================================\n")
cat("Done running BANKSY clustering on cell2location results.\n")
cat("Results saved to:", output_dir_cell2loc, "\n")
cat("Number of clusters found:", length(unique(merged_obj_diet$banksy_cluster_cell2loc)), "\n")
cat("=====================================\n")
