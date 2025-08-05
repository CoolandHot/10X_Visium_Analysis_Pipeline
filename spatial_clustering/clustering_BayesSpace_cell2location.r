# ================================================
# BayesSpace clustering on cell2location results
# Using cell abundances from cell2location as features for spatial clustering
# ================================================

source("util_headers.r")
source(paste0("./cluster_labelling/cluster_labelling_utils.r"))

# Read cell2location results
cell2loc_data <- read.csv(paste0(output_dirs$cell_type_cell2loc, "cell_abundances_and_clusters.csv"),
    stringsAsFactors = FALSE, row.names = 1
)

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

    # Subset data to common cells
    cell2loc_data_subset <- cell2loc_data[common_cells, , drop = FALSE]
    coords_subset <- coords[common_cells, ]

    cat("Successfully aligned", length(common_cells), "cells with spatial coordinates\n")
} else {
    stop("Failed to extract spatial coordinates from original object")
}

# ==========================================
# BayesSpace clustering on cell abundance data
# ==========================================
library(SingleCellExperiment)

# Create SingleCellExperiment object directly
sce <- SingleCellExperiment(
    assays = list(counts = t(as.matrix(cell2loc_data_subset))),
    colData = data.frame(
        cell_id = common_cells,
        array_row = coords_subset$x,
        array_col = coords_subset$y,
        row.names = common_cells
    )
)

# Add batch information if available
if (!is.null(batch_info)) {
    colData(sce)$batch <- batch_info[common_cells, "batch"]
}

# Add batch information from original data if available
if ("batch" %in% colnames(merged_obj_original@meta.data)) {
    batch_info_origin <- merged_obj_original@meta.data[common_cells, "batch", drop = FALSE]
    if (!("batch" %in% colnames(colData(sce)))) {
        colData(sce)$batch <- batch_info_origin$batch
    } else {
        cat("Batch information already exists in the SCE object, skipping addition.\n")
    }
}

# Run PCA on cell abundance data
library(scater)
sce <- runPCA(sce, exprs_values = "counts", ncomponents = min(30, nrow(sce) - 1))

# BayesSpace preprocessing
sce <- BayesSpace::spatialPreprocess(sce, platform = "Visium", skip.PCA = TRUE, log.normalize = FALSE)

# inspect the number of clusters for BayesSpace
BayesSpace::qTune(sce, qs = seq(10, 30), platform = "Visium", d = min(15, ncol(reducedDim(sce, "PCA"))), cores = 8) |>
    BayesSpace::qPlot()
# BayesSpace_k_clusters <- 50

# BayesSpace clustering
sce <- BayesSpace::spatialCluster(sce,
    q = BayesSpace_k_clusters, d = min(15, ncol(reducedDim(sce, "PCA"))),
    nrep = 2e+05, burn.in = 10000, platform = "Visium"
)

# Export BayesSpace cluster assignments
bayesspace_results <- data.frame(
    Cell = rownames(colData(sce)),
    BayesSpace_cluster_cell2loc = sce$spatial.cluster,
    stringsAsFactors = FALSE
)

write.csv(bayesspace_results,
    file = paste0(output_dirs$clustering, "BayesSpace_cluster_cell2location.csv"),
    row.names = FALSE
)

# Create Seurat object for visualization only
merged_obj_diet <- CreateSeuratObject(
    counts = t(as.matrix(cell2loc_data_subset)),
    project = "cell2location_BayesSpace",
    assay = "CellAbundance"
)

# Add metadata
merged_obj_diet$spatial_x <- coords_subset$x
merged_obj_diet$spatial_y <- coords_subset$y
merged_obj_diet$BayesSpace_cluster_cell2loc <- sce$spatial.cluster
if (!is.null(batch_info)) {
    merged_obj_diet$batch <- batch_info[common_cells, "batch"]
}
merged_obj_diet@images <- merged_obj_original@images

# Add PCA from SCE for UMAP
pca_embeddings <- reducedDim(sce, "PCA")
colnames(pca_embeddings) <- paste0("PC_", seq_len(ncol(pca_embeddings)))
merged_obj_diet[["pca"]] <- CreateDimReducObject(embeddings = pca_embeddings, key = "PC_", assay = "CellAbundance")

# Create output directory for cell2location BayesSpace results
output_dir_cell2loc <- output_dirs$clustering_spatialAware
dir.create(output_dir_cell2loc, recursive = TRUE, showWarnings = FALSE)

# Generate visualizations
Idents(merged_obj_diet) <- "BayesSpace_cluster_cell2loc"

# Spatial plot
if ("batch" %in% colnames(merged_obj_diet@meta.data)) {
    batch_names <- unique(merged_obj_diet$batch)
    total_cells_per_batch <- table(merged_obj_diet$batch)
} else {
    batch_names <- "All"
    total_cells_per_batch <- length(Cells(merged_obj_diet))
    names(total_cells_per_batch) <- "All"
}

# Plot spatial clusters
cluster_labelling_env$plot_spatial_clusters_with_cell_percentage(merged_obj_diet,
    group_by = "BayesSpace_cluster_cell2loc",
    file_name = "Cell2Location_BayesSpace_clusters",
    plot_title = "BayesSpace clusters (Cell2Location)",
    output_dir = output_dir_cell2loc,
    color_map = color_map,
    total_cells_per_batch = total_cells_per_batch,
    batch_names = batch_names
)

# Run UMAP for visualization
merged_obj_diet <- RunUMAP(merged_obj_diet,
    reduction = "pca", dims = seq_len(min(15, ncol(pca_embeddings))),
    reduction.name = "cell2loc_umap", n.components = 3
)

# Plot UMAP
cluster_labelling_env$plot_umap(
    merged_obj_diet, "cell2loc_umap", "BayesSpace_cluster_cell2loc",
    "Cell2Location_BayesSpace_clusters", "BayesSpace clusters (Cell2Location)",
    output_dir_cell2loc, color_map
)

cat("=====================================\n")
cat("Done running BayesSpace clustering on cell2location results.\n")
cat("Results saved to:", output_dir_cell2loc, "\n")
cat("Number of clusters found:", length(unique(merged_obj_diet$BayesSpace_cluster_cell2loc)), "\n")
cat("=====================================\n")
quit("no")
