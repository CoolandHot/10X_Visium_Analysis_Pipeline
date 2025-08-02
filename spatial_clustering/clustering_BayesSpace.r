# ================================================
# Identifying spatially-defined tissue domains using BayesSpace
# https://benchmarkst-reproducibility.readthedocs.io/en/latest/BayesSpace.html
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
} else {
    stop("Failed to extract 'x' and 'y' from spatial coordinates. Check image data in merged_obj_diet.")
}

# ==========================================
# spatial domain clustering using BayesSpace
# https://www.ezstatconsulting.com/BayesSpace/articles/BayesSpace.html
# ==========================================
library(SingleCellExperiment)

merged_obj_diet <- DietSeurat(merged_obj,
    assays = ifelse(VisiumHD, "sketch", "Spatial"),
    dimreducs = "PCA",
    graphs = NULL
)
cols_to_keep <- c("orig.ident", "nCount_Spatial", "nFeature_Spatial")
merged_obj_diet@meta.data <- merged_obj_diet@meta.data[, cols_to_keep, drop = FALSE]

sce <- Seurat::as.SingleCellExperiment(merged_obj_diet)
# confirm spatial info are in SCE
if (!all(c("array_row", "array_col") %in% colnames(colData(sce)))) {
    # Ensure coords_aligned is in the same order as sce cells
    coords_for_sce <- coords_aligned[rownames(colData(sce)), , drop = FALSE]
    colData(sce)$array_row <- coords_for_sce$spatial_x
    colData(sce)$array_col <- coords_for_sce$spatial_y
}

# platform If "Visium", select six neighboring spots around center;
# if "ST" (Spatial Transcriptomics), select four adjacent spots.
sce <- BayesSpace::spatialPreprocess(sce, platform = "Visium", skip.PCA = T, log.normalize = FALSE)

# Selecting the number of clusters
sce <- BayesSpace::qTune(sce, qs = seq(5, 30), platform = "Visium", d = 7)
BayesSpace::qPlot(sce)

# Spatial clustering using BayesSpace
sce <- BayesSpace::spatialCluster(sce,
    q = BayesSpace_k_clusters, d = min(15, ncol(reducedDim(sce, "PCA"))), # Number of PCs to use
    nrep = 2e+05, burn.in = 10000, platform = "Visium"
)
# Export BayesSpace cluster assignments to CSV
data.frame(Cell = rownames(merged_obj@meta.data), BayesSpace_cluster = sce$spatial.cluster) |>
    write.csv(file = paste0(output_dirs$clustering, "BayesSpace_cluster.csv"), row.names = FALSE)

# Read the CSV and set 'Cell' as rownames
BayesSpace_cluster_df <- read.csv(paste0(output_dirs$clustering, "BayesSpace_cluster.csv"), stringsAsFactors = FALSE, check.names = FALSE)
rownames(BayesSpace_cluster_df) <- BayesSpace_cluster_df$Cell
BayesSpace_cluster_df$Cell <- NULL

# enhance the resolution of the PCs, and add these PCs as well as predicted cluster labels at subspot resolution
# sce <- BayesSpace::spatialEnhance(sce,
#   q = BayesSpace_k_clusters, d = 15, # d is the PCs to use
#   nrep = 2e+05, burn.in = 10000, platform = "ST" # only ST allows for array_row and array_col in spatialEnhance
# )

# Visualize BayesSpace clusters using BayesSpace's plot
# BayesSpace::clusterPlot(sce)

merged_obj$BayesSpace_cluster <- BayesSpace_cluster_df
Idents(merged_obj) <- "BayesSpace_cluster" # Set idents on the full object for plotting

cluster_labelling_env$plot_spatial_clusters_with_cell_percentage(merged_obj,
    group_by = "BayesSpace_cluster", file_name = "AllCells_BayesSpace_clusters", plot_title = "All cells",
    output_dir = output_dirs$clustering_spatialAware, color_map = color_map, total_cells_per_batch = table(merged_obj$batch), batch_names = batch_names
)
cluster_labelling_env$plot_umap(merged_obj, ifelse(VisiumHD, "sketch_umap", "umap"), "BayesSpace_cluster", "AllCells_BayesSpace_clusters", "All cells", output_dirs$clustering_spatialAware, color_map)

cat("=====================================\n")
cat("Done inferring spatial domain clustering using BayesSpace.\n")
cat("=====================================\n")
quit("no")
