library(Seurat)
library(ggplot2)


config <- yaml::read_yaml(paste0("config/batch_config.yaml"))
color_map <- config$color_map |> unlist()
batch_file_names <- config$batch_file_names
batch_file_h5_names <- config$batch_file_h5_names # for non VisiumHD
batch_names <- config$batch_names
output.file.prefix <- config$output_file_prefix
raw_data_dir <- config$raw_data_dir
project_dir <- config$project_dir
rds_data_dir <- config$rds_data_dir

inspect_replicates <- config$inspect_replicates

cluster_method <- config$cluster_method

on_disk <- config$on_disk

cell2location_h5ad <- config$cell2location_h5ad
spaceXR_ref_rds <- config$spaceXR_ref_rds

VisiumHD <- config$VisiumHD
Cluster_resolution <- config$Cluster_resolution
Cluster_n_neighbors <- config$Cluster_n_neighbors


outTumour_cluster_nums_vector <- config$Differential_Gene_Analysis$outTumour_cluster_nums
inTumour_cluster_nums_vector <- config$Differential_Gene_Analysis$inTumour_cluster_nums
edgeTumour_cluster_nums_vector <- config$Differential_Gene_Analysis$edgeTumour_cluster_nums
across_batch_comparisons <- config$Differential_Gene_Analysis$across_batch_comparisons


# only for spatial domain clustering
BayesSpace_k_clusters <- config$BayesSpace_k_clusters


run_dim_reduction <- function(seurat_obj, reduction_name, assay_use = NULL) {
    if (is.null(assay_use)) {
        assay_use <- ifelse(VisiumHD, "sketch", "Spatial")
    }

    # Check if we need to run dimension reduction
    if (!reduction_name %in% names(seurat_obj@reductions)) {
        # Run PCA first if it doesn't exist
        if (!"pca" %in% names(seurat_obj@reductions)) {
            seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, assay = assay_use, selection.method = "vst", nfeatures = 2000)
            seurat_obj <- Seurat::ScaleData(seurat_obj, assay = assay_use, features = Seurat::VariableFeatures(seurat_obj))
            seurat_obj <- Seurat::RunPCA(seurat_obj, assay = assay_use, features = Seurat::VariableFeatures(seurat_obj))
        }

        # Run UMAP
        seurat_obj <- Seurat::RunUMAP(seurat_obj,
            dims = 1:30,
            reduction = "pca",
            assay = assay_use,
            reduction.name = reduction_name,
            n.components = 3, # 3D UMAP
            verbose = TRUE
        )
    }

    return(seurat_obj)
}

export_spatial_cluster_result_to_tiff <- function(seurat_obj, group_by, image_alpha, file_name, output_dir) {
    levels(seurat_obj@meta.data[[group_by]]) <- c(levels(seurat_obj@meta.data[[group_by]]), "Unknown")
    seurat_obj@meta.data[[group_by]][is.na(seurat_obj@meta.data[[group_by]])] <- "Unknown"

    spatial_plots <- Seurat::SpatialDimPlot(seurat_obj,
        group.by = group_by,
        cols = color_map, alpha = 1, image.alpha = image_alpha,
        label = FALSE
    )


    lapply(seq_along(spatial_plots), function(i) {
        p <- spatial_plots[[i]] +
            theme_void() + # Remove all theme elements
            theme(
                legend.position = "none", # Remove legend
                plot.title = element_blank(), # Remove title
                plot.subtitle = element_blank(), # Remove subtitle
                axis.text = element_blank(), # Remove axis text
                axis.title = element_blank(), # Remove axis titles
                axis.ticks = element_blank(), # Remove axis ticks
                panel.border = element_blank(), # Remove panel border
                plot.margin = margin(0, 0, 0, 0) # Remove plot margins
            )

        ggsave(paste0(output_dir, file_name, "_spatial_", batch_names[[i]], ".tiff"),
            plot = p,
            device = "tiff",
            width = 10, # Adjust to match your H&E/IMC image
            height = 10, # Adjust to match your H&E/IMC image
            units = "in",
            dpi = 300, # High resolution for overlay
            compression = "lzw"
        )
    })
}

# Define comprehensive centralized output directory structure from config
output_dirs <- config$output_dirs

# Create all directories
lapply(output_dirs, function(dir) dir.create(dir, recursive = TRUE, showWarnings = FALSE))
