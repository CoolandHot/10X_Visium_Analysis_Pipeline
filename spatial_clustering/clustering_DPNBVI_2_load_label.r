source("util_headers.r")

source(paste0(project_dir, "spatial_clustering", "/", "clustering_DPNBVI_0_utils.r"))

csv_filenames <- yaml::read_yaml(paste0(project_dir, "spatial_clustering", "/", "cluster_DPNBVI_config.yaml"))


merged_obj <- readRDS(paste0(
    project_dir, "rds_data/",
    output.file.prefix,
    csv_filenames$merged_batch$rds_suffix
))
if (VisiumHD) {
    DefaultAssay(merged_obj) <- "sketch"
}


# ==============================
# save the cluster information back to the original object
# ==============================
clusters_ext <- cluster_labelling_env$load_inferred_labels(
    paste0(output_dirs$manual_assignments, csv_filenames$merged_batch$label),
    paste0(output_dirs$manual_counts, csv_filenames$merged_batch$raw_count),
    merged_obj
)
merged_obj$clusters_ext <- clusters_ext
# Export cluster assignments to CSV
data.frame(Cell = rownames(merged_obj@meta.data), clusters_ext = clusters_ext) |>
    write.csv(file = paste0(output_dirs$clustering, "clusters_ext.csv"), row.names = FALSE)


# ==============================
# Load and Plot the clusters
# ==============================
source(paste0(project_dir, "cluster_labelling/cluster_labelling_utils.r"))

# Read the CSV and set 'Cell' as rownames
clusters_ext_df <- read.csv(paste0(output_dirs$clustering, "clusters_ext.csv"), stringsAsFactors = FALSE, check.names = FALSE)
rownames(clusters_ext_df) <- clusters_ext_df$Cell
clusters_ext_df$Cell <- NULL

merged_obj$clusters_ext <- clusters_ext_df
Idents(merged_obj) <- "clusters_ext" # Set idents on the full object for plotting

cluster_labelling_env$plot_spatial_clusters_with_cell_percentage(merged_obj,
    group_by = "clusters_ext", file_name = "AllCells_clusters_exts", plot_title = "All cells",
    output_dir = output_dirs$clustering_spatialAware, color_map = color_map, total_cells_per_batch = table(merged_obj$batch), batch_names = batch_names
)

cluster_labelling_env$plot_umap(merged_obj, "umap", "clusters_ext", "AllCells_clusters_exts", "All cells", output_dirs$clustering_spatialAware, color_map)

# ==============================







# ==============================
# save the cluster information back to the original object
# ==============================
# Load the Ptprc_pos and Ptprc_neg objects
Ptprc_pos$cd45_pos_clusters_ext <- cluster_labelling_env$load_inferred_labels(
    paste0(output_dirs$manual_assignments, csv_filenames$CD45_pos$label),
    paste0(output_dirs$manual_counts, csv_filenames$CD45_pos$raw_count),
    Ptprc_pos
)
Ptprc_pos <- run_dim_reduction(Ptprc_pos, "cd45_pos_umap")
saveRDS(Ptprc_pos, paste0(project_dir, "rds_data/", output.file.prefix, csv_filenames$CD45_pos$rds_suffix))

Ptprc_neg$cd45_neg_clusters_ext <- cluster_labelling_env$load_inferred_labels(
    paste0(output_dirs$manual_assignments, csv_filenames$CD45_neg$label),
    paste0(output_dirs$manual_counts, csv_filenames$CD45_neg$raw_count),
    Ptprc_neg
)
Ptprc_neg <- run_dim_reduction(Ptprc_neg, "cd45_neg_umap")
saveRDS(Ptprc_neg, paste0(project_dir, "rds_data/", output.file.prefix, csv_filenames$CD45_neg$rds_suffix))
