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

if (cluster_method %in% colnames(merged_obj@meta.data)) {
    valid_cluster_method <- cluster_method
} else {
    # Load cluster assignments from CSV file
    cluster_csv_path <- paste0(output_dirs$clustering, cluster_method, ".csv")
    if (file.exists(cluster_csv_path)) {
        cluster_data <- read.csv(cluster_csv_path, row.names = 1)
        cluster_vec <- setNames(cluster_data[[1]], rownames(cluster_data))
        merged_obj@meta.data[[cluster_method]] <- cluster_vec[rownames(merged_obj@meta.data)]
        valid_cluster_method <- cluster_method
    } else {
        stop(paste("Cluster method", cluster_method, "not found in metadata or CSV."))
    }
}

Idents(merged_obj) <- valid_cluster_method



# =============================================
# Export raw count for DPNBVI to cluster
# ===========================================
if (!VisiumHD) {
    merged_obj <- JoinLayers(object = merged_obj, layers = "counts")
}
cluster_labelling_env$export_top_variable_raw_count(merged_obj, output_dirs$manual_counts, csv_filenames$merged_batch$raw_count, VisiumHD, csv_filenames$topVar_num)


# ===========================================
#  non-immune/immune cells individual Clustering
# ===========================================
# we first split the cells into two groups based on Ptprc expression
# (CD45+ is the immune cells and CD45- is the non-immune cells)
# and then cluster each group separately
# ===========================================
if (!VisiumHD) {
    merged_obj <- JoinLayers(object = merged_obj, layers = "data")
}
ptprc_expr <- GetAssayData(
    merged_obj,
    assay = ifelse(VisiumHD, "sketch", "Spatial"),
    layer = "data"
)[grep("^Ptprc$", rownames(merged_obj), ignore.case = TRUE), ]
Ptprc_pos <- subset(merged_obj, cells = names(ptprc_expr)[ptprc_expr > 0])
Ptprc_neg <- subset(merged_obj, cells = names(ptprc_expr)[ptprc_expr == 0])

if (!VisiumHD) {
    Ptprc_pos <- JoinLayers(Ptprc_pos, assay = "Spatial")
    Ptprc_neg <- JoinLayers(Ptprc_neg, assay = "Spatial")
}
cluster_labelling_env$export_top_variable_raw_count(Ptprc_pos, output_dirs$manual_counts, csv_filenames$CD45_pos$raw_count, VisiumHD, csv_filenames$topVar_num)
cluster_labelling_env$export_top_variable_raw_count(Ptprc_neg, output_dirs$manual_counts, csv_filenames$CD45_neg$raw_count, VisiumHD, csv_filenames$topVar_num)
