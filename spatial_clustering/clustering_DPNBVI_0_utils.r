cluster_labelling_env <- new.env()

# Export the top variable raw count matrix to a CSV file.
# Args:
#   seurat_obj: A Seurat object containing the data.
#   output_dir: Directory where the output file will be saved.
#   file_name: Name of the output CSV file.
#   visium_hd: Boolean indicating whether VisiumHD data is used.
cluster_labelling_env$export_top_variable_raw_count <- function(seurat_obj, output_dir, file_name, visium_hd = FALSE, nfeatures = 1000) {
    raw_gene_cnt <- GetAssayData(seurat_obj, assay = ifelse(visium_hd, "sketch", "Spatial"), layer = "counts")
    seurat_obj <- FindVariableFeatures(seurat_obj, assay = ifelse(visium_hd, "sketch", "Spatial"), selection.method = "vst", nfeatures = nfeatures)
    variable_features <- VariableFeatures(seurat_obj)
    feature_variances <- apply(raw_gene_cnt[variable_features, ], 1, var) |> sort(decreasing = TRUE)
    ranked_features <- names(feature_variances)
    top_var_gene_cnt <- raw_gene_cnt[ranked_features, ]

    # Remove rows and columns with all zeros
    top_var_gene_cnt <- top_var_gene_cnt[rowSums(top_var_gene_cnt) > 0, colSums(top_var_gene_cnt) > 1]

    top_var_gene_cnt |>
        t() |>
        write.csv(file = paste0(output_dir, file_name), row.names = TRUE)
}

# Load inferred cluster labels and map them to the Seurat object.
# Args:
#   output_dir: Directory containing the label and raw count files.
#   label_file: Name of the file containing inferred labels.
#   raw_count_file: Name of the file containing raw counts.
#   seurat_obj: A Seurat object to which the labels will be applied.
# Returns:
#   A factor of cluster labels mapped to the Seurat object.
cluster_labelling_env$load_inferred_labels <- function(label_file, raw_count_file, seurat_obj) {
    labels <- read.csv(label_file, header = TRUE)$label
    clusters_ext <- factor(labels, labels = seq_along(1:length(unique(labels))))
    names(clusters_ext) <- read.csv(raw_count_file, row.names = 1) |> rownames()
    clusters_ext <- clusters_ext[colnames(seurat_obj)]
    names(clusters_ext) <- colnames(seurat_obj)
    return(clusters_ext)
}
