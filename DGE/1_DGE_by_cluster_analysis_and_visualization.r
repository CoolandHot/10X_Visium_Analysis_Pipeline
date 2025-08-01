# ===============================================================================
# DIFFERENTIAL GENE EXPRESSION ANALYSIS AND VISUALIZATION BY CLUSTERS
# ===============================================================================
# This script combines DGE analysis and visualization functionality
# Sections:
# 1. Setup and Data Loading
# 2. Region Assignment
# 3. Differential Gene Expression Analysis
# 4. Data Export and Merging
# 5. Across Sample Spatial Visualization
# ===============================================================================

source("util_headers.r")

# ===============================================================================
# SECTION 1: SETUP AND DATA LOADING
# ===============================================================================

merged_obj <- readRDS(paste0(project_dir, "rds_data/", output.file.prefix, "_clustered_12k.rds"))

# Validate cluster method and load if necessary
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

# ===============================================================================
# SECTION 2: REGION ASSIGNMENT
# ===============================================================================

# Assign region labels based on cluster numbers
cell_region_cluster <- Idents(merged_obj)
merged_obj$cell_region_cluster <- dplyr::case_when(
    cell_region_cluster %in% outTumour_cluster_nums_vector ~ "outTumour",
    cell_region_cluster %in% inTumour_cluster_nums_vector ~ "inTumour",
    cell_region_cluster %in% edgeTumour_cluster_nums_vector ~ "edgeTumour",
    .default = as.character(cell_region_cluster)
) |> as.factor()

# Create batch_cluster identifier for comparisons
merged_obj$batch_cluster <- paste0(merged_obj$batch, "_", merged_obj$cell_region_cluster)

# Split object by batch for visualization
merged_obj_list <- SplitObject(merged_obj, split.by = "batch")

# ===============================================================================
# SECTION 3: DIFFERENTIAL GENE EXPRESSION ANALYSIS
# ===============================================================================

# Helper function for differential expression analysis
diff_expr_DE_csv <- function(group1, group2, cluster_id_str, output_file_path) {
    if (sum(merged_obj[[cluster_id_str]] == group1) > 3 && sum(merged_obj[[cluster_id_str]] == group2) > 3) {
        # In Seurat::FindMarkers: ident.1 is the test group, ident.2 is the reference group
        # Positive avg_log2FC means higher expression in ident.1 (group1) vs ident.2 (group2)
        Seurat::FindMarkers(merged_obj,
            assay = ifelse(VisiumHD, "Spatial.008um", "Spatial"),
            ident.1 = group1, # test group - what we're testing
            ident.2 = group2, # reference group - what we're comparing against
            group.by = cluster_id_str,
            verbose = FALSE,
            recorrect_umi = FALSE # use the existing corrected counts
        ) |>
            # Positive avg_log2FC indicates higher expression in group1 vs group2
            dplyr::arrange(dplyr::desc(avg_log2FC)) |>
            dplyr::relocate(p_val, .after = last_col()) |>
            write.csv(file = output_file_path)
    }
}

# -------------------------
# Across samples comparison
# -------------------------
across_sample_dir <- paste0(output_dirs$deg_clusters, "dge_results_across_sample/")
if (!dir.exists(across_sample_dir)) {
    dir.create(across_sample_dir, recursive = TRUE)
}

for (compare_pair in across_batch_comparisons) {
    batch_id1 <- compare_pair[1]
    batch_id2 <- compare_pair[2]

    # For each region, compare as a whole group
    for (region in c("inTumour", "outTumour", "edgeTumour")) {
        group1 <- paste0(batch_id1, "_", region)
        group2 <- paste0(batch_id2, "_", region)

        # Export scaled data to CSV if not already exported
        group1_cells <- rownames(merged_obj@meta.data)[merged_obj$batch_cluster == group1]
        group2_cells <- rownames(merged_obj@meta.data)[merged_obj$batch_cluster == group2]
        group1_csv <- paste0(output_dirs$deg_clusters, "gene_counts_", region, "_", batch_id1, "_scaled.csv")
        group2_csv <- paste0(output_dirs$deg_clusters, "gene_counts_", region, "_", batch_id2, "_scaled.csv")

        if (!file.exists(group1_csv)) {
            group1_scaled <- as.data.frame(t(as.matrix(Seurat::GetAssayData(merged_obj, slot = "scale.data")[, group1_cells, drop = FALSE])))
            group1_scaled$barcode <- rownames(group1_scaled)
            write.csv(group1_scaled, file = group1_csv, row.names = TRUE)
        }
        if (!file.exists(group2_csv)) {
            group2_scaled <- as.data.frame(t(as.matrix(Seurat::GetAssayData(merged_obj, slot = "scale.data")[, group2_cells, drop = FALSE])))
            group2_scaled$barcode <- rownames(group2_scaled)
            write.csv(group2_scaled, file = group2_csv, row.names = TRUE)
        }

        # Run differential expression analysis
        output_file_path <- paste0(across_sample_dir, region, "_", batch_id1, "_vs_", batch_id2, ".csv")
        diff_expr_DE_csv(
            group1 = group1, group2 = group2,
            cluster_id_str = "batch_cluster",
            output_file_path = output_file_path
        )
    }
}

# -------------------------
# Within samples comparison
# -------------------------
within_sample_dir <- paste0(output_dirs$deg_clusters, "dge_results_within_sample/")
if (!dir.exists(within_sample_dir)) {
    dir.create(within_sample_dir, recursive = TRUE)
}

for (batch_group in batch_names) {
    # For each region pair, compare as a whole group
    for (region1 in c("inTumour", "edgeTumour")) {
        for (region2 in c("outTumour")) {
            group1 <- paste0(batch_group, "_", region1)
            group2 <- paste0(batch_group, "_", region2)

            # Export group1 scaled data to CSV if not already exported
            group1_cells <- rownames(merged_obj@meta.data)[merged_obj$batch_cluster == group1]
            group1_csv <- paste0(output_dirs$deg_clusters, "gene_counts_", region1, "_", batch_group, "_scaled.csv")
            if (!file.exists(group1_csv)) {
                group1_scaled <- as.data.frame(t(as.matrix(Seurat::GetAssayData(merged_obj, slot = "scale.data")[, group1_cells, drop = FALSE])))
                group1_scaled$barcode <- rownames(group1_scaled)
                write.csv(group1_scaled, file = group1_csv, row.names = TRUE)
            }

            diff_expr_DE_csv(
                group1 = group1, group2 = group2,
                cluster_id_str = "batch_cluster",
                output_file_path = paste0(
                    within_sample_dir,
                    region1, "_against_", region2, "_", batch_group, ".csv"
                )
            )
        }
    }
}

# ===============================================================================
# SECTION 4: DATA EXPORT AND MERGING
# ===============================================================================

# Recursively list all CSV files under output_dirs$deg_clusters
all_cluster_de_files <- list.files(
    path = output_dirs$deg_clusters,
    pattern = "\\.csv$",
    recursive = TRUE,
    full.names = TRUE
)
# Only include files in dge_results_across_sample or dge_results_within_sample
all_cluster_de_files <- all_cluster_de_files[
    grepl("dge_results_across_sample|dge_results_within_sample", all_cluster_de_files)
]

merged_data_cluster <- data.frame()
for (de_output in all_cluster_de_files) {
    if (!grepl("\\.csv$", de_output)) {
        cat("Skipping non-CSV file:", de_output, "\n")
        next
    }
    data <- read.csv(de_output)
    data$File_Name <- basename(de_output)
    if (nrow(merged_data_cluster) == 0) {
        merged_data_cluster <- data
    } else {
        merged_data_cluster <- rbind(merged_data_cluster, data)
    }
}

# Adjust parsing of group, region, and comparison columns for new file naming
merged_data_cluster <- dplyr::mutate(merged_data_cluster,
    group = dplyr::case_when(
        grepl("_against_", File_Name) ~ "within_sample",
        TRUE ~ "across_sample"
    ),
    region = dplyr::case_when(
        group == "within_sample" ~ stringr::str_extract(File_Name, "^(inTumour|edgeTumour)_against_outTumour"),
        group == "across_sample" ~ stringr::str_extract(File_Name, "^(inTumour|outTumour|edgeTumour)"),
        TRUE ~ NA_character_
    ),
    comparison = dplyr::case_when(
        group == "within_sample" ~ stringr::str_remove(File_Name, "^(inTumour|edgeTumour)_against_(outTumour)_"),
        group == "across_sample" ~ stringr::str_remove(File_Name, "^(inTumour|outTumour|edgeTumour)_"),
        TRUE ~ File_Name
    ) |> stringr::str_remove("\\.csv$")
)

merged_data_cluster |>
    dplyr::rename(gene = X) |>
    dplyr::select(-File_Name) |>
    dplyr::arrange(gene, dplyr::desc(avg_log2FC)) |>
    dplyr::relocate(gene, avg_log2FC, group, region, comparison, .before = everything()) |>
    write.csv(
        file = paste0(output_dirs$deg_clusters, "merged_dge_on_clusters.csv"),
        row.names = FALSE
    )

# ===============================================================================
# SECTION 5: ACROSS SAMPLE SPATIAL VISUALIZATION
# ===============================================================================

# Helper function to plot and save SpatialDimPlot for a given region
plot_region <- function(subset_obj, region_name, output_path) {
    region_cells <- rownames(subset_obj@meta.data)[subset_obj$cell_region_cluster == region_name]
    if (length(region_cells) == 0) {
        message(paste("No cells found for region", region_name, "in this sample."))
        return()
    }
    p <- SpatialDimPlot(subset_obj,
        cells.highlight = region_cells,
        cols.highlight = c("#FFFF00", "grey50"), facet.highlight = TRUE, combine = TRUE
    ) + NoLegend()
    ggsave(paste0(output_path, region_name, ".pdf"), plot = p)
}

# Across sample spatial visualization
for (compare_pair in across_batch_comparisons) {
    batch_id1 <- compare_pair[1]
    batch_id2 <- compare_pair[2]
    output.compare <- paste0(output_dirs$deg_clusters)
    if (!dir.exists(output.compare)) {
        dir.create(output.compare, recursive = TRUE)
    }
    batch_id1_subset <- merged_obj_list[[batch_id1]]
    DefaultAssay(batch_id1_subset) <- ifelse(VisiumHD, "Spatial.008um", "Spatial")
    Idents(batch_id1_subset) <- valid_cluster_method
    batch_id2_subset <- merged_obj_list[[batch_id2]]
    DefaultAssay(batch_id2_subset) <- ifelse(VisiumHD, "Spatial.008um", "Spatial")
    Idents(batch_id2_subset) <- valid_cluster_method
    for (region in c("outTumour", "inTumour", "edgeTumour")) {
        # Plot for batch_id1
        plot_region(batch_id1_subset, region, paste0(output.compare, batch_id1, "_"))
        # Plot for batch_id2
        plot_region(batch_id2_subset, region, paste0(output.compare, batch_id2, "_"))
    }
}

quit("no")
