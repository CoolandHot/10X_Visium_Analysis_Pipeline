source("util_headers.r")

source(paste0(project_dir, "cluster_labelling/cluster_labelling_utils.r"))

csv_filenames <- yaml::read_yaml(paste0(project_dir, "cluster_labelling", "/", "manual_label_settings.yaml"))


merged_obj <- readRDS(paste0(rds_data_dir, output.file.prefix, "_merged.rds"))
if (VisiumHD) {
    DefaultAssay(merged_obj) <- "sketch"
}

if (!VisiumHD) {
    merged_obj <- JoinLayers(object = merged_obj, layers = "counts")
}

clusters_ext <- cluster_labelling_env$load_inferred_labels(
    paste0(output_dirs$manual_assignments, csv_filenames$merged_batch$label),
    paste0(output_dirs$manual_counts, csv_filenames$merged_batch$raw_count),
    merged_obj
)
merged_obj$clusters_ext <- clusters_ext

Ptprc_neg <- readRDS(paste0(project_dir, "rds_data/", output.file.prefix, csv_filenames$CD45_neg$rds_suffix))
Ptprc_pos <- readRDS(paste0(project_dir, "rds_data/", output.file.prefix, csv_filenames$CD45_pos$rds_suffix))



# ===========================================
#  Cluster Visualization: Tumour/Non-Tumour Expression on the Spatial Plot
# ===========================================
# Egfr is the tumour marker
# https://www.nature.com/articles/s41590-022-01215-0
# ===========================================
merged_obj_list <- SplitObject(merged_obj, split.by = "batch")

# Iterate over batches and generate plots
for (batch_id in batch_names) {
    one_subset <- merged_obj_list[[batch_id]]
    Seurat::Idents(one_subset) <- "clusters_ext"

    # Generate PDF and Plotly plots for the tumor marker (e.g., Egfr)
    cluster_labelling_env$pdf_spatial_plot(
        seurat_obj = one_subset,
        batch_id = batch_id,
        marker = "Egfr",
        group_by = "clusters_ext",
        output_dir = output_dirs$manual_markers,
        color_map = color_map
    )
    cluster_labelling_env$plotly_spatial_plot(
        seurat_obj = one_subset,
        batch_id = batch_id,
        marker = "Egfr",
        group_by = "clusters_ext",
        output_dir = output_dirs$manual_markers,
        color_map = color_map
    )
}





# ===========================================
#  Cluster Visualization: non-immune/immune cells on spatial plot
# ===========================================
# based on Ptprc expression
# (CD45+ is the immune cells and CD45- is the non-immune cells)
# and then cluster each group separately
# ===========================================
total_cells_per_batch <- table(merged_obj$batch)

cluster_labelling_env$plot_umap(Ptprc_pos, "cd45_pos_umap", "cd45_pos_clusters_ext", "cd45_pos_clusters", "CD45+ cells", output_dirs$manual_immune, color_map)
cluster_labelling_env$plot_spatial_clusters_with_cell_percentage(Ptprc_pos,
    group_by = "cd45_pos_clusters_ext", file_name = "cd45_pos_clusters", plot_title = "CD45+ cells",
    output_dir = output_dirs$manual_immune, color_map = color_map, total_cells_per_batch = total_cells_per_batch, batch_names = batch_names
)
cluster_labelling_env$plot_umap(Ptprc_neg, "cd45_neg_umap", "cd45_neg_clusters_ext", "cd45_neg_clusters", "CD45- cells", output_dirs$manual_immune, color_map)
cluster_labelling_env$plot_spatial_clusters_with_cell_percentage(Ptprc_neg,
    group_by = "cd45_neg_clusters_ext", file_name = "cd45_neg_clusters", plot_title = "CD45- cells",
    output_dir = output_dirs$manual_immune, color_map = color_map, total_cells_per_batch = total_cells_per_batch, batch_names = batch_names
)

cluster_labelling_env$plot_umap(merged_obj, ifelse(VisiumHD, "sketch_umap", "umap"), "clusters_ext", "AllCells_clusters", "All cells", output_dirs$manual_immune, color_map)
cluster_labelling_env$plot_spatial_clusters_with_cell_percentage(merged_obj,
    group_by = "clusters_ext", file_name = "AllCells_clusters", plot_title = "All cells",
    output_dir = output_dirs$manual_immune, color_map = color_map, total_cells_per_batch = total_cells_per_batch, batch_names = batch_names
)










# ===========================================
Ptprc_neg <- readRDS(paste0(project_dir, "rds_data/", output.file.prefix, csv_filenames$CD45_neg$rds_suffix))
Ptprc_pos <- readRDS(paste0(project_dir, "rds_data/", output.file.prefix, csv_filenames$CD45_pos$rds_suffix))
Seurat::Idents(Ptprc_neg) <- "cd45_neg_clusters_ext"
Seurat::Idents(Ptprc_pos) <- "cd45_pos_clusters_ext"
# ===========================================


# ==============================
# find cluster markers
# ==============================
# find markers for each cluster
Ptprc_neg <- JoinLayers(Ptprc_neg, assay = ifelse(VisiumHD, "sketch", "Spatial"))
Ptprc_pos <- JoinLayers(Ptprc_pos, assay = ifelse(VisiumHD, "sketch", "Spatial"))

Ptprc_neg_markers <- Seurat::FindAllMarkers(
    object = Ptprc_neg,
    assay = ifelse(VisiumHD, "sketch", "Spatial"),
    min.pct = 0.1,
    logfc.threshold = 0.01,
    test.use = "LR",
    verbose = TRUE
)
Ptprc_pos_markers <- Seurat::FindAllMarkers(
    object = Ptprc_pos,
    assay = ifelse(VisiumHD, "sketch", "Spatial"),
    min.pct = 0.1,
    logfc.threshold = 0.01,
    test.use = "LR",
    verbose = TRUE
)
write.csv(Ptprc_neg_markers, file = paste0(output_dirs$manual_immune, "cd45_neg_cluster_markers.csv"))
write.csv(Ptprc_pos_markers, file = paste0(output_dirs$manual_immune, "cd45_pos_cluster_markers.csv"))

Idents(merged_obj) <- "clusters_ext"
merged_obj <- JoinLayers(merged_obj, assay = ifelse(VisiumHD, "sketch", "Spatial"))
cluster_markers <- Seurat::FindAllMarkers(
    object = merged_obj,
    assay = ifelse(VisiumHD, "Spatial.008um", "Spatial"),
    min.pct = 0.1,
    logfc.threshold = 0.01,
    test.use = "LR",
    verbose = TRUE
)
write.csv(cluster_markers, file = paste0(output_dirs$manual_immune, "AllCells_cluster_markers.csv"))


# # ===========================================
# #           Tree Clustering
# # ===========================================
# if (requireNamespace("ape", quietly = TRUE)) {
#     # the buildClusterTree function may not work in some environments
#     Idents(merged_obj) <- "clusters_ext"
#     if (any(is.na(Idents(merged_obj)))) {
#         warning("Some cells have NA identities, removing them before building tree")
#         merged_obj <- subset(merged_obj, cells = names(Idents(merged_obj))[!is.na(Idents(merged_obj))])
#     }
#     merged_obj <- Seurat::BuildClusterTree(merged_obj, reorder = TRUE, reorder.numeric = TRUE)
#     p1 <- Seurat::PlotClusterTree(object = merged_obj)
#     pdf(paste0(output_dirs$manual_immune, "cluster_tree.pdf"))
#     print(p1)
#     dev.off()
# }

# =========================================================================
#          Visualize the expression of every biomarker on UMAP
# =========================================================================
# Define cell types based on known markers
# https://www.cellsignal.com/pathways/immune-cell-markers-mouse

# Load gene patterns from YAML file
gene_patterns <- csv_filenames$gene_patterns

cell_type_genes <- lapply(gene_patterns, function(pattern) {
    grep(pattern, rownames(merged_obj), ignore.case = TRUE, value = TRUE)
})

biomarkers <- unlist(cell_type_genes)

# plot each marker on UMAP, facet by marker
chunks <- split(biomarkers, ceiling(seq_along(biomarkers) / 4))
pdf(paste0(output_dirs$manual_immune, "cd45_neg_expression_on_umap_per_marker.pdf"))
for (chunk in chunks) {
    umap_per_marker <- Seurat::FeaturePlot(Ptprc_neg, features = chunk)
    print(umap_per_marker)
}
dev.off()
pdf(paste0(output_dirs$manual_immune, "cd45_pos_expression_on_umap_per_marker.pdf"))
for (chunk in chunks) {
    umap_per_marker <- Seurat::FeaturePlot(Ptprc_pos, features = chunk)
    print(umap_per_marker)
}
dev.off()



# =============================================
# =========== Clusters -> cell types =========
# =============================================
cluster_markers <- read.csv(paste0(output_dirs$manual_immune, "AllCells_cluster_markers.csv"))
cluster_markers_matches_biomarkers <- cluster_markers |>
    dplyr::filter(gene %in% biomarkers) |>
    dplyr::mutate(possible_cell_types = sapply(gene, function(g) {
        paste(gsub("\\d+$", "", names(biomarkers)[biomarkers == g]), collapse = ", ")
    }))

write.csv(cluster_markers_matches_biomarkers, file = paste0(output_dirs$manual_immune, "cluster_markers_matches_biomarkers.csv"), row.names = FALSE)


# ************************************************
# Identify clusters potentially associated with each cell type
# ************************************************
# categorize cluster numbers to known cell types
cluster_markers_matches_biomarkers <- read.csv(paste0(output_dirs$manual_immune, "cluster_markers_matches_biomarkers.csv"))

cell_type_clusters <- lapply(cell_type_genes, function(gene_list) {
    cluster_markers_matches_biomarkers$cluster[cluster_markers_matches_biomarkers$gene %in% gene_list] |>
        as.vector() |>
        unique()
})

# Create a function to find unique numbers
find_unique_numbers <- function(cluster_list) {
    # Create a frequency table of all numbers
    all_numbers <- unlist(cluster_list)
    number_counts <- table(all_numbers)

    # Find which clusters each number belongs to
    unique_results <- list()

    for (num in names(number_counts)) {
        present_in <- names(cluster_list)[sapply(cluster_list, function(x) num %in% x)]
        if (length(present_in) == 1) {
            unique_results[[num]] <- present_in
        }
    }

    return(unique_results)
}

# Find the unique cluster that belonging to a cell type
unique_numbers <- find_unique_numbers(cell_type_clusters)

print(unique_numbers)



# ====================== Assign cell types to clusters ======================
# Assign cell types to clusters
# ============================================================================
# =====.......... because the above markers can't be individually distinguished,
# =====.......... below jobs do not work. Stop here.
# ============================================================================

# seurat_clusters <- Idents(merged_obj)
# collapsed <- rep(NA, length(seurat_clusters))

# for (i in seq_along(cell_type_clusters)) {
#     ind <- which(seurat_clusters %in% cell_type_clusters[[i]])
#     collapsed[ind] <- names(cell_type_clusters)[i]
# }
# collapsed <- factor(collapsed, levels = names(cell_type_clusters), ordered = TRUE)
# names(collapsed) <- names(seurat_clusters)

# merged_obj$manual_cell_type <- collapsed
# Idents(merged_obj) <- "manual_cell_type"

# p1 <- DimPlot(merged_obj, reduction = "umap", label = TRUE) + NoAxes() + NoLegend() + ggtitle("Cell Types made by manual labelling from biomarkers")
# ggsave(paste0(output_dirs$manual_immune, "cell_types_umap.pdf"), p1, width = 10, height = 10)




# ************************************************
# The above method returns too few identified cell types on clusters
# ************************************************
# =============================================
#  Manual study the cluster markers
#   with the associated cell types
#  determined by voting
# =============================================
cluster_markers_matches_biomarkers <- read.csv(paste0(output_dirs$manual_immune, "cluster_markers_matches_biomarkers.csv"))


# # check an individual cluster
# cluster_id <- 10
# cluster_markers_matches_biomarkers |>
#     dplyr::filter(abs(avg_log2FC) > 1 & p_val_adj < 0.05) |>
#     dplyr::filter(cluster == cluster_id) |>
#     dplyr::select(possible_cell_types) |>
#     dplyr::pull() |>
#     determine_cell_type_by_votes()


all_clusters_cell_types <- cluster_markers_matches_biomarkers |>
    dplyr::filter(abs(avg_log2FC) > 1 & p_val_adj < 0.05) |>
    dplyr::group_by(cluster) |>
    dplyr::summarise(
        most_likely_cell_type = cluster_labelling_env$determine_cell_type_by_votes(possible_cell_types)
    )

print(all_clusters_cell_types, n = Inf)



# ====================== Assign cell types to clusters ======================
# Assign cell types to clusters
# ============================================================================

seurat_clusters <- Idents(merged_obj)

collapsed <- sapply(seurat_clusters, function(cluster) {
    match <- all_clusters_cell_types$most_likely_cell_type[all_clusters_cell_types$cluster == cluster]
    if (length(match) > 0) {
        return(match)
    } else {
        return("Unassigned")
    }
})

unique_cell_types <- c(unique(all_clusters_cell_types$most_likely_cell_type), "Unassigned")

collapsed <- factor(collapsed, levels = unique_cell_types, ordered = TRUE)
names(collapsed) <- names(seurat_clusters)

merged_obj$manual_cell_type <- collapsed
Idents(merged_obj) <- "manual_cell_type"
data.frame(Cell = rownames(merged_obj@meta.data), manual_cell_type = merged_obj$manual_cell_type) |>
    write.csv(file = paste0(output_dirs$manual_assignments, "manual_cell_type.csv"), row.names = FALSE)

# Ensure color_map includes all levels of manual_cell_type
if (!"Unassigned" %in% unique_cell_types) {
    unique_cell_types <- c(unique_cell_types, "Unassigned")
}
color_map <- setNames(color_map[seq_along(unique_cell_types)], unique_cell_types)

p1 <- DimPlot(merged_obj, cols = color_map, reduction = "umap", label = TRUE) +
    NoAxes() +
    ggtitle("Cell Types made by manual labelling from biomarkers") +
    theme(
        legend.position = "bottom",
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 8),
        legend.box.margin = margin(t = 10)
    ) +
    guides(color = guide_legend(ncol = 3, override.aes = list(size = 3)))

# Create separate spatial plots for each batch with specific legend settings
batch_ids <- names(merged_obj@images)
p2_list <- list()

for (i in seq_along(batch_ids)) {
    p2_list[[i]] <- SpatialDimPlot(merged_obj,
        images = batch_ids[i],
        group.by = "manual_cell_type",
        cols = color_map,
        alpha = 1,
        image.alpha = 0.7,
        label = TRUE,
        repel = TRUE,
        label.size = 4
    ) +
        theme(
            legend.position = "bottom",
            legend.key.width = unit(0.5, "cm"),
            legend.text = element_text(size = 8),
            legend.box.margin = margin(t = 10)
        ) +
        guides(fill = guide_legend(ncol = 3, override.aes = list(size = 3)))
}


ggpubr::ggexport(
    plotlist = c(list(p1), p2_list),
    filename = paste0(output_dirs$manual_immune, "manual_labelling_cell_types.pdf"),
    width = 15,
    height = 15
)
