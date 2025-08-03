source("util_headers.r")

merged_obj <- readRDS(paste0(rds_data_dir, output.file.prefix, "_merged_clustered.rds"))



# ================= SPECIFIC to KD PROJECT =================
# Remove all cells that are not in the "0" cluster
merged_obj <- subset(merged_obj, cells = colnames(merged_obj)[merged_obj$seurat_clusters == "0"])
cat(length(colnames(merged_obj)), " cells remain after subsetting to cluster 0", "\n")
# Plot total counts by batch
batch_counts <- data.frame(
    batch = merged_obj$batch,
    nCount = merged_obj$nCount_Spatial.008um
)
total_counts_by_batch <- aggregate(nCount ~ batch, data = batch_counts, sum)
p_counts <- ggplot(total_counts_by_batch, aes(x = reorder(batch, nCount), y = nCount, fill = batch)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(title = "Total counts by batch (only on cluster 0)", x = "Batch", y = "Total counts") +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)
    ) +
    NoLegend()
ggsave(paste0(output_dirs$clustering, "batch_total_counts_only_on_cluster_0.pdf"), plot = p_counts)

# ================= SPECIFIC to KD PROJECT =================





# ===========================================
#      visualization of all cells clustering
# ===========================================
#  plots UMAP in 2D
Seurat::DefaultAssay(merged_obj) <- ifelse(VisiumHD, "sketch", "Spatial")
Seurat::Idents(merged_obj) <- "seurat_clusters"
p1 <- Seurat::DimPlot(merged_obj, cols = color_map, reduction = "umap", label = F) +
    ggtitle("Sketched clustering") +
    theme(legend.position = "bottom")
ggplot2::ggsave(paste0(output_dirs$clustering_umap_spatial, output.file.prefix, "_clusters_sketch.pdf"), plot = p1)

#  plots UMAP in 3D
library(plotly)
umap_3d <- merged_obj[["umap"]]@cell.embeddings
Idents(merged_obj) <- "seurat_clusters"
idents <- Idents(merged_obj) |> na.omit()
plotly::plot_ly(
    x = umap_3d[, 1], y = umap_3d[, 2], z = umap_3d[, 3],
    type = "scatter3d", mode = "markers",
    color = idents,
    colors = color_map,
    marker = list(size = 2)
) |>
    htmlwidgets::saveWidget(paste0(output_dirs$clustering_umap_spatial, output.file.prefix, "_clusters_sketch_3D.html"))

if (VisiumHD) {
    # ===========================================
    # switch to full dataset
    # ===========================================
    Seurat::DefaultAssay(merged_obj) <- "Spatial.008um"
    Seurat::Idents(merged_obj) <- "seurat_cluster_full"
    p2 <- Seurat::DimPlot(merged_obj, cols = color_map, reduction = "full.umap", label = F) +
        ggtitle("Projected clustering (full dataset)") +
        theme(legend.position = "bottom")
    ggplot2::ggsave(paste0(output_dirs$clustering_umap_spatial, output.file.prefix, "_clusters_projected.pdf"), plot = p2)

    # ===========================================
    #           plots UMAP in 3D
    # ===========================================
    umap_3d <- merged_obj[["full.umap"]]@cell.embeddings
    idents <- merged_obj$seurat_cluster_full |> na.omit()
    plotly::plot_ly(
        x = umap_3d[, 1], y = umap_3d[, 2], z = umap_3d[, 3],
        type = "scatter3d", mode = "markers",
        color = idents,
        colors = color_map,
        marker = list(size = 2)
    ) |>
        htmlwidgets::saveWidget(paste0(output_dirs$clustering_umap_spatial, output.file.prefix, "_clusters_projected_3D.html"))
}



# ======================================
# spatial plots for individual inspection
# ======================================
merged_obj_list <- SplitObject(merged_obj, split.by = "batch")

for (batch_id in batch_names) {
    output.batch <- paste0(output_dirs$clustering_umap_spatial, batch_id, "/")
    if (!dir.exists(output.batch)) {
        dir.create(output.batch, recursive = TRUE)
    }

    one_subset <- merged_obj_list[[batch_id]]
    #######################
    # make cluster tree for all top 5 markers across clusters within sample.
    #######################
    # downsample to make visualization easier
    object_subset <- subset(one_subset, cells = Cells(one_subset[[ifelse(VisiumHD, "Spatial.008um", "Spatial")]]), downsample = 1000)
    # Order clusters by similarity
    # object_subset <- Seurat::BuildClusterTree(object_subset, assay = "Spatial.008um", reduction = "full.sketch_integrated.rpca", reorder = T)

    # markers <- Seurat::FindAllMarkers(object_subset, assay = ifelse(VisiumHD, "Spatial.008um", "Spatial"))
    # write.csv(markers, paste0(output.batch, "all_markers_", batch_id, ".csv"))

    # top5 <- markers |>
    #     dplyr::group_by(cluster) |>
    #     dplyr::arrange(gene, dplyr::desc(avg_log2FC)) |>
    #     dplyr::filter(p_val_adj < 0.05) |>
    #     dplyr::slice_head(n = 5) |>
    #     dplyr::ungroup()

    # clusters_with_markers <- unique(top5$cluster)
    # object_subset <- Seurat::ScaleData(object_subset, assay = ifelse(VisiumHD, "Spatial.008um", "Spatial"), features = top5$gene)
    # object_subset <- subset(object_subset, idents = clusters_with_markers)
    # p5 <- Seurat::DoHeatmap(object_subset, assay = ifelse(VisiumHD, "Spatial.008um", "Spatial"), features = top5$gene, size = 2.5) +
    #     theme(axis.text = element_text(size = 5.5)) +
    #     NoLegend() +
    #     ggtitle(paste0("Top 5 markers for each cluster in ", batch_id))
    # ggsave(paste0(output.batch, "top5_markers_heatmap", batch_id, ".pdf"), plot = p5)

    #######################
    # spatial plot with histology image as the background
    #######################
    if (VisiumHD) {
        Seurat::DefaultAssay(one_subset) <- "sketch"
        Seurat::Idents(one_subset) <- "seurat_clusters"
        one_subset <- Seurat::DietSeurat(one_subset, assays = "sketch")
    }

    p3 <- Seurat::SpatialDimPlot(one_subset,
        group.by = "seurat_clusters",
        cols = color_map, alpha = 0.7, image.alpha = 0.7,
        label = T, repel = T,
        label.size = 4, pt.size.factor = ifelse(VisiumHD, 4, 1.6)
    ) + ggtitle("cluster by SNN")
    ggsave(paste0(output.batch, "clusters_spatial_", batch_id, ".pdf"), plot = p3)


    #######################
    # interactive plots
    #######################
    idents <- Idents(one_subset) |> na.omit()
    st_location <- Seurat::GetTissueCoordinates(one_subset)[names(idents), ]
    st_location$cell <- NULL
    plotly::plot_ly(
        x = st_location[, 1], y = st_location[, 2],
        type = "scatter", mode = "markers",
        color = idents,
        colors = color_map,
        marker = list(size = 2)
    ) |>
        plotly::layout(paper_bgcolor = "black", plot_bgcolor = "black", font = list(color = "white")) |>
        htmlwidgets::saveWidget(paste0(output.batch, "clusters_spatial_sketch_", batch_id, ".html"))
}

cat("=====================================\n")
cat("Done visualization\n")
cat("=====================================\n")
q("no")
