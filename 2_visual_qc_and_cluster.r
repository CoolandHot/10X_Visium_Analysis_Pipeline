source("util_headers.r")

merged_obj <- readRDS(paste0(rds_data_dir, output.file.prefix, "_merged_clustered.rds"))

target_assay <- ifelse(VisiumHD, "sketch", "Spatial")

# ================= SPECIFIC to KD PROJECT =================
# 1. Remove noise spots from the merged object
# 2. Re-cluster the remaining cells in cluster "0" (the spots that have cells)
# 3. Save the re-clustered object without noise spots
no_noise_subset_path <- paste0(rds_data_dir, "KD_3Runs_merged_no_noiseSpots", "_merged.rds")
if (file.exists(no_noise_subset_path)) {
    cat("Loading merged object without noise spots (only cells previously in Cluster0)\n")
    merged_sub_obj <- readRDS(no_noise_subset_path)
} else {
    # Remove all cells that are not in the "0" cluster
    merged_sub_obj <- subset(merged_obj, cells = colnames(merged_obj)[merged_obj$seurat_clusters == "0"])
    cat(length(colnames(merged_sub_obj)), " cells remain after subsetting to cluster 0", "\n")
    merged_sub_obj$seurat_clusters <- NULL

    # Plot total counts by batch again, but only for cells in cluster 0
    batch_counts <- data.frame(
        batch = merged_sub_obj$batch,
        nCount = merged_sub_obj$nCount_Spatial.008um
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
    rm(batch_counts, total_counts_by_batch, p_counts)

    # re-cluster on the remaining cells
    merged_sub_obj <- FindNeighbors(
        merged_sub_obj,
        assay = target_assay,
        reduction = "PCA",
        dims = 1:30,
        k.param = 30
    ) |>
        FindClusters(cluster.name = "seurat_clusters", resolution = 0.3) |>
        RunUMAP(
            reduction = "PCA",
            reduction.name = "umap",
            return.model = TRUE,
            dims = 1:30,
            n.components = 3
        )

    saveRDS(merged_sub_obj, no_noise_subset_path)
}

# export cluster results to csv
data.frame(
    cell = Cells(merged_sub_obj),
    cluster = merged_sub_obj$seurat_clusters
) |>
    write.csv(paste0(output_dirs$clustering, "seurat_clusters_no_noiseSpots.csv"), row.names = FALSE)

merged_obj <- merged_sub_obj
rm(merged_sub_obj)
# ================= SPECIFIC to KD PROJECT =================

# ===========================================
#      QC INSPECTION PLOTS
# ===========================================
cat("Creating QC inspection plots...\n")

# Calculate QC metrics if not already present
if (!"nCount_Spatial" %in% colnames(merged_obj@meta.data)) {
    merged_obj$nCount_Spatial <- merged_obj$nCount_Spatial.008um
}
if (!"nFeature_Spatial" %in% colnames(merged_obj@meta.data)) {
    merged_obj$nFeature_Spatial <- merged_obj$nFeature_Spatial.008um
}

# Create QC plots for the full merged object
qc_data <- data.frame(
    total_counts = merged_obj$nCount_Spatial,
    n_genes = merged_obj$nFeature_Spatial,
    batch = merged_obj$batch
)

# Create histogram plots
library(gridExtra)
p_hist1 <- ggplot(qc_data, aes(x = total_counts)) +
    geom_histogram(bins = 50, alpha = 0.7, fill = "steelblue") +
    labs(x = "Total counts per cell", y = "Number of cells", title = "Total Counts Distribution") +
    theme_minimal()

p_hist2 <- ggplot(qc_data, aes(x = n_genes)) +
    geom_histogram(bins = 50, alpha = 0.7, fill = "darkgreen") +
    labs(x = "Number of genes per cell", y = "Number of cells", title = "Genes per Cell Distribution") +
    theme_minimal()

p_hist3 <- ggplot(qc_data, aes(x = total_counts, fill = batch)) +
    geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
    labs(x = "Total counts per cell", y = "Number of cells", title = "Total Counts by Batch") +
    theme_minimal() +
    theme(legend.position = "bottom")

# Create scatter plots
p_scatter1 <- ggplot(qc_data, aes(x = total_counts, y = n_genes)) +
    geom_point(alpha = 0.5, size = 0.5) +
    labs(x = "Total counts", y = "Number of genes", title = "Counts vs Genes") +
    theme_minimal()

p_scatter2 <- ggplot(qc_data, aes(x = total_counts, y = n_genes, color = batch)) +
    geom_point(alpha = 0.5, size = 0.5) +
    labs(x = "Total counts", y = "Number of genes", title = "Counts vs Genes by Batch") +
    theme_minimal() +
    theme(legend.position = "bottom")

# Combine histogram plots
qc_histograms <- grid.arrange(p_hist1, p_hist2, p_hist3, ncol = 3)
ggsave(paste0(output_dirs$clustering, "qc_histograms.pdf"), plot = qc_histograms, width = 15, height = 5)

# Combine scatter plots
qc_scatters <- grid.arrange(p_scatter1, p_scatter2, ncol = 2)
ggsave(paste0(output_dirs$clustering, "qc_scatter_plots.pdf"), plot = qc_scatters, width = 12, height = 5)

# Create spatial QC plots for each batch
pdf(paste0(output_dirs$clustering, "qc_spatial_plots.pdf"), width = 16, height = 6)
for (batch_id in batch_names) {
    one_subset <- subset(merged_obj, subset = batch == batch_id)

    # Spatial plot with total counts
    p_spatial_counts <- SpatialFeaturePlot(one_subset,
        features = "nCount_Spatial",
        alpha = 0.7, image.alpha = 0.5,
        pt.size.factor = ifelse(VisiumHD, 4, 1.6)
    ) +
        ggtitle(paste0("Total Counts - ", batch_id))

    # Spatial plot with number of genes
    p_spatial_genes <- SpatialFeaturePlot(one_subset,
        features = "nFeature_Spatial",
        alpha = 0.7, image.alpha = 0.5,
        pt.size.factor = ifelse(VisiumHD, 4, 1.6)
    ) +
        ggtitle(paste0("Number of Genes - ", batch_id))

    grid.arrange(p_spatial_counts, p_spatial_genes, ncol = 2)
}
dev.off()

# Create summary statistics table
qc_summary <- qc_data |>
    dplyr::group_by(batch) |>
    dplyr::summarise(
        n_cells = dplyr::n(),
        mean_counts = round(mean(total_counts), 2),
        median_counts = round(median(total_counts), 2),
        mean_genes = round(mean(n_genes), 2),
        median_genes = round(median(n_genes), 2),
        .groups = "drop"
    )

qc_summary_long <- reshape2::melt(qc_summary, id.vars = "batch")

# Create bar plots for summary statistics
p_summary_counts <- ggplot(qc_summary, aes(x = batch, y = mean_counts, fill = batch)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste0("n=", n_cells)), vjust = -0.5, size = 3) +
    labs(title = "Mean Total Counts by Batch", x = "Batch", y = "Mean Total Counts") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

p_summary_genes <- ggplot(qc_summary, aes(x = batch, y = mean_genes, fill = batch)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste0("n=", n_cells)), vjust = -0.5, size = 3) +
    labs(title = "Mean Number of Genes by Batch", x = "Batch", y = "Mean Number of Genes") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

# Combine summary plots
qc_summary_plot <- grid.arrange(p_summary_counts, p_summary_genes, ncol = 2)
ggsave(paste0(output_dirs$clustering, "qc_summary_plot.pdf"), plot = qc_summary_plot, width = 12, height = 5)

write.csv(qc_summary, paste0(output_dirs$clustering, "qc_summary_statistics.csv"), row.names = FALSE)

cat("QC inspection plots completed and saved.\n")

# ===========================================
#      MITOCHONDRIAL GENE INSPECTION
# ===========================================
cat("Creating mitochondrial gene inspection plots...\n")

# Calculate mitochondrial gene percentage
mito_genes <- grep("^[Mm][Tt]-|^MT-", rownames(merged_obj), value = TRUE)
if (length(mito_genes) > 0) {
    target_assay_mito <- ifelse(VisiumHD, "Spatial.008um", "Spatial")
    merged_obj <- PercentageFeatureSet(merged_obj, pattern = "^[Mm][Tt]-|^MT-", col.name = "percent.mt", assay = target_assay_mito)

    # Create mitochondrial QC data
    mito_qc_data <- data.frame(
        percent_mt = merged_obj$percent.mt,
        total_counts = merged_obj$nCount_Spatial,
        n_genes = merged_obj$nFeature_Spatial,
        batch = merged_obj$batch
    )

    # Mitochondrial percentage distribution
    p_mito_hist <- ggplot(mito_qc_data, aes(x = percent_mt)) +
        geom_histogram(bins = 50, alpha = 0.7, fill = "red") +
        labs(x = "Mitochondrial Gene %", y = "Number of cells", title = "Mitochondrial Gene % Distribution") +
        theme_minimal() +
        geom_vline(xintercept = c(10, 20, 30), linetype = "dashed", alpha = 0.5, color = "blue")

    # Mitochondrial percentage by batch
    p_mito_batch <- ggplot(mito_qc_data, aes(x = batch, y = percent_mt, fill = batch)) +
        geom_violin(alpha = 0.7) +
        geom_boxplot(width = 0.1, fill = "white", alpha = 0.8) +
        labs(title = "Mitochondrial Gene % by Batch", x = "Batch", y = "Mitochondrial Gene %") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
        geom_hline(yintercept = c(10, 20, 30), linestyle = "dashed", alpha = 0.5, color = "red")

    # Scatter plots with mitochondrial percentage
    p_mito_scatter1 <- ggplot(mito_qc_data, aes(x = total_counts, y = percent_mt)) +
        geom_point(alpha = 0.5, size = 0.5) +
        labs(x = "Total counts", y = "Mitochondrial Gene %", title = "Counts vs Mitochondrial %") +
        theme_minimal()

    p_mito_scatter2 <- ggplot(mito_qc_data, aes(x = n_genes, y = percent_mt)) +
        geom_point(alpha = 0.5, size = 0.5) +
        labs(x = "Number of genes", y = "Mitochondrial Gene %", title = "Genes vs Mitochondrial %") +
        theme_minimal()

    # Combine mitochondrial plots
    mito_plots <- grid.arrange(p_mito_hist, p_mito_batch, p_mito_scatter1, p_mito_scatter2, ncol = 2)
    ggsave(paste0(output_dirs$clustering, "mitochondrial_qc_plots.pdf"), plot = mito_plots, width = 12, height = 10)

    # Spatial mitochondrial plots for each batch
    pdf(paste0(output_dirs$clustering, "mitochondrial_spatial_plots.pdf"), width = 16, height = 6)
    for (batch_id in batch_names) {
        one_subset <- subset(merged_obj, subset = batch == batch_id)

        # Spatial plot with mitochondrial percentage
        p_spatial_mito <- SpatialFeaturePlot(one_subset,
            features = "percent.mt",
            alpha = 0.7, image.alpha = 0.5,
            pt.size.factor = ifelse(VisiumHD, 4, 1.6)
        ) +
            ggtitle(paste0("Mitochondrial % - ", batch_id))

        # Spatial plot with total counts for comparison
        p_spatial_counts <- SpatialFeaturePlot(one_subset,
            features = "nCount_Spatial",
            alpha = 0.7, image.alpha = 0.5,
            pt.size.factor = ifelse(VisiumHD, 4, 1.6)
        ) +
            ggtitle(paste0("Total Counts - ", batch_id))

        grid.arrange(p_spatial_mito, p_spatial_counts, ncol = 2)
    }
    dev.off()

    cat("Mitochondrial gene inspection plots completed and saved.\n")
} else {
    cat("No mitochondrial genes found in the dataset.\n")
}



# ==============================================
# ===========================================
#      visualization of all cells clustering
# ===========================================
# ==============================================
#  plots UMAP in 2D
Seurat::DefaultAssay(merged_obj) <- target_assay
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

spatial_dimplots <- list()

for (batch_id in batch_names) {
    # output.batch <- paste0(output_dirs$clustering_umap_spatial, batch_id, "/")
    output.batch <- output_dirs$clustering_umap_spatial
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

    spatial_dimplots[[batch_id]] <- Seurat::SpatialDimPlot(one_subset,
        group.by = "seurat_clusters",
        cols = color_map, alpha = 0.7, image.alpha = 0.7,
        label = T, repel = T,
        label.size = 4, pt.size.factor = ifelse(VisiumHD, 4, 1.6)
    ) + ggtitle(paste0("cluster by SNN in ", batch_id))


    #######################
    # interactive plots
    #######################
    idents <- Idents(one_subset) |> na.omit()
    st_location <- Seurat::GetTissueCoordinates(one_subset)[names(idents), ]
    st_location$cell <- NULL
    plotly::plot_ly(
        x = st_location[, 2], y = -st_location[, 1],
        type = "scatter", mode = "markers",
        color = idents,
        colors = color_map,
        marker = list(size = 2)
    ) |>
        plotly::layout(paper_bgcolor = "black", plot_bgcolor = "black", font = list(color = "white")) |>
        htmlwidgets::saveWidget(paste0(output.batch, "clusters_spatial_sketch_", batch_id, ".html"))
}

ggpubr::ggexport(
    plotlist = spatial_dimplots,
    filename = paste0(output_dirs$clustering_umap_spatial, output.file.prefix, "_clusters_spatial.pdf")
)

cat("=====================================\n")
cat("Done visualization\n")
cat("=====================================\n")
q("no")
