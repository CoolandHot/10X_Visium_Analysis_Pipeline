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

# Generate and save a spatial plot as a PDF.
# Args:
#   seurat_obj: A Seurat object containing the data.
#   batch_id: Identifier for the batch being plotted.
#   marker: Gene marker to visualize.
#   group_by: Metadata column to group by.
#   output_dir: Directory where the output file will be saved.
#   color_map: Color map for the clusters.
cluster_labelling_env$pdf_spatial_plot <- function(seurat_obj, batch_id, marker, group_by, output_dir, color_map) {
    marker_expr_plot <- Seurat::SpatialFeaturePlot(
        seurat_obj,
        features = marker,
        image.alpha = 0.7
    ) + ggplot2::ggtitle(paste(marker, "Expression"))

    cluster_plot <- Seurat::SpatialDimPlot(
        seurat_obj,
        group.by = group_by,
        cols = color_map,
        image.alpha = 0.7
    ) + ggplot2::ggtitle("Clusters")

    combined_plot <- cowplot::plot_grid(marker_expr_plot, cluster_plot, nrow = 1)

    ggplot2::ggsave(
        filename = paste0(output_dir, "spatial_heatmap_", batch_id, "_", marker, ".pdf"),
        plot = combined_plot,
        width = 15,
        height = 7
    )
}

# Generate and save an interactive spatial plot as an HTML file using Plotly.
# Args:
#   seurat_obj: A Seurat object containing the data.
#   batch_id: Identifier for the batch being plotted.
#   marker: Gene marker to visualize.
#   group_by: Metadata column to group by.
#   output_dir: Directory where the output file will be saved.
#   color_map: Color map for the clusters.
cluster_labelling_env$plotly_spatial_plot <- function(seurat_obj, batch_id, marker, group_by, output_dir, color_map) {
    marker_expr <- Seurat::GetAssayData(seurat_obj, layer = "data")[grep(paste0("^", marker, "$"), rownames(seurat_obj), ignore.case = TRUE), ]
    st_location <- Seurat::GetTissueCoordinates(seurat_obj)[names(marker_expr), ]
    st_location$cell <- NULL

    # Plot marker expression
    plot1 <- plotly::plot_ly(
        x = st_location[, 2], y = -st_location[, 1],
        type = "scatter", mode = "markers",
        color = marker_expr,
        marker = list(size = 5)
    ) |>
        plotly::layout(
            paper_bgcolor = "black",
            plot_bgcolor = "black",
            font = list(color = "white")
        )

    # Plot cluster identities
    idents <- Seurat::Idents(seurat_obj) |> na.omit()
    st_location <- Seurat::GetTissueCoordinates(seurat_obj)[names(idents), ]
    st_location$cell <- NULL
    plot2 <- plotly::plot_ly(
        x = st_location[, 2], y = -st_location[, 1],
        type = "scatter", mode = "markers",
        color = idents,
        colors = color_map,
        marker = list(size = 3)
    ) |>
        plotly::layout(
            paper_bgcolor = "black",
            plot_bgcolor = "black",
            font = list(color = "white")
        )

    # Combine plots
    combined_plot <- plotly::subplot(
        plot1 |> plotly::colorbar(title = paste(marker, "Expression")),
        plot2,
        nrows = 1,
        shareX = TRUE,
        shareY = TRUE
    ) |>
        plotly::layout(
            legend = list(orientation = "h", y = -0.2, x = 0.5, xanchor = "center", yanchor = "top"),
            annotations = list(
                list(x = 0.25, y = 1, text = paste(marker, "Expression"), showarrow = FALSE, xref = "paper", yref = "paper", font = list(color = "white")),
                list(x = 0.75, y = 1, text = "Clusters", showarrow = FALSE, xref = "paper", yref = "paper", font = list(color = "white"))
            )
        )

    # Save as HTML
    htmlwidgets::saveWidget(
        combined_plot,
        file = paste0(output_dir, "clusters_spatial_", batch_id, "_", marker, ".html")
    )
}

# Generate and save UMAP plots (2D and 3D) for a Seurat object.
# Args:
#   seurat_obj: A Seurat object containing the data.
#   reduction: Dimensionality reduction method (e.g., "umap").
#   group_by: Metadata column to group by.
#   file_name: Name of the output file.
#   plot_title: Title of the plot.
#   output_dir: Directory where the output file will be saved.
#   color_map: Color map for the clusters.
cluster_labelling_env$plot_umap <- function(seurat_obj, reduction, group_by, file_name, plot_title, output_dir, color_map) {
    if (!reduction %in% names(seurat_obj@reductions)) {
        seurat_obj <- run_dim_reduction(seurat_obj, reduction)
    }

    # 2D UMAP
    p <- Seurat::DimPlot(
        seurat_obj,
        reduction = reduction,
        group.by = group_by,
        label = TRUE,
        cols = color_map
    ) + Seurat::NoAxes() +
        ggplot2::ggtitle(plot_title) +
        ggplot2::theme_minimal()

    ggplot2::ggsave(
        filename = paste0(output_dir, file_name, "_umap.pdf"),
        plot = p,
        width = 10,
        height = 10
    )

    # 3D UMAP
    umap_data <- Seurat::Embeddings(seurat_obj, reduction = reduction)
    clusters <- seurat_obj[[group_by, drop = TRUE]]
    plot_data <- data.frame(
        UMAP_1 = umap_data[, 1],
        UMAP_2 = umap_data[, 2],
        UMAP_3 = umap_data[, 3],
        cluster = clusters
    )

    plotly::plot_ly(
        data = plot_data,
        x = ~UMAP_1,
        y = ~UMAP_2,
        z = ~UMAP_3,
        color = ~cluster,
        colors = color_map,
        type = "scatter3d",
        mode = "markers",
        marker = list(size = 3)
    ) |>
        plotly::layout(title = file_name) |>
        htmlwidgets::saveWidget(
            file = paste0(output_dir, file_name, "_umap_3d.html"),
            selfcontained = TRUE
        )
}

# Generate spatial cluster plots with cell percentage annotations.
# Args:
#   seurat_obj: A Seurat object containing the data.
#   group_by: Metadata column to group by.
#   file_name: Name of the output file.
#   plot_title: Title of the plot.
#   output_dir: Directory where the output file will be saved.
#   color_map: Color map for the clusters.
#   total_cells_per_batch: Total number of cells per batch.
#   batch_names: Names of the batches.
cluster_labelling_env$plot_spatial_clusters_with_cell_percentage <- function(seurat_obj, group_by, file_name, plot_title, output_dir, color_map, total_cells_per_batch, batch_names) {
    levels(seurat_obj@meta.data[[group_by]]) <- c(levels(seurat_obj@meta.data[[group_by]]), "Unknown")
    seurat_obj@meta.data[[group_by]][is.na(seurat_obj@meta.data[[group_by]])] <- "Unknown"

    spatial_plots <- Seurat::SpatialDimPlot(seurat_obj,
        group.by = group_by,
        cols = color_map, alpha = 1, image.alpha = 0.7,
        label = TRUE, repel = TRUE,
        label.size = 4, pt.size.factor = ifelse(VisiumHD, 4, 1.6)
    )

    spatial_plots_list <- lapply(seq_along(spatial_plots), function(i) {
        current_batch <- batch_names[[i]]
        slice_cells <- Seurat::WhichCells(seurat_obj, expression = batch == current_batch)
        slice_clusters <- seurat_obj@meta.data[slice_cells, group_by]
        cluster_counts <- table(slice_clusters)
        cluster_counts <- cluster_counts[order(as.numeric(names(cluster_counts)))]
        count_text <- "Cells per cluster:\n"
        for (j in seq_along(cluster_counts)) {
            count_text <- paste0(
                count_text, names(cluster_counts)[j], ": ",
                cluster_counts[j], " (", sprintf("%.1f%%", 100 * cluster_counts[j] / total_cells_per_batch[[current_batch]]), ")\n"
            )
        }
        spatial_plots[[i]] +
            ggplot2::ggtitle(paste0(plot_title, " in ", current_batch)) +
            ggplot2::annotation_custom(
                grid::textGrob(count_text,
                    x = 0.95, y = 0.95,
                    hjust = 1, vjust = 1,
                    gp = grid::gpar(fontsize = 10)
                )
            ) +
            ggplot2::labs(fill = "Cluster IDs")
    })
    names(spatial_plots_list) <- paste0("Slice_", seq_along(spatial_plots))

    ggpubr::ggexport(
        plotlist = spatial_plots_list,
        filename = paste0(output_dir, file_name, "_spatial.pdf"),
        width = 10, height = 10
    )
}

# Determine the most likely cell type based on voting from possible cell types.
# Args:
#   possible_cell_types: A vector of possible cell types for a cluster.
# Returns:
#   The most likely cell type or a combination of cell types if tied.
cluster_labelling_env$determine_cell_type_by_votes <- function(possible_cell_types) {
    cell_type_counts <- table(unlist(strsplit(possible_cell_types, ", ")))
    max_count <- max(cell_type_counts)
    most_common_cell_types <- names(cell_type_counts[cell_type_counts == max_count])
    if (length(most_common_cell_types) > 1) {
        return(paste(most_common_cell_types, collapse = "/"))
    } else {
        return(most_common_cell_types)
    }
}
