# Only works on Seurat5!!!!


# spatialTrans_data_prefix <- "/vol/research/scratch1/NOBACKUP/hh01116/10X_Genomics_KD_data/"
spatialTrans_data_prefix <- "/app/data/"
output.dir <- paste0(spatialTrans_data_prefix, "cellType_deconvolute/output/")


library(ggplot2)
custom_colors <- c(
  "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f",
  "#e5c494", "#b3b3b3", "#1b9e77", "#d95f02", "#7570b3", "#e7298a",
  "#66a61e", "#e6ab02", "#a6761d", "#666666", "#a6cee3", "#1f78b4",
  "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
  "#cab2d6", "#6a3d9a", "#ffff99", "#b15928"
)

## interested cell types
immune_cells <- c(
  "Microglial cells", "Macrophages", "NK cells", "Dendritic cells",
  "T cells", "Neutrophils", "B cells", "Mast cells", "Monocytes"
)


make_areaPlot_spatialPiePlot <- function(aggregated_data, spatial_plot, batch) {
  # returns a list of two plots, $area and $spatial

  # Create a consistent color scale for all cell types
  all_cell_types <- colnames(aggregated_data)[-1]
  color_scale <- setNames(custom_colors[1:length(all_cell_types)], all_cell_types)

  ##############################################
  # area plot
  ##############################################
  pie_data <- aggregated_data |>
    tidyr::pivot_longer(-cluster, names_to = "cell_type", values_to = "count") |>
    dplyr::group_by(cluster) |>
    dplyr::mutate(percentage = count / sum(count))

  area_plot <- ggplot(pie_data, aes(x = cluster, y = percentage, fill = cell_type)) +
    geom_bar(stat = "identity", position = "fill") +
    labs(
      x = "Cluster", y = "Relative Percentage", fill = "Cell Type",
      title = paste("Cell type percentages on each cluster in", batch)
    ) +
    scale_fill_manual(values = color_scale) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ##############################################
  # pie chart on top of spatial map
  ##############################################
  # Step 1: Create the base spatial plot
    #   spatial_plot <- Seurat::SpatialDimPlot(st_object,
    #     group.by = "seurat_clusters",
    #     label = TRUE,
    #     alpha = .6
    #   ) + Seurat::NoLegend()

  # Step 2: Extract the cluster labels coordinates
  label_positions <- ggplot_build(spatial_plot)$data[[2]] |> dplyr::select(label, x, y)

  # Function to create a pie chart grob
  create_pie <- function(cluster_label) {
    cluster_data <- pie_data[pie_data$cluster == cluster_label, -1]
    if(nrow(cluster_data)!=0){
        pie_chart <- ggplot(cluster_data, aes(x = "", y = count, fill = cell_type)) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar("y") +
        theme_void() +
        theme(legend.position = "none") +
        scale_fill_manual(values = color_scale) +
        geom_text(aes(x = 1, y = sum(count) / 2, label = cluster_label), size = 6, vjust = -0.5) # add cluster label on pie chart

        return(pie_chart)
    }

  }

  # Step 3: Overlay pie charts onto the spatial plot
  base_plot <- spatial_plot + labs(title = paste("Cell type percentages on each cluster in", batch))
  for (i in 1:nrow(label_positions)) {
    cluster <- label_positions$label[i]
    x <- label_positions$x[i]
    y <- label_positions$y[i]

    pie_grob <- create_pie(cluster)

    if(!is.null(pie_grob)){
        base_plot <- base_plot + annotation_custom(
        grob = ggplotGrob(pie_grob),
        xmin = x - 30, xmax = x + 30,
        ymin = y - 30, ymax = y + 30
        )
    }
  }

  # Step 4: add a legend
  # Create a dummy legend grob
  legend_plot <- ggplot(
    data.frame(
      cell_type = all_cell_types,
      dummy = 1
    ),
    aes(x = 1, y = dummy, fill = cell_type)
  ) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = color_scale, name = "Cell Type") +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 6),
      legend.direction = "horizontal"
    )

  legend_grob <- cowplot::get_plot_component(legend_plot, "guide-box-bottom", return_all = TRUE)

  # Add the legend to the base plot
  final_plot <- cowplot::ggdraw(base_plot) + cowplot::draw_grob(legend_grob, y = -0.4)

  return(list(area = area_plot, spatial = final_plot))
}


merged_obj <- readRDS(paste0(spatialTrans_data_prefix,"rds_data/",  "four_merged_clustered.rds"))
merged_obj_list <- Seurat::SplitObject(merged_obj, split.by="batch")
CARD_prop_list <- readRDS(paste0(spatialTrans_data_prefix, "cellType_deconvolute/deconvoluted/CARD_cellType_prop_fourBatches.Rds"))


# ================
for(batch in c("SAL", "KD", "RAD", "KD_RAD")){
    st_object <- merged_obj_list[[batch]]
    Seurat::DefaultAssay(st_object) <- "Spatial.008um"
    spatial_plot <- Seurat::SpatialDimPlot(
        st_object,
        group.by = "projected_seurat_cluster",
        label = TRUE,
        alpha = .6
    ) + Seurat::NoLegend()

    n_cells_keep <- st_object[['Spatial.008um']]$counts |> ncol()
    Seurat::Idents(st_object) <- "projected_seurat_cluster"
    clusters <- Seurat::Idents(st_object)[1:n_cells_keep]

    # extract the cluster numbers and combine with cell type proportions as a dataframe
    # A tibble: 2,187 × 19
    #    cell_id            cluster `Ciliated cells` `Oligodendrocyte precursor cells`
    #    <chr>              <chr>              <dbl>                             <dbl>
    #  1 AAACAAGTATCTCCCA-1 0                 0.0976                           0.0345
    #  2 AAACAGAGCGACTCCT-1 3                 0.0536                           0.00968
    #  3 AAACCCGAACGAAATC-1 0                 0.103                            0.0262
    result_df <- tibble::as_tibble(CARD_prop_list[[batch]], rownames = "cell_id") |>
    dplyr::left_join(
        tibble::tibble(
        cell_id = gsub(paste0("^", batch,"_"), "", names(clusters)),
        cluster = as.vector(clusters)
        ),
        by = "cell_id"
    ) |>
    dplyr::select(cell_id, cluster, dplyr::everything())

    if(exists("immune_cells")){
    # if defined immune cell types, filter to only immune cells
    result_df <- dplyr::select(result_df, cell_id, cluster, dplyr::all_of(immune_cells))
    }

    aggregated_data <- result_df |>
    dplyr::select(-cell_id) |>
    dplyr::group_by(cluster) |>
    dplyr::summarise(dplyr::across(dplyr::everything(), sum))

    write.csv(aggregated_data, 
    file = paste0(output.dir, "GBM_HD_cellType_proportion_on_clusters_", batch, ".csv"),
    row.names=FALSE
    )

    result_plot_list <- make_areaPlot_spatialPiePlot(aggregated_data, spatial_plot, batch)

    ggpubr::ggexport(
    filename = paste0(output.dir, "GBM_HD_cellType_proportion_on_clusters_", batch, ".pdf"),
    plotlist = result_plot_list
    )
}