source("util_headers.r")


merged_obj <- readRDS(paste0(rds_data_dir, output.file.prefix, "_merged.rds"))

# ===========================================
# load cell type prediction results
# ===========================================
cellType_1_cluster_df <- read.csv(paste0(output_dirs$clustering, "predict_GBM_scSeq_cellType_1.csv"), stringsAsFactors = FALSE, check.names = FALSE)
rownames(cellType_1_cluster_df) <- cellType_1_cluster_df$Cell
cellType_1_cluster_df$Cell <- NULL

cellType_2_cluster_df <- read.csv(paste0(output_dirs$clustering, "predict_GBM_scSeq_cellType_2.csv"), stringsAsFactors = FALSE, check.names = FALSE)
rownames(cellType_2_cluster_df) <- cellType_2_cluster_df$Cell
cellType_2_cluster_df$Cell <- NULL

merged_obj$predict_GBM_scSeq_cellType_1 <- cellType_1_cluster_df
merged_obj$predict_GBM_scSeq_cellType_2 <- cellType_2_cluster_df

merged_obj <- JoinLayers(merged_obj)
# ===========================================


merged_obj_list <- SplitObject(merged_obj, split.by = "batch")


# ===========================================
# predicted cell types under each cluster
# export csv & visualize
# ===========================================
export_csv_and_plot <- function(aggregated_data, cellType_cat, batch_id) {
  aggregated_data |>
    write.csv(paste0(output_dirs$cell_type_spacexr, "spaceXR_", cellType_cat, "_", batch_id, ".csv"),
      row.names = T
    )

  clean_data <- aggregated_data[!is.na(rownames(aggregated_data)), !is.na(colnames(aggregated_data))] |> as.data.frame.matrix()
  all_cell_types <- colnames(clean_data)
  color_scale <- setNames(custom_colors[1:length(all_cell_types)], all_cell_types)
  clean_data$cluster <- rownames(clean_data)

  ##############################################
  # area plot
  ##############################################
  pie_data <- clean_data |>
    tidyr::pivot_longer(-cluster, names_to = "cell_type", values_to = "count") |>
    dplyr::group_by(cluster) |>
    dplyr::mutate(percentage = count / sum(count))

  area_plot <- ggplot(pie_data, aes(x = cluster, y = percentage, fill = cell_type)) +
    geom_bar(stat = "identity", position = "fill") +
    labs(
      x = "Cluster",
      y = "Relative Percentage",
      fill = "Cell Type",
      title = paste(cellType_cat, " percentages on each cluster in", batch_id)
    ) +
    scale_fill_manual(values = color_scale) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  return(area_plot)
}


plot_list <- list()

# ===========================================
# export predicted cell types under each cluster csv
# ===========================================
for (batch_id in batch_names) {
  one_subset <- subset(merged_obj_list[[batch_id]],
    cells = Cells(merged_obj_list[[batch_id]][[ifelse(VisiumHD, "Spatial.008um", "Spatial")]])
  )

  for (cellType_cat in c(
    "predict_GBM_scSeq_cellType_1",
    "predict_GBM_scSeq_cellType_2"
  )) {
    meta_data <- one_subset[[]]
    aggregated_data <- table(meta_data[[cellType_cat]],
      meta_data[[valid_cluster_method]],
      useNA = "ifany"
    ) |> t()

    plot_list[[length(plot_list) + 1]] <- export_csv_and_plot(aggregated_data, cellType_cat, batch_id)
  }
}


ggpubr::ggexport(
  filename = paste0(output_dirs$cell_type_spacexr, "spaceXR_predict_cellType_proportion.pdf"),
  plotlist = plot_list
)

quit("no")
