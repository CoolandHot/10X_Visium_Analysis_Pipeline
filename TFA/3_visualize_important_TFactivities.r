library(ggplot2)
raw.data.dir <- "/app/data/"
data.dir <- "/app/data/TFA"
setwd(data.dir)
output.dir <- paste0(data.dir, "/output/")
cell_types <- read.csv(paste0(output.dir, "ordered_cell_types.csv"), row.names = 1)

config <- yaml::read_yaml(paste0(raw.data.dir, "config/batch_config.yaml"))
batch_names <- config$batch_names

for (batch_id in batch_names) {
    TFA <- read.csv(paste0(output.dir, batch_id, "_TF_activities.csv"), row.names = 1)
    TFA$cell_types <- cell_types$x[match(TFA$cell_types, rownames(cell_types))]
    TFA <- TFA[!is.na(TFA$cell_types), ] # Filter out rows with NA cell types

    top_features_df <- read.csv(paste0(output.dir, batch_id, "_top_features_by_cell_type_RandomForest.csv"))
    top10_per_celltype <- top_features_df |>
        dplyr::mutate(feature = make.names(feature)) |>
        dplyr::group_by(cell_type) |>
        dplyr::slice_head(n = 10) |>
        dplyr::ungroup()

    selected_TF_Z <- dplyr::select(TFA, top10_per_celltype$feature)
    selected_TF_Z <- apply(selected_TF_Z, 2, function(x) (x - min(x)) / (max(x) - min(x)))

    ## Build the dataframe for heatmap
    data_fr <- as.data.frame(selected_TF_Z)
    data_fr$cells <- rownames(selected_TF_Z)
    data_fr$cell_type <- TFA$cell_types

    data_fr_long <- tidyr::pivot_longer(data_fr,
        cols = -c(cells, cell_type),
        names_to = "TF",
        values_to = "activity"
    )


    # the cell types & genes are too many, split into 4 cell types & 40 genes each plot
    unique_cell_types <- unique(data_fr_long$cell_type)
    cell_type_groups <- split(unique_cell_types, ceiling(seq_along(unique_cell_types) / 4))
    unique_TFs <- unique(data_fr_long$TF)
    tf_groups <- split(unique_TFs, ceiling(seq_along(unique_TFs) / 40))

    plots <- list()
    for (cell_group_index in seq_along(cell_type_groups)) {
        cell_group <- cell_type_groups[[cell_group_index]]

        # Set the order of cell_type based on the current cell group
        filtered_data <- data_fr_long |>
            dplyr::filter(cell_type %in% cell_group)
        filtered_data$cell_type <- factor(filtered_data$cell_type, levels = cell_group)

        for (tf_group_index in seq_along(tf_groups)) {
            tf_group <- tf_groups[[tf_group_index]]

            # Further filter data for the current TF group
            filtered_data_group <- filtered_data |>
                dplyr::filter(TF %in% tf_group)

            # Generate heatmap
            plots[[length(plots) + 1]] <- ggplot(filtered_data_group, mapping = aes(x = cells, y = TF, fill = activity)) +
                ggtitle(paste0("plot ", tf_group_index, " for cell types:"),
                    subtitle = paste(cell_group, collapse = ", ")
                ) +
                geom_tile() +
                scale_fill_gradient(low = "white", high = "red", name = "Activity") +
                scale_y_discrete(position = "right") +
                facet_grid(cols = vars(cell_type), drop = TRUE, space = "free", scales = "free") +
                theme(
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                    axis.title.y = element_blank(),
                    axis.text.y = element_text(size = 12, face = "bold"),
                    strip.text = element_text(size = 10)
                )
        }
    }

    ggpubr::ggexport(
        filename = paste0(output.dir, batch_id, "_Top10_TFA_heatmap_RandomForest.pdf"),
        plotlist = plots
    )
}


cat("========================\n")
cat("done\n")
cat("========================\n")
quit("no")
