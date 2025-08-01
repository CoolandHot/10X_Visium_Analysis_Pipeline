source("util_headers.r")
# ===============================================================================
# CELL TYPE ANALYSIS AND VISUALIZATION
# ===============================================================================
# Load cell type abundance and cluster assignments
cell_type_abundance_path <- paste0(output_dirs$cell_type_cell2loc, "cell_abundances_and_clusters.csv")
cell_type_df <- read.csv(cell_type_abundance_path, row.names = 1)

cluster_csv_path <- paste0(output_dirs$clustering, cluster_method, ".csv")
if (file.exists(cluster_csv_path)) {
    cluster_data <- read.csv(cluster_csv_path, row.names = 1)
    cluster_vec <- setNames(cluster_data[[1]], rownames(cluster_data))
} else {
    stop(paste("Cluster method", cluster_method, "not found in CSV."))
}

# Add region assignments to cell type dataframe
cell_type_df$region <- dplyr::case_when(
    cluster_data[rownames(cell_type_df), 1] %in% inTumour_cluster_nums_vector ~ "inTumour",
    cluster_data[rownames(cell_type_df), 1] %in% outTumour_cluster_nums_vector ~ "outTumour",
    cluster_data[rownames(cell_type_df), 1] %in% edgeTumour_cluster_nums_vector ~ "edgeTumour",
    TRUE ~ NA_character_
)

# Export cell type data for each region as a whole
for (region in c("inTumour", "outTumour", "edgeTumour")) {
    region_df <- cell_type_df[cell_type_df$region == region, ]
    write.csv(
        region_df,
        file = paste0(output_dirs$deg_celltypes, region, "_cellTypes.csv"),
        row.names = TRUE
    )
}

# Handle combined cell types if available
cell_type_abundance_combined_path <- paste0(output_dirs$cell_type_cell2loc, "cell_abundances_and_clusters_combined.csv")
if (file.exists(cell_type_abundance_combined_path)) {
    combine_celltype <- read.csv(cell_type_abundance_combined_path, row.names = 1)
    combine_celltype$region <- dplyr::case_when(
        cluster_data[rownames(combine_celltype), 1] %in% inTumour_cluster_nums_vector ~ "inTumour",
        cluster_data[rownames(combine_celltype), 1] %in% outTumour_cluster_nums_vector ~ "outTumour",
        cluster_data[rownames(combine_celltype), 1] %in% edgeTumour_cluster_nums_vector ~ "edgeTumour",
        TRUE ~ NA_character_
    )
    for (region in c("inTumour", "outTumour", "edgeTumour")) {
        region_combined <- combine_celltype[combine_celltype$region == region, ]
        write.csv(
            region_combined,
            file = paste0(output_dirs$deg_combine_celltypes, region, "_cellTypes_combined.csv"),
            row.names = TRUE
        )
    }
}

# ===============================================================================
# HELPER FUNCTIONS
# ===============================================================================

# Function to prepare proportion data for visualization
prepare_proportion_data <- function(cell_type_data, celltype_columns) {
    batch_region_props_list <- list()
    for (batch in unique(cell_type_data$batch)) {
        batch_df <- cell_type_data[cell_type_data$batch == batch, ]
        region_props <- lapply(
            c("inTumour", "outTumour", "edgeTumour"),
            function(region) {
                df <- batch_df[batch_df$region == region, celltype_columns, drop = FALSE]
                if (nrow(df) == 0) {
                    return(NULL)
                }
                colSums(df, na.rm = TRUE)
            }
        )
        names(region_props) <- c("inTumour", "outTumour", "edgeTumour")
        region_props_df <- do.call(rbind, region_props)
        if (!is.null(region_props_df)) {
            region_props_norm_df <- sweep(region_props_df, 1, rowSums(region_props_df, na.rm = TRUE), FUN = "/") |> as.data.frame()
            region_props_norm_df$region <- rownames(region_props_norm_df)
            region_props_norm_df$batch <- batch
            batch_region_props_list[[batch]] <- region_props_norm_df
        }
    }
    all_region_props_norm_df <- do.call(rbind, batch_region_props_list)
    region_props_long <- reshape2::melt(all_region_props_norm_df, id.vars = c("batch", "region"), variable.name = "cell_type", value.name = "proportion")
    # Filter out rows with zero or NA proportions to avoid warnings
    region_props_long[!is.na(region_props_long$proportion) & region_props_long$proportion > 0, ]
}

# Function to create static PDF plots
create_static_plots <- function(region_props_long, output_file, title_prefix) {
    pdf(file = output_file, width = 8, height = 5)
    for (b in unique(region_props_long$batch)) {
        batch_data <- region_props_long[region_props_long$batch == b, ]
        if (nrow(batch_data) == 0) next

        p <- ggplot(batch_data, aes(x = region, y = proportion, fill = cell_type)) +
            geom_bar(stat = "identity", position = "fill") +
            labs(title = paste(title_prefix, "-", b), x = "Region", y = "Proportion") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            scale_fill_manual(values = rainbow(length(unique(batch_data$cell_type))))
        print(p)
    }
    dev.off()
}

# Function to create interactive plotly plots
create_interactive_plot <- function(region_props_long, output_file, title) {
    plotly_fig <- plot_ly(
        region_props_long,
        x = ~region,
        y = ~proportion,
        color = ~cell_type,
        type = "bar",
        text = ~ paste(cell_type, "<br>Proportion:", round(proportion, 3)),
        hovertemplate = "%{text}<extra></extra>",
        transforms = list(
            list(
                type = "filter",
                target = ~batch,
                operation = "=",
                value = unique(region_props_long$batch)[1]
            )
        )
    ) |>
        layout(
            barmode = "stack",
            title = title,
            xaxis = list(title = "Region"),
            yaxis = list(title = "Proportion"),
            updatemenus = list(
                list(
                    type = "dropdown",
                    active = 0,
                    x = 1.15,
                    y = 1,
                    buttons = lapply(
                        unique(region_props_long$batch),
                        function(b) {
                            list(
                                method = "restyle",
                                args = list("transforms[0].value", b),
                                label = paste("Batch:", b)
                            )
                        }
                    )
                )
            )
        )
    htmlwidgets::saveWidget(plotly_fig, file = output_file, selfcontained = TRUE)
}

# ===============================================================================
# CELL TYPE PROPORTION VISUALIZATION
# ===============================================================================

library(ggplot2)
library(reshape2)
library(plotly)

config_yaml <- yaml::read_yaml("./cellType/config.yaml")
combine_cell_types <- config_yaml$shared$combine_cell_types
cell_type_combinations <- config_yaml$shared$cell_type_combinations

# If combine_cell_types is TRUE, aggregate subtypes into combined cell types
if (isTRUE(combine_cell_types)) {
    # For each combination, sum the columns and create combined columns
    combined_mat <- sapply(names(cell_type_combinations), function(comb_name) {
        subtypes <- cell_type_combinations[[comb_name]]
        valid_subtypes <- intersect(subtypes, colnames(cell_type_df))
        if (length(valid_subtypes) == 0) {
            rep(NA, nrow(cell_type_df))
        } else {
            rowSums(cell_type_df[, valid_subtypes, drop = FALSE], na.rm = TRUE)
        }
    })
    combined_df <- as.data.frame(combined_mat, row.names = rownames(cell_type_df))
    # Keep subtypes not included in any combination
    all_combined_subtypes <- unique(unlist(cell_type_combinations))
    remaining_subtypes <- setdiff(colnames(cell_type_df), c(all_combined_subtypes, "region", "cluster", "batch", "region_cluster"))
    # Keep region, cluster, batch, region_cluster columns if present
    meta_cols <- intersect(c("region", "cluster", "batch", "region_cluster"), colnames(cell_type_df))
    cell_type_df_combined <- cbind(
        combined_df,
        cell_type_df[, remaining_subtypes, drop = FALSE],
        cell_type_df[, meta_cols, drop = FALSE]
    )
}

# --- Original cell type visualization ---
celltype_cols <- setdiff(colnames(cell_type_df), c("region", "cluster", "batch", "region_cluster"))
region_props_long <- prepare_proportion_data(cell_type_df, celltype_cols)

create_static_plots(
    region_props_long,
    paste0(output_dirs$deg_celltypes, "cellType_proportion_by_region_by_batch.pdf"),
    "Cell type proportion by region"
)

create_interactive_plot(
    region_props_long,
    paste0(output_dirs$deg_celltypes, "cellType_proportion_by_region_by_batch.html"),
    "Cell type proportion by region (split by batch)"
)

# Save original cell type proportion data to CSV
write.csv(
    region_props_long,
    file = paste0(output_dirs$deg_celltypes, "cellType_proportion_data.csv"),
    row.names = FALSE
)

# --- Combined cell type visualization ---
if (isTRUE(combine_cell_types)) {
    meta_cols <- intersect(c("region", "cluster", "batch", "region_cluster"), colnames(cell_type_df))
    celltype_combined_cols <- setdiff(colnames(cell_type_df_combined), meta_cols)
    region_props_combined_long <- prepare_proportion_data(cell_type_df_combined, celltype_combined_cols)

    create_static_plots(
        region_props_combined_long,
        paste0(output_dirs$deg_celltypes, "cellTypeCombined_proportion_by_region_by_batch.pdf"),
        "Combined cell type proportion by region"
    )

    create_interactive_plot(
        region_props_combined_long,
        paste0(output_dirs$deg_celltypes, "cellTypeCombined_proportion_by_region_by_batch.html"),
        "Combined cell type proportion by region (split by batch)"
    )

    # Save combined cell type proportion data to CSV
    write.csv(
        region_props_combined_long,
        file = paste0(output_dirs$deg_celltypes, "cellTypeCombined_proportion_data.csv"),
        row.names = FALSE
    )
}

quit("no")
