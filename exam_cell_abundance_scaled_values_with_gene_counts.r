# ===================================================
# Plot line plots using ggplot2 with two y-axes
library(ggplot2)
library(reshape2)
require(Matrix)

scdata <- readRDS("rds_data/GBM_1st_four_merged_spatial_forTrailmaker.rds")
gene_counts <- colSums(scdata@assays$RNA@counts)
cell_abundance <- read.csv("output/cell_type_prediction/cell2location_results/cell_abundances_and_clusters.csv", row.names = 1) |> dplyr::select(-c("batch", "region_cluster"))
cell_counts <- rowSums(cell_abundance)
names(cell_counts) <- rownames(cell_abundance)

# Select first 100 matched rownames
matched_names <- intersect(names(gene_counts), names(cell_counts))
matched_names <- matched_names[1:min(100, length(matched_names))]
gene_counts_matched <- gene_counts[matched_names]
cell_counts_matched <- cell_counts[matched_names]

# export the matched_names of gene_counts and cell_abundance to two separate csv files
write.csv(scdata@assays$RNA@counts[, matched_names] |> t(),
    file = "output/gene_counts_matched.csv",
    row.names = TRUE, quote = FALSE
)
write.csv(gene_counts_matched,
    file = "output/gene_counts_sum_matched.csv",
    row.names = TRUE, quote = FALSE
)
write.csv(cell_abundance[matched_names, ],
    file = "output/cell_abundance_matched.csv",
    row.names = TRUE, quote = FALSE
)
write.csv(cell_counts_matched,
    file = "output/cell_abundance_sum_matched.csv",
    row.names = TRUE, quote = FALSE
)


df_plot <- data.frame(
    Cell_Index = seq_along(gene_counts_matched),
    Gene_Counts = as.numeric(gene_counts_matched),
    Cell_Counts = as.numeric(cell_counts_matched)
)

# Scale Cell_Counts to match Gene_Counts range for plotting
scale_factor <- max(df_plot$Gene_Counts, na.rm = TRUE) / max(df_plot$Cell_Counts, na.rm = TRUE)
df_plot$Cell_Counts_scaled <- df_plot$Cell_Counts * scale_factor

ggplot(df_plot, aes(x = Cell_Index)) +
    geom_line(aes(y = Gene_Counts, color = "Gene_Counts"), size = 1.2) +
    geom_line(aes(y = Cell_Counts_scaled, color = "Cell_Counts"), size = 1.2) +
    scale_y_continuous(
        name = "Gene Counts",
        sec.axis = sec_axis(~ . / scale_factor, name = "Cell Counts")
    ) +
    scale_color_manual(values = c("Gene_Counts" = "blue", "Cell_Counts" = "red")) +
    labs(
        title = "Gene Counts and Cell Counts (first 100 matched)",
        x = "Cell Index",
        color = "Type"
    ) +
    theme_minimal()
