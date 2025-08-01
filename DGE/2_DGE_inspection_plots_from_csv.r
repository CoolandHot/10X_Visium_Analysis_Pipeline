# ===============================================================================
# CLUSTER DGE INSPECTION PLOTS SCRIPT
# ===============================================================================
# Purpose:
#   This script generates inspection plots for merged differential gene expression
#   (DGE) results on clusters, cell subtypes, and optionally combined cell types.
#   It is intended to be run after generating merged DGE CSV files.
#
# What are "Inspection Plots"?
#   In this context, "inspection plots" are a set of diagnostic and summary visualizations
#   designed to help users quickly assess, compare, and interpret the results of DGE analyses.
#   These plots provide visual summaries of gene expression changes across different groups
#   or comparisons, highlight the most significant genes, and reveal patterns or outliers
#   that may warrant further investigation.
#
#   Specifically, this script produces:
#     - Scatterplots: Compare log2 fold changes for shared genes between all pairs of comparisons,
#       allowing visual assessment of concordance or divergence in gene regulation.
#     - Stacked Bar Charts: Show the number of up- and down-regulated genes per comparison,
#       summarizing the overall direction and magnitude of differential expression.
#     - Volcano Plots: Display the relationship between fold change and statistical significance
#       for each comparison, highlighting the most strongly regulated genes.
#
# Usage:
#   - Source this script in R after running the DGE analysis and merging steps.
#   - The script expects merged DGE CSV files to be present in the output directories.
#   - Plots and intermediate data will be saved in subdirectories of the output folders.
#
# Main Features:
#   1. Standardizes DGE data columns from different formats.
#   2. Generates scatterplots comparing log2 fold changes for all pairs of comparisons.
#   3. Produces stacked bar charts showing the number of up/down-regulated genes per comparison.
#   4. Creates volcano plots for each comparison, highlighting top up- and down-regulated genes.
#   5. Handles DGE results for clusters, cell types, and optionally combined cell types.
#
# Output:
#   - PDF plots and CSV data files are saved in 'inspection_plots_pdf' and
#     'inspection_plot_data' directories under each DGE output folder.
#
# Dependencies:
#   - Requires ggplot2, ggrepel, dplyr, gtools, and RColorBrewer packages.
#   - Assumes 'util_headers.r' is available for loading project-specific variables.
#
# Author: [Your Name or Lab]
# Date: [YYYY-MM-DD]
# ===============================================================================

# Inspection plots for merged DGE on clusters
# Usage: source this script after generating merged_dge_on_clusters.csv
source("util_headers.r")

library(ggplot2)
library(RColorBrewer)

cluster_dge_inspect_env <- new.env()

# Standardize column names for different DGE file formats
cluster_dge_inspect_env$standardize_dge_columns <- function(dge) {
    # Check if this is the cellTypes format
    if ("names" %in% colnames(dge) && "logfoldchanges" %in% colnames(dge)) {
        # Map cellTypes format to standard format
        dge <- dge |>
            dplyr::rename(
                gene = names,
                avg_log2FC = logfoldchanges,
                p_val_adj = pvals_adj,
                p_val = pvals
            )
        # Add missing columns with default values if they don't exist
        if (!"pct.1" %in% colnames(dge)) dge$pct.1 <- NA
        if (!"pct.2" %in% colnames(dge)) dge$pct.2 <- NA
        # Rename cell_type to region for consistency (if region doesn't exist)
        if ("cell_type" %in% colnames(dge) && !"region" %in% colnames(dge)) {
            dge <- dge |> dplyr::rename(region = cell_type)
        }
    }

    # Ensure avg_log2FC is numeric
    dge$avg_log2FC <- as.numeric(dge$avg_log2FC)
    dge
}


# Load merged DGE data
cluster_dge_inspect_env$load_dge_data <- function(dge_csv_path) {
    dge <- read.csv(dge_csv_path, stringsAsFactors = FALSE)
    # Standardize column names based on format
    dge <- cluster_dge_inspect_env$standardize_dge_columns(dge)
    dge
}


# 1. Scatterplots for all pairs of comparisons (shared genes)
cluster_dge_inspect_env$plot_comparison_scatterplots <- function(dge_subset, lfc_threshold = 0.25, pval_threshold = 0.05) {
    comparisons <- unique(dge_subset$comparison)
    combn_pairs <- gtools::combinations(length(comparisons), 2, comparisons)
    for (i in seq_len(nrow(combn_pairs))) {
        comp1 <- combn_pairs[i, 1]
        comp2 <- combn_pairs[i, 2]
        d1 <- dge_subset |> dplyr::filter(comparison == comp1, p_val_adj < pval_threshold)
        d2 <- dge_subset |> dplyr::filter(comparison == comp2, p_val_adj < pval_threshold)
        shared_genes <- intersect(d1$gene, d2$gene)
        if (length(shared_genes) == 0) next

        d1 <- d1 |>
            dplyr::filter(gene %in% shared_genes) |>
            dplyr::select(gene, avg_log2FC) |>
            dplyr::rename(avg_log2FC_1 = avg_log2FC) |>
            dplyr::distinct(gene, .keep_all = TRUE)

        d2 <- d2 |>
            dplyr::filter(gene %in% shared_genes) |>
            dplyr::select(gene, avg_log2FC) |>
            dplyr::rename(avg_log2FC_2 = avg_log2FC) |>
            dplyr::distinct(gene, .keep_all = TRUE)

        plot_df <- dplyr::left_join(d1, d2, by = "gene")

        plot_df$size <- abs(plot_df$avg_log2FC_1) + abs(plot_df$avg_log2FC_2)
        write.csv(plot_df, file = file.path(cluster_dge_inspect_env$plot_data_dir, paste0("scatter_", comp1, "_against_", comp2, ".csv")), row.names = FALSE)
        p <- ggplot(plot_df, aes(x = avg_log2FC_1, y = avg_log2FC_2, label = gene)) +
            ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.7) +
            ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.7) +
            ggplot2::geom_point(aes(size = size), alpha = 0.6, color = "darkred") +
            ggrepel::geom_text_repel(max.overlaps = 10, size = 2) +
            ggplot2::labs(
                title = paste("Shared DE Genes:", comp1, "vs", comp2),
                x = paste0("avg_log2FC (", comp1, ")"),
                y = paste0("avg_log2FC (", comp2, ")")
            ) +
            ggplot2::theme_minimal()
        ggplot2::ggsave(file.path(cluster_dge_inspect_env$plot_pdf_dir, paste0("scatter_", comp1, "_against_", comp2, ".pdf")), p, width = 7, height = 7)
    }
}

# 2. Stacked bar chart: number of up/down-regulated genes per comparison
cluster_dge_inspect_env$plot_stacked_barchart <- function(dge_subset, lfc_threshold = 0.25, pval_threshold = 0.05) {
    dge_subset <- dge_subset |>
        dplyr::mutate(regulation = dplyr::case_when(
            avg_log2FC >= lfc_threshold & p_val_adj < pval_threshold ~ "Up",
            avg_log2FC <= -lfc_threshold & p_val_adj < pval_threshold ~ "Down",
            TRUE ~ "NS"
        ))
    bar_df <- dge_subset |>
        dplyr::filter(regulation != "NS") |>
        dplyr::group_by(comparison, regulation) |>
        dplyr::summarise(n_genes = dplyr::n(), .groups = "drop")
    write.csv(bar_df, file = file.path(cluster_dge_inspect_env$plot_data_dir, "stacked_barchart_data.csv"), row.names = FALSE)
    p <- ggplot(bar_df, aes(x = comparison, y = n_genes, fill = regulation)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::labs(title = "Number of Up/Down-regulated Genes per Comparison", x = "Comparison", y = "Number of Genes") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    ggplot2::ggsave(file.path(cluster_dge_inspect_env$plot_pdf_dir, "stacked_barchart.pdf"), p, width = 10, height = 6)
}

# 3. Volcano plots for each comparison
cluster_dge_inspect_env$plot_volcano_plots <- function(dge_subset, lfc_threshold = 0.25, pval_threshold = 0.05) {
    comparisons <- unique(dge_subset$comparison)
    for (comp in comparisons) {
        comp_df <- dge_subset |> dplyr::filter(comparison == comp)
        comp_df <- comp_df |>
            dplyr::mutate(
                regulation = dplyr::case_when(
                    avg_log2FC >= lfc_threshold & p_val_adj < pval_threshold ~ "Up",
                    avg_log2FC <= -lfc_threshold & p_val_adj < pval_threshold ~ "Down",
                    TRUE ~ "NS"
                )
            )

        # Identify top 10 up- and down-regulated genes
        top_up <- comp_df |>
            dplyr::filter(regulation == "Up") |>
            dplyr::arrange(desc(avg_log2FC)) |>
            dplyr::slice_head(n = 10)

        top_down <- comp_df |>
            dplyr::filter(regulation == "Down") |>
            dplyr::arrange(avg_log2FC) |>
            dplyr::slice_head(n = 10)

        # Create label column
        comp_df <- comp_df |>
            dplyr::mutate(
                label = dplyr::case_when(
                    gene %in% top_up$gene ~ gene,
                    gene %in% top_down$gene ~ gene,
                    TRUE ~ ""
                )
            )

        write.csv(comp_df, file = file.path(cluster_dge_inspect_env$plot_data_dir, paste0("volcano_", comp, ".csv")), row.names = FALSE)
        p <- ggplot(comp_df, aes(x = avg_log2FC, y = -log10(p_val_adj), color = regulation)) +
            ggplot2::geom_point(alpha = 0.6) +
            ggrepel::geom_text_repel(aes(label = label), max.overlaps = 20, size = 3) +
            ggplot2::scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
            ggplot2::labs(title = paste("Volcano Plot:", comp), x = "avg_log2FC", y = "-log10(p_val_adj)") +
            ggplot2::theme_minimal()
        ggplot2::ggsave(file.path(cluster_dge_inspect_env$plot_pdf_dir, paste0("volcano_", comp, ".pdf")), p, width = 7, height = 6)
    }
}

# Main function to run all plots
cluster_dge_inspect_env$run_all_plots <- function(dge_csv_path) {
    dge <- cluster_dge_inspect_env$load_dge_data(dge_csv_path)
    # Split DGE data by group
    if (!"group" %in% colnames(dge)) {
        stop("Column 'group' not found in DGE data. Please ensure the data contains a 'group' column.")
    }

    groups <- unique(dge$group)
    cat("Found groups:", paste(groups, collapse = ", "), "\n")

    for (group_name in groups) {
        cat("Processing group:", group_name, "\n")
        dge_subset <- dge |> dplyr::filter(group == group_name)

        # Update output directories for this group
        group_plot_data_dir <- file.path(cluster_dge_inspect_env$plot_data_dir, group_name)
        group_plot_pdf_dir <- file.path(cluster_dge_inspect_env$plot_pdf_dir, group_name)
        dir.create(group_plot_data_dir, showWarnings = FALSE, recursive = TRUE)
        dir.create(group_plot_pdf_dir, showWarnings = FALSE, recursive = TRUE)

        # Temporarily update environment paths
        original_data_dir <- cluster_dge_inspect_env$plot_data_dir
        original_pdf_dir <- cluster_dge_inspect_env$plot_pdf_dir
        cluster_dge_inspect_env$plot_data_dir <- group_plot_data_dir
        cluster_dge_inspect_env$plot_pdf_dir <- group_plot_pdf_dir

        # Run plots for this group subset
        cluster_dge_inspect_env$plot_comparison_scatterplots(dge_subset)
        cluster_dge_inspect_env$plot_stacked_barchart(dge_subset)
        cluster_dge_inspect_env$plot_volcano_plots(dge_subset)

        # Restore original paths
        cluster_dge_inspect_env$plot_data_dir <- original_data_dir
        cluster_dge_inspect_env$plot_pdf_dir <- original_pdf_dir
    }
}


# ==========================================
# Run all plots for the DGE on clusters, cell subtypes, and [optionally] combined cell types
# ==========================================

dge_types <- list(
    list(
        csv_path = paste0(output_dirs$deg_clusters, "merged_dge_on_clusters.csv"),
        plot_data_dir = paste0(output_dirs$deg_clusters, "inspection_plot_data/"),
        plot_pdf_dir = paste0(output_dirs$deg_clusters, "inspection_plots_pdf/")
    ),
    list(
        csv_path = paste0(output_dirs$deg_celltypes, "merged_dge_on_cellTypes.csv"),
        plot_data_dir = paste0(output_dirs$deg_celltypes, "inspection_plot_data/"),
        plot_pdf_dir = paste0(output_dirs$deg_celltypes, "inspection_plots_pdf/")
    )
)

# Optionally add combine_celltypes if present
if (!is.null(output_dirs$deg_combine_celltypes)) {
    dge_types[[length(dge_types) + 1]] <- list(
        csv_path = paste0(output_dirs$deg_combine_celltypes, "merged_dge_on_cellTypes.csv"),
        plot_data_dir = paste0(output_dirs$deg_combine_celltypes, "inspection_plot_data/"),
        plot_pdf_dir = paste0(output_dirs$deg_combine_celltypes, "inspection_plots_pdf/")
    )
}

for (dge_info in dge_types) {
    dge_csv_path <- dge_info$csv_path
    if (!file.exists(dge_csv_path)) {
        warning("Merged DGE CSV file not found. Please generate it before running this script: ", dge_csv_path)
        next
    }
    cat("Processing DGE file:", dge_csv_path, "\n")
    cluster_dge_inspect_env$plot_data_dir <- dge_info$plot_data_dir
    cluster_dge_inspect_env$plot_pdf_dir <- dge_info$plot_pdf_dir
    dir.create(cluster_dge_inspect_env$plot_data_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(cluster_dge_inspect_env$plot_pdf_dir, showWarnings = FALSE, recursive = TRUE)
    cluster_dge_inspect_env$run_all_plots(dge_csv_path)
}
