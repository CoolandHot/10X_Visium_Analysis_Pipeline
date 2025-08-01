# BiocManager::install("ComplexHeatmap")
# install.packages("dendextend")
# install.packages("pheatmap")
# install.packages("magick")

setwd("./projects")
sapply(c("ggplot2", "Seurat", "patchwork", "RColorBrewer"), require, character.only = TRUE) |>
    suppressPackageStartupMessages()
# maintain same color scale in the SpatialFeaturePlot & SpatialDimPlot
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))


# https://www.datanovia.com/en/lessons/heatmap-in-r-static-and-interactive-visualization/
# https://www.youtube.com/watch?v=ht1r34-ifVI

diff_genes_fc <- read.csv("./composite_markers/de_output/DE_allTumorPOS_allGroups_unfiltered.csv") |>
    dplyr::filter(p_val_adj < 0.05) |>
    dplyr::select(gene, avg_log2FC, position, group)

diff_genes_fc_vs_sal <- diff_genes_fc |>
    dplyr::filter(grepl("_vs_SAL", group)) |>
    tidyr::pivot_wider(
        id_cols = gene,
        names_from = c(position, group),
        values_from = avg_log2FC
    ) |>
    dplyr::relocate(gene, dplyr::starts_with("inTumour_"), dplyr::starts_with("edgeTumour_"))


# -------------------
# + avg_log2FC +
# -------------------
# plot avg_log2FC of ADI, RAD, ADI+RAD vs SAL
# across inTumour, edgeTumour, outTumour

# convert to matrix with
# gene as rownames, columns: inTumour_ADI_RAD_vs_SAL, inTumour_ADI_vs_SAL, ..., outTumour_ADI_vs_SAL
diff_genes_fc_vs_sal_mtx <- diff_genes_fc_vs_sal |>
    # select genes: exists in ALL DE groups
    tidyr::drop_na() |>
    tibble::column_to_rownames(var = "gene") |>
    as.matrix()

# unify color scale
col_logFC <- circlize::colorRamp2(
    c(min(diff_genes_fc_vs_sal_mtx, na.rm = T), 0, max(diff_genes_fc_vs_sal_mtx, na.rm = T)),
    c("#FFCC70", "#C850C0", "#4158D0")
)

# start plots
pdf("./composite_markers/output/avg_log2FC_vs_SAL_heatmap_across_allGroups.pdf",
    width = 12, height = 15
)
ComplexHeatmap::Heatmap(
    diff_genes_fc_vs_sal_mtx,
    # color hierarchical clustering lines
    cluster_rows = dendextend::color_branches(hclust(dist(diff_genes_fc_vs_sal_mtx)), k = 2),
    row_title = "genes\ngene selected criteria: exists in ALL DE groups",
    cluster_columns = F,
    column_title = "average log2FC against SAL",
    column_names_rot = 45, # label rotate 45
    name = "log2FC",
    col = col_logFC,
    layer_fun = function(j, i, x, y, w, h, fill) { # add text to each grid using layer_fun
        grid::grid.text(
            sprintf(
                "%.3f",
                ComplexHeatmap::pindex(diff_genes_fc_vs_sal_mtx, i, j)
            ),
            x, y
        )
    }
)
dev.off()



# -------------------
# + gene expression +
# -------------------
gbm_combined_compr <- readRDS(file = "./composite_markers/gbm_combine_.rds")


#' wrap up function for plotting heatmap of the selected genes
#' with normalized gene expression on the left, fold change on the right
#'
#' @selected_genes
#' @geneExpr_y_title title for y.axis of the left heatmap
#' @geneExpr_row_cluster bool, cluster the rows?
#' @fc_mtx the matrix to be plotted on the right
#'          rownames are the selected genes
#'          columns are inTumour_ADI_vs_SAL, inTumour_RAD_vs_SAL, inTumour_ADI+RAD_vs_SAL
heatmap_geneExpr_FC <- function(selected_genes, geneExpr_y_title, geneExpr_row_cluster, fc_mtx, add_gene_count = FALSE) {
    tumor_pos_str <- gbm_combined_compr$tumor_pos_num
    tumor_pos_list <- dplyr::case_when(
        tumor_pos_str == "0" ~ "inTumour",
        tumor_pos_str == "1" ~ "edgeTumour",
        tumor_pos_str == "2" ~ "outTumour",
        .default = as.character(tumor_pos_str)
    )

    # extract the normalized gene expression from the integrated Seurat object
    # matrix for heatmap on the left
    # rownames is gene, columns are inTumour_ADI, edgeTumour_ADI, ..., outTumour_SAL
    gbm_mean_expr_mtx <- Seurat::FetchData(gbm_combined_compr,
        vars = selected_genes
    ) |>
        dplyr::mutate(
            group = gbm_combined_compr$batch,
            position = tumor_pos_list
        ) |>
        # compute the mean across `group, position`
        dplyr::group_by(group, position) |>
        dplyr::summarise(dplyr::across(dplyr::everything(), mean),
            .groups = "drop"
        ) |>
        tidyr::unite(pos_group, position, group, sep = "_") |>
        dplyr::mutate(order_column = dplyr::case_when(
            startsWith(pos_group, "inTumour_") ~ 1,
            startsWith(pos_group, "edgeTumour_") ~ 2,
            TRUE ~ 3 # If there are other rows not matching the conditions, place them last
        )) |>
        dplyr::arrange(order_column) |>
        dplyr::select(-order_column) |>
        tibble::column_to_rownames(var = "pos_group") |>
        as.matrix() |>
        t()

    h1 <- ComplexHeatmap::Heatmap(
        gbm_mean_expr_mtx,
        cluster_rows = geneExpr_row_cluster,
        row_title = geneExpr_y_title,
        cluster_columns = F,
        column_names_rot = 45, # label rotate 45
        column_title = "Mean gene expression",
        name = "Gene expression",
        col = circlize::colorRamp2(
            c(quantile(gbm_mean_expr_mtx)[1], quantile(gbm_mean_expr_mtx)[4]),
            c("#C850C0", "#4158D0")
        ),
        use_raster = TRUE,
        layer_fun = ifelse(add_gene_count, function(j, i, x, y, w, h, fill) { # add text to each grid using layer_fun
            grid::grid.text(
                sprintf(
                    "%.3f",
                    ComplexHeatmap::pindex(gbm_mean_expr_mtx, i, j)
                ),
                x, y,
                gp = grid::gpar(fontsize = 8)
            )
        },
        NULL
        )
    )

    h2 <- ComplexHeatmap::Heatmap(
        fc_mtx,
        row_labels = row.names(fc_mtx),
        cluster_rows = F,
        cluster_columns = F,
        col = circlize::colorRamp2(
            c(min(fc_mtx, na.rm = T), 0, max(fc_mtx, na.rm = T)),
            c("#FFCC70", "#C850C0", "#4158D0")
        ),
        name = "log2FC",
        column_title = "average log2FC",
        column_names_rot = 45,
        layer_fun = function(j, i, x, y, w, h, fill) { # add text to each grid using layer_fun
            grid::grid.text(
                sprintf(
                    "%.3f",
                    ComplexHeatmap::pindex(fc_mtx, i, j)
                ),
                x, y,
                gp = grid::gpar(fontsize = 8)
            )
        }
    )

    return(h1 + h2)
}


# extract normalized gene expression of the given genes
# the genes where rowSums(is.na(diff_genes_fc_vs_sal)) <= 3
diff_genes_fc_vs_sal_filtered <- diff_genes_fc_vs_sal[rowSums(is.na(diff_genes_fc_vs_sal)) <= 3, ] |>
    dplyr::filter(!grepl("mt\\.", gene) & !grepl("^Rp[sl]", gene))

# matrix for heatmap on the right
gbm_FC_values <- diff_genes_fc_vs_sal_filtered |>
    dplyr::select(gene, matches("inTumour_\\w+_vs_SAL")) |>
    tibble::column_to_rownames(var = "gene") |>
    as.matrix()

# start plots
pdf("./composite_markers/output/gene_express_heatmap_across_allGroups.pdf",
    width = 12, height = 30
)
p1 <- heatmap_geneExpr_FC(
    selected_genes = diff_genes_fc_vs_sal_filtered$gene,
    geneExpr_y_title = "genes\ngene selected criteria: exists in at least 75% DE groups & not mito-/ribo-gene",
    geneExpr_row_cluster = TRUE,
    fc_mtx = gbm_FC_values
)
p1
dev.off()


# only select the genes of interested immune cells
selected_genes <- c("Stat1", "Irf7", "Irf8", "Cx3cr1", "Ccl7", "Ccl2", "Ccl5", "Cxcl10")

selected_genes <- c("KLF2", "KLF6", "IL10", "VEGFA", "CCL4", "ARG2", "CXCL2", "CXCL8", "MARCO ", "HIF1A", "MAFB") |>
    sapply(
        function(gene) {
            result <- grep(gene, rownames(gbm_combined_compr), ignore.case = TRUE, value = TRUE)
            if (length(result) != 0) {
                return(result)
            }
        },
        simplify = FALSE
    ) |>
    unlist()


gbm_FC_values <- diff_genes_fc_vs_sal |>
    dplyr::select(gene, dplyr::matches("inTumour_")) |>
    dplyr::filter(gene %in% selected_genes) |>
    # add the mismatched genes in selected_genes
    {
        \(filter_subset) {
            missing_genes <- setdiff(selected_genes, filter_subset$gene)
            missing_rows <- matrix(NA,
                nrow = length(missing_genes),
                ncol = ncol(filter_subset)
            ) |> as.data.frame()
            colnames(missing_rows) <- colnames(filter_subset)
            missing_rows$gene <- missing_genes

            # combine the mismatched gene dataframe with the matched dataframe
            rbind(missing_rows, filter_subset)
        }
    }() |>
    dplyr::arrange(match(gene, selected_genes)) |> # make sure heatmap left and right are in the same order
    tibble::column_to_rownames(var = "gene") |>
    as.matrix()


pdf("./composite_markers/output/gene_express_heatmap_interested_immune_celss.pdf",
    width = 12, height = 8
)
p1 <- heatmap_geneExpr_FC(
    selected_genes = selected_genes,
    geneExpr_y_title = "genes",
    geneExpr_row_cluster = FALSE,
    fc_mtx = gbm_FC_values,
    add_gene_count = TRUE
)
p1
dev.off()

# -------------------
# +    pathways     +
# -------------------
