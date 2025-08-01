# spatialTrans_data_prefix <- "/vol/research/scratch1/NOBACKUP/hh01116/10X_Genomics_KD_data/"
spatialTrans_data_prefix <- "/app/data/"
output.dir <- "/app/output/"

###################################################################
# read the pre-deconvoluted datasets, plot
###################################################################
colors <- c(
  "#9e0142", "#d62828", "#d53e4f", "#f46d43", "#fca311", "#fdae61", "#fee08b",
  "#ffffbf", "#e5e5e5", "#eeb1d5", "#a9def9", "#abdda4", "#66c2a5", "#3288bd", "#5e4fa2",
  "#003049", "#14213d"
)
# barplot(rep(1, 17), col=colors, border=NA, space=0)

## interested cell types
ct.visualize <- c(
  # "Ciliated cells", "Oligodendrocyte precursor cells", "Neuroendocrine cells",
  # "Oligodendrocytes", "Endothelial cells", "Tuft cells", "Microglial cells",
  # "Astrocytes", "Macrophages", "NK cells", "Fibroblasts", "Dendritic cells",
  # "T cells", "Neutrophils", "B cells", "Mast cells", "Monocytes"
  "Microglial cells", "Macrophages", "NK cells", "Dendritic cells",
  "T cells", "Neutrophils", "B cells", "Mast cells", "Monocytes"
)

p_list <- list()
proportion_list <- list()

for (batch in c("SAL", "KD", "RAD", "KD_RAD")) {
  CARD_obj <- readRDS(paste0(spatialTrans_data_prefix, "cellType_deconvolute/deconvoluted/", batch, "_CARD_minGene500_minSpot10.Rds"))

  p_list[[length(p_list) + 1]] <- CARD::CARD.visualize.pie(
    proportion = CARD_obj@Proportion_CARD[, ct.visualize],
    spatial_location = CARD_obj@spatial_location,
    radius = 2.5,
    colors = colors
  ) + ggplot2::theme(
    legend.text = ggplot2::element_text(size = 8),
    legend.title = ggplot2::element_text(size = 8, face = "bold")
  ) + ggplot2::labs(title = paste(batch, "spatial distribution of the cell type proportion"))

  # ggplot2::ggsave(paste0(
  #     "/user/HS502/hh01116/jobs/YMa_CARD_docker/output/",
  #     batch, "_pie_chart.pdf"
  #     ), plot = p1)

  p_list[[length(p_list) + 1]] <- CARD::CARD.visualize.prop(
    proportion = CARD_obj@Proportion_CARD,
    spatial_location = CARD_obj@spatial_location,
    ct.visualize = ct.visualize,
    colors = c("lightblue", "lightyellow", "red"),
    NumCols = 5, ### number of columns in the figure panel
    pointSize = 3.0
  ) + ggplot2::theme(
    strip.text = ggplot2::element_text(size = 8), # adjust facet subtitle
  ) + ggplot2::labs(title = paste(batch, "spatial distribution on each cell type"))

  p_list[[length(p_list) + 1]] <- CARD::CARD.visualize.Cor(CARD_obj@Proportion_CARD[, ct.visualize], colors = NULL) +
    ggplot2::theme(
      legend.text = ggplot2::element_text(size = 8),
      legend.title = ggplot2::element_text(size = 8, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 10),
      axis.text.y = ggplot2::element_text(size = 10),
      plot.title = ggplot2::element_text(size = 12)
    ) +
    ggplot2::labs(title = paste(batch, "cell type correlation"))

    proportion_list[[batch]] <- CARD_obj@Proportion_CARD
}

ggpubr::ggexport(p_list, filename = paste0(
  output.dir,
  "GBM_HD_deconvolution_analysis.pdf"
))

saveRDS(proportion_list, paste0(spatialTrans_data_prefix, "cellType_deconvolute/deconvoluted/CARD_cellType_prop_fourBatches.Rds"))


# two cell types comparison
# p3 <- CARD::CARD.visualize.prop.2CT(
#     proportion = CARD_obj@Proportion_CARD,
#     spatial_location = CARD_obj@spatial_location,
#     ct2.visualize = c("T cells", "Macrophages"),
#     colors = list(c("lightblue", "lightyellow", "red"),
#                   c("lightblue", "lightyellow", "black"))
# )
# print(p3)

