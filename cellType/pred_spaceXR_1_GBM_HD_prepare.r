# https://satijalab.org/seurat/articles/spatial_vignette#spatial-deconvolution-using-rctd
# https://satijalab.org/seurat/articles/seurat5_spatial_vignette_2#cb29

# devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)

source("util_headers.r")

merged_obj <- readRDS(paste0(project_dir, "rds_data/", output.file.prefix, "_clustered_12k.rds"))
merged_obj_list <- SplitObject(merged_obj, split.by = "batch")

# # =================================================
# #  scRNA-seq dataset reference 1
# # https://www.dropbox.com/s/cs6pii5my4p3ke3/mouse_hippocampus_reference.rds?dl=0
# ref <- readRDS(paste0(project_dir,"rds_data/",  "mouse_hippocampus_reference.rds"))
# ref <- Seurat::UpdateSeuratObject(ref)
# Idents(ref) <- "celltype"
# counts <- ref[["RNA"]]$counts
# cluster <- as.factor(ref$celltype)
# names(cluster) <- colnames(ref)
# nUMI <- ref$nCount_RNA
# names(nUMI) <- colnames(ref)
# reference <- spacexr::Reference(counts, cluster, nUMI)
# saveRDS(reference, paste0(project_dir,"rds_data/",  "ref_mouse_hippocampus_reference.rds"))
# reference <- readRDS(paste0(project_dir,"rds_data/",  "ref_mouse_hippocampus_reference.rds"))
# # =================================================


# # =================================================
# #  scRNA-seq dataset reference 2
# # allen.corted.ref can be downloaded here:
# # https://www.dropbox.com/s/cuowvm4vrf65pvq/allen_cortex.rds?dl=1
# allen.cortex.ref <- readRDS(paste0(project_dir,"rds_data/",  "allen_cortex.rds"))
# allen.cortex.ref <- Seurat::UpdateSeuratObject(allen.cortex.ref)
# Idents(allen.cortex.ref) <- "subclass"
# # remove CR cells because there aren't enough of them for annotation
# allen.cortex.ref <- subset(allen.cortex.ref, subset = subclass != "CR")
# counts <- allen.cortex.ref[['RNA']]$counts
# cluster <- as.factor(allen.cortex.ref$subclass)
# names(cluster) <- colnames(allen.cortex.ref)
# nUMI <- allen.cortex.ref$nCount_RNA
# names(nUMI) <- colnames(allen.cortex.ref)
# nUMI <- colSums(counts)
# levels(cluster) <- gsub("/", "-", levels(cluster))
# reference <- spacexr::Reference(counts, cluster, nUMI)
# saveRDS(reference, paste0(project_dir,"rds_data/",  "ref_allen_cortex.rds"))
# reference <- readRDS(paste0(project_dir,"rds_data/",  "ref_allen_cortex.rds"))
# # =================================================


# =================================================
#  scRNA-seq dataset reference 3
# https://doi.org/10.1038/s41590-022-01215-0
# sc_meta <- readRDS(paste0(old_project_dir, "meta_data.rds"))
# sc_count <- readRDS(paste0(old_project_dir, "raw_counts_sparse.rds"))
# cluster <- sc_meta$cell_type
# names(cluster) <- colnames(sc_count)
# nUMI <- colSums(sc_count)
# names(nUMI) <- colnames(sc_count)
# reference <- spacexr::Reference(sc_count, cluster, nUMI)
# saveRDS(reference, paste0(project_dir,"rds_data/",  "ref_GBM_scSeq.rds"))
reference <- readRDS(spaceXR_ref_rds)
# =================================================

for (batch_id in batch_names) {
    # slide.seq <- subset(
    #     merged_obj_list[[batch_id]],
    #     cells = Cells(merged_obj_list[[batch_id]][["Spatial.008um"]]),
    #     downsample = 15000)
    slide.seq <- subset(
        merged_obj_list[[batch_id]],
        cells = Cells(merged_obj_list[[batch_id]][[ifelse(VisiumHD, "sketch", "Spatial")]])
    )
    # set up query with the RCTD function SpatialRNA
    counts <- slide.seq[[ifelse(VisiumHD, "sketch", "Spatial")]]$counts
    coords <- GetTissueCoordinates(slide.seq)
    colnames(coords) <- c("x", "y")
    coords[is.na(colnames(coords))] <- NULL
    query <- spacexr::SpatialRNA(coords, counts, colSums(counts))

    RCTD <- spacexr::create.RCTD(query, reference, max_cores = 8)
    RCTD <- spacexr::run.RCTD(RCTD, doublet_mode = "doublet")
    slide.seq <- AddMetaData(slide.seq, metadata = RCTD@results$results_df)

    saveRDS(slide.seq, paste0(output_dirs$cell_type_deconv, "spaceXR_GBM_scSeq_", batch_id, "_sketch.rds"))
    # slide.seq <- readRDS(paste0(output_dirs$cell_type_deconv, "spaceXR_GBM_scSeq_", batch_id, "_sketch.rds"))

    p1 <- ImageDimPlot(slide.seq,
        group.by = "first_type",
        size = 1.5, cols = "polychrome",
        dark.background = F
    ) + ggplot2::ggtitle("Cell type - first_type")
    p2 <- ImageDimPlot(slide.seq,
        group.by = "second_type",
        size = 1.5, cols = "polychrome",
        dark.background = F
    ) + ggplot2::ggtitle("Cell type - second_type")

    # Save plots to cell_type_prediction/visualizations/
    ggpubr::ggexport(
        filename = paste0(output_dirs$cell_type_viz, "spaceXR_GBM_scSeq_", batch_id, "_predict_cellType_sketch.pdf"),
        plotlist = list(p1, p2)
    )
}


# ==================================
# map cell type to clusters and write to merged_obj
# ==================================

Idents(merged_obj) <- valid_cluster_method
projected_seurat_cluster <- Idents(merged_obj)

first_type_combined <- NULL
second_type_combined <- NULL

for (batch_id in batch_names) {
    slide.seq <- readRDS(paste0(output_dirs$cell_type_deconv, "spaceXR_GBM_scSeq_", batch_id, "_sketch.rds"))
    pred_first_type <- slide.seq$first_type
    names(pred_first_type) <- names(slide.seq$first_type)
    pred_second_type <- slide.seq$second_type
    names(pred_second_type) <- names(slide.seq$second_type)

    if (is.null(first_type_combined)) {
        first_type_combined <- pred_first_type[names(projected_seurat_cluster)]
        second_type_combined <- pred_second_type[names(projected_seurat_cluster)]
    } else {
        na_positions <- is.na(first_type_combined)
        first_type_combined[na_positions] <- pred_first_type[names(projected_seurat_cluster)[na_positions]]
        na_positions <- is.na(second_type_combined)
        second_type_combined[na_positions] <- pred_second_type[names(projected_seurat_cluster)[na_positions]]
    }
}

names(first_type_combined) <- names(projected_seurat_cluster)
names(second_type_combined) <- names(projected_seurat_cluster)

table(first_type_combined, useNA = "ifany") |> print()
table(second_type_combined, useNA = "ifany") |> print()

# Export cluster assignments to CSV
data.frame(Cell = rownames(merged_obj@meta.data), predict_GBM_scSeq_cellType_1 = first_type_combined) |>
    write.csv(file = paste0(output_dirs$clustering, "predict_GBM_scSeq_cellType_1.csv"), row.names = FALSE)

data.frame(Cell = rownames(merged_obj@meta.data), predict_GBM_scSeq_cellType_2 = second_type_combined) |>
    write.csv(file = paste0(output_dirs$clustering, "predict_GBM_scSeq_cellType_2.csv"), row.names = FALSE)


quit("no")
