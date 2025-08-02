# https://www.youtube.com/watch?v=y0xDk_P2V7A
# https://github.com/jaleesr/BITFAM
raw.data.dir <- "/app/data/"
output.dir <- paste0(raw.data.dir, "TFA/output/")
dir.create(output.dir, recursive = TRUE)
output_rds.dir <- paste0(raw.data.dir, "TFA/tmp_rds/")
dir.create(output_rds.dir, recursive = TRUE)


library(Seurat)
library(ggplot2)

config <- yaml::read_yaml(paste0(raw.data.dir, "config/batch_config.yaml"))
batch_names <- config$batch_names
output.file.prefix <- config$output_file_prefix
VisiumHD <- config$VisiumHD

merged_obj <- readRDS(paste0(raw.data.dir, "rds_data/", output.file.prefix, "_merged.rds"))
merged_obj_list <- SplitObject(merged_obj, split.by = "batch")

for (batch_id in batch_names) {
    one_subset <- merged_obj_list[[batch_id]]

    normalized_mtx <- LayerData(one_subset, assay = ifelse(VisiumHD, "sketch", "Spatial"), layer = "data")
    data_matrix_normalized <- BITFAM::BITFAM_preprocess(raw_data = normalized_mtx)
    # extract predict_GBM_scSeq_cellType_1 and map only those exist in `Z`
    predict_GBM_scSeq_cellType_1 <- one_subset$predict_GBM_scSeq_cellType_1

    # for mouse, the gene names should be only the first letter uppercase
    # head(rownames(data_matrix_normalized), n = 5)
    BITFAM_res <- BITFAM::BITFAM(data = data_matrix_normalized, species = "mouse", scATAC_obj = NA, ncores = parallel::detectCores())
    # attach seurat clusters in the model result
    BITFAM_res$cell_types <- predict_GBM_scSeq_cellType_1[BITFAM_res$Cell_names]
    saveRDS(BITFAM_res, paste0(output_rds.dir, "BITFAM_res_", batch_id, "_12k.rds"))
}


# plotting the TFA with UMAP
for (batch_id in batch_names) {
    BITFAM_res <- readRDS(paste0(output_rds.dir, "BITFAM_res_", batch_id, "_12k.rds"))
    infer_TF_activities <- BITFAM::BITFAM_activities(BITFAM_res)
    infer_TF_activities |>
        cbind(cell_types = BITFAM_res$cell_types) |>
        write.csv(paste0(output.dir, batch_id, "_TF_activities", ".csv"), row.names = TRUE)

    # plot a UMAP on the inferred transcription factor activities and assign seurat clusters to the spots
    umap_df <- uwot::umap(infer_TF_activities, n_neighbors = 15, min_dist = 0.1, metric = "euclidean") |> as.data.frame()
    colnames(umap_df) <- c("UMAP1", "UMAP2")
    umap_df$clusters <- BITFAM_res$cell_types

    p1 <- ggplot(umap_df, aes(UMAP1, UMAP2, color = clusters)) +
        geom_point(size = 0.8)
    ggsave(paste0(output.dir, batch_id, "_TF_activities", ".pdf"), plot = p1)
}


# export cell type factors and corresponding value items
BITFAM_res <- readRDS(paste0(output_rds.dir, "BITFAM_res_", batch_names[length(batch_names)], "_12k.rds"))
all_cell_types <- BITFAM_res$cell_types
for (batch_id in batch_names[-length(batch_names)]) {
    BITFAM_res <- readRDS(paste0(output_rds.dir, "BITFAM_res_", batch_id, "_12k.rds"))
    # combine factors. Factor vectors can't be combined with c()
    all_cell_types <- unlist(list(all_cell_types, BITFAM_res$cell_types))
}
all_cell_types_uni <- unique(all_cell_types)

write.csv(all_cell_types_uni, paste0(output.dir, "ordered_cell_types.csv"))


cat("========================\n")
cat("TFA done\n")
cat("========================\n")
quit("no")
