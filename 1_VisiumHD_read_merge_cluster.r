# https://satijalab.org/seurat/articles/seurat5_integration
# https://satijalab.org/seurat/articles/essential_commands.html#merge-objects-without-integration
# https://satijalab.org/seurat/articles/visiumhd_analysis_vignette


# ===========================================
#           Setup and Configuration
# ===========================================
source("util_headers.r")

lapply(output_dirs, function(dir) dir.create(dir, recursive = TRUE, showWarnings = FALSE))

options(future.globals.maxSize = 45 * 1024^3) # 45GB

dir.create(paste0(project_dir, "rds_data/"), recursive = TRUE)


# ===========================================
#           Read from 10X_spatial folder to rds
# ===========================================
for (i in seq_along(batch_file_names)) {
    prefix <- batch_names[i]
    if (VisiumHD) {
        gbm_subset <- Load10X_Spatial(
            data.dir = paste0(raw_data_dir, batch_file_names[i], "/outs"),
            bin.size = c(8, 16),
            slice = paste0(prefix, "_slice")
        )
    } else {
        gbm_subset <- Seurat::Load10X_Spatial(
            data.dir = paste0(raw_data_dir, batch_file_names[i]),
            filename = batch_file_h5_names[i],
            assay = "Spatial",
            slice = paste0(prefix, "_slice"),
            filter.matrix = TRUE,
            to.upper = FALSE
        )
    }
    gbm_subset[["batch"]] <- prefix
    saveRDS(gbm_subset, paste0(project_dir, "rds_data/", prefix, "_raw.rds"))
}

# ===========================================
#           Sketch sampling to reduce size
# ===========================================
brain_data_list <- list()

for (i in seq_along(batch_file_names)) {
    prefix <- batch_names[i]
    gbm_subset <- readRDS(paste0(project_dir, "rds_data/", prefix, "_raw.rds"))
    colnames(gbm_subset) <- paste0(prefix, "_", colnames(gbm_subset))
    gbm_subset@meta.data$orig.ident <- gbm_subset@meta.data$batch

    if (VisiumHD) {
        Seurat::DefaultAssay(gbm_subset) <- "Spatial.008um"
    }
    gbm_subset <- Seurat::NormalizeData(gbm_subset) |> Seurat::FindVariableFeatures()

    if (VisiumHD) {
        gbm_subset <- Seurat::SketchData(
            object = gbm_subset,
            ncells = 12000,
            method = "LeverageScore",
            sketched.assay = "sketch"
        )

        # switch analysis to sketched cells
        Seurat::DefaultAssay(gbm_subset) <- "sketch"
    }

    # perform clustering workflow (it's normalised)
    # gbm_subset <- Seurat::FindVariableFeatures(gbm_subset) |>
    #   Seurat::ScaleData() |>
    #   Seurat::RunPCA(assay = "sketch",reduction.name = "sketch_PCA") |>
    #   Seurat::FindNeighbors(assay = "sketch", reduction = "sketch_PCA", dims = 1:50) |>
    #   Seurat::FindClusters(cluster.name = "sketch_seurat_cluster", resolution = 3) |>
    #   Seurat::RunUMAP(reduction = "sketch_PCA", reduction.name = "sketch_umap", return.model = T, dims = 1:50)


    brain_data_list[[i]] <- gbm_subset
}


# merge
merged_obj <- merge(
    x = brain_data_list[[1]],
    y = brain_data_list[-1]
)

# ===========================================
#           PCA & Integrate->Join & UMAP
# ===========================================
if (VisiumHD) {
    Seurat::DefaultAssay(merged_obj) <- "sketch"
}
merged_obj <- NormalizeData(merged_obj) |>
    FindVariableFeatures() |>
    ScaleData() |>
    RunPCA(assay = ifelse(VisiumHD, "sketch", "Spatial"), reduction.name = ifelse(VisiumHD, "sketch_PCA", "PCA"))

merged_obj <- IntegrateLayers(
    merged_obj,
    method = RPCAIntegration,
    orig.reduction = ifelse(VisiumHD, "sketch_PCA", "PCA"),
    new.reduction = ifelse(VisiumHD, "sketch_integrated.rpca", "integrated.rpca"),
    verbose = FALSE
)

# ===========================================
#           Join Layers
# ===========================================
if (VisiumHD) {
    merged_obj <- JoinLayers(merged_obj, assay = "sketch")
    merged_obj <- JoinLayers(merged_obj, assay = "Spatial.008um")
    merged_obj <- JoinLayers(merged_obj, assay = "Spatial.016um")
}

# ===========================================
#           Clustering
# ===========================================
if (VisiumHD) {
    DefaultAssay(merged_obj) <- "sketch"
}
merged_obj <- FindNeighbors(
    merged_obj,
    assay = ifelse(VisiumHD, "sketch", "Spatial"),
    reduction = ifelse(VisiumHD, "sketch_integrated.rpca", "integrated.rpca"),
    dims = 1:30,
    k.param = Cluster_n_neighbors
) |>
    FindClusters(cluster.name = ifelse(VisiumHD, "sketch_seurat_cluster", "seurat_cluster"), resolution = Cluster_resolution) |>
    RunUMAP(
        reduction.name = ifelse(VisiumHD, "sketch_umap", "umap"),
        return.model = T, dims = 1:30,
        reduction = ifelse(VisiumHD, "sketch_integrated.rpca", "integrated.rpca"), n.components = 3
    )

saveRDS(merged_obj, paste0(project_dir, "rds_data/", "tmp.rds"))
# Remove all variables except some
rm(list = setdiff(ls(), c("merged_obj", "VisiumHD", "project_dir", "output.file.prefix")))
# Run garbage collection to free up memory
gc()

cat("===========================================\n")
cat("Done the first preprocessing and clustering\n")
cat("Now proceed to projecting data which requires huge memory\n")
cat("If it's killed, please manually run from the below commands.\n")
cat("===========================================\n")
# merged_obj <- readRDS(paste0(project_dir, "rds_data/", "tmp.rds"))

# ===========================================
#           Project Data
# ===========================================
if (VisiumHD) {
    merged_obj <- Seurat::ProjectData(
        object = merged_obj,
        assay = "Spatial.008um",
        full.reduction = "full.sketch_integrated.rpca",
        sketched.assay = "sketch",
        sketched.reduction = "sketch_integrated.rpca",
        umap.model = "sketch_umap",
        dims = 1:50,
        refdata = list(projected_seurat_cluster = "sketch_seurat_cluster")
    )
}

# cell_region_cluster <- merged_obj$projected_seurat_cluster
# merged_obj$cell_region_cluster <- dplyr::case_when(
#       cell_region_cluster %in% c("1", "2", "5", "7", "9") ~ "outTumour",
#       cell_region_cluster %in% c("3", "18", "11", "14") ~ "inTumour",
#       cell_region_cluster %in% c("4", "13") ~ "edgeTumour",
#       .default = as.character(cell_region_cluster)
#     ) |> as.factor()

saveRDS(merged_obj, paste0(project_dir, "rds_data/", output.file.prefix, "_clustered_12k.rds"))
# str(merged_obj_list[['SAL']]) |> capture.output() |>writeLines(con=paste0(project_dir,"rds_data/",  "two_merged_clustered_12k_SAL.txt"))

file.remove(paste0(project_dir, "rds_data/", "tmp.rds"))

cat("===========================================\n")
cat("Done preprocessing and clustering\n")
cat("===========================================\n")
quit("no")
