# ===========================================
#           Setup and Configuration
# ===========================================
message("=== Setup and Configuration ===")

source("util_headers.r") # Loads config and output_dirs

# Create output directories
lapply(output_dirs, function(dir) dir.create(dir, recursive = TRUE, showWarnings = FALSE))
dir.create(rds_data_dir, recursive = TRUE)

# Set memory options
options(future.globals.maxSize = on_disk$max_memory_gb * 1024^3)

# ===========================================
#           Read and Preprocess Batches (On-Disk)
# ===========================================
message("=== Reading and Preprocessing Batches ===")
target_assay <- ifelse(VisiumHD, "sketch", "Spatial")
assay_name <- ifelse(VisiumHD, "Spatial.008um", "Spatial")

low_count_threshold <- 50000

batch_counts_df <- data.frame(
    batch = character(),
    total_counts = numeric(),
    n_cells = numeric(),
    low_counts = logical(),
    stringsAsFactors = FALSE
)

filter.genes.cells <- function(obj, min.value, min.cells, min.genes) {
    data.slot <- Seurat::GetAssayData(obj, layer = "counts")
    num.cells <- Matrix::rowSums(data.slot >= min.value)
    genes.use <- names(num.cells[which(num.cells >= min.cells)])
    obj <- subset(obj, features = genes.use)
    num.genes <- Matrix::colSums(data.slot >= min.value)
    obj$num.genes <- num.genes
    obj <- subset(obj, subset = num.genes >= min.genes)
    obj$num.genes <- NULL
    return(obj)
}

# --- Process all batches ---
for (i in seq_along(batch_file_names)) {
    prefix <- batch_names[i]
    message("Processing batch: ", prefix)
    rds_path <- paste0(rds_data_dir, prefix, "_raw.rds")

    if (file.exists(rds_path)) {
        gbm_subset <- readRDS(rds_path)
    } else {
        if (VisiumHD) {
            gbm_subset <- Seurat::Load10X_Spatial(
                data.dir = paste0(raw_data_dir, batch_file_names[i], "/outs"),
                bin.size = c(8),
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
        colnames(gbm_subset) <- paste0(prefix, "_", colnames(gbm_subset))
        gbm_subset@meta.data$orig.ident <- gbm_subset@meta.data$batch
        gbm_subset <- filter.genes.cells(gbm_subset, min.value = 0, min.cells = on_disk$qc_min_cells, min.genes = on_disk$qc_min_genes)
        # https://satijalab.org/seurat/articles/seurat5_bpcells_interaction_vignette
        # --- Convert counts layer to on-disk matrix ---

        counts_dir <- paste0(rds_data_dir, prefix, "_counts_bpcells")
        counts_matrix <- gbm_subset[[assay_name]]$counts
        if (!is.integer(counts_matrix[1, 1])) {
            counts_matrix <- BPCells::convert_matrix_type(counts_matrix, "uint32_t")
        }
        BPCells::write_matrix_dir(mat = counts_matrix, dir = counts_dir)
        gbm_subset[[assay_name]]$counts <- BPCells::open_matrix_dir(dir = counts_dir)
        saveRDS(gbm_subset, rds_path)
    }

    # --- Inspect counts ---
    total_counts <- sum(Matrix::colSums(Seurat::GetAssayData(gbm_subset, layer = "counts")))
    n_cells <- ncol(gbm_subset)
    is_low_count <- total_counts < low_count_threshold
    batch_counts_df <- rbind(batch_counts_df, data.frame(
        batch = prefix,
        total_counts = total_counts,
        n_cells = n_cells,
        low_counts = is_low_count
    ))
    gc()
}

p_counts <- ggplot(batch_counts_df, aes(x = reorder(batch, total_counts), y = total_counts)) +
    geom_col(aes(fill = low_counts), alpha = 0.8) +
    geom_hline(yintercept = low_count_threshold, linetype = "dashed", color = "red", alpha = 0.7) +
    scale_fill_manual(
        values = c("FALSE" = "steelblue", "TRUE" = "orange"),
        labels = c("Normal", "Low counts"), name = "Count Status"
    ) +
    scale_y_continuous(labels = scales::comma_format()) +
    labs(
        title = "Total Counts per Batch",
        subtitle = paste("Red dashed line indicates low count threshold (", scales::comma(low_count_threshold), ")"),
        x = "Batch", y = "Total Counts"
    ) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)
    ) +
    geom_text(aes(label = paste0(scales::comma(total_counts), "\n(", scales::comma(n_cells), " cells)")),
        hjust = 0.5, size = 3, angle = 0
    )
ggsave(paste0(output_dirs$clustering, "batch_total_counts.pdf"), p_counts, width = 12, height = 8)
message("Saved total counts visualization to: ", paste0(output_dirs$clustering, "batch_total_counts.pdf"))


# ===========================================
#           Mitochondrial QC Analysis
# ===========================================
message("=== Mitochondrial QC Analysis ===")
plot_mitochondrial_qc <- function(batch_names, rds_data_dir, output_dirs) {
    mito_qc_list <- list()
    for (batch_name in batch_names) {
        gbm_subset <- readRDS(paste0(rds_data_dir, batch_name, "_raw.rds"))
        mito_genes <- grep("^[Mm][Tt]-", rownames(gbm_subset), value = TRUE)
        if (length(mito_genes) > 0) {
            gbm_subset <- PercentageFeatureSet(gbm_subset, pattern = "^[Mm][Tt]-|^MT-", col.name = "percent.mt", assay = assay_name)
            mito_qc_list[[batch_name]] <- data.frame(
                batch = batch_name,
                cell_id = colnames(gbm_subset),
                percent_mt = gbm_subset$percent.mt,
                n_genes = ifelse(VisiumHD, gbm_subset$nFeature_Spatial.008um, gbm_subset$nFeature_Spatial),
                n_counts = ifelse(VisiumHD, gbm_subset$nCount_Spatial.008um, gbm_subset$nCount_Spatial)
            )
        } else {
            mito_qc_list[[batch_name]] <- data.frame(
                batch = batch_name,
                cell_id = character(0),
                percent_mt = numeric(0),
                n_genes = numeric(0),
                n_counts = numeric(0)
            )
        }
    }
    mito_qc_df <- do.call(rbind, mito_qc_list)
    if (nrow(mito_qc_df) > 0) {
        p_mito <- ggplot(mito_qc_df, aes(x = batch, y = percent_mt, fill = batch)) +
            geom_violin(alpha = 0.7) +
            geom_boxplot(width = 0.1, fill = "white", alpha = 0.8) +
            labs(
                title = "Mitochondrial Gene Expression Percentage by Batch",
                subtitle = "Use this plot to determine mitochondrial cutoff threshold",
                x = "Batch", y = "Mitochondrial Gene %"
            ) +
            theme_minimal() +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1),
                plot.title = element_text(hjust = 0.5),
                plot.subtitle = element_text(hjust = 0.5),
                legend.position = "none"
            ) +
            geom_hline(yintercept = c(10, 20, 30), linetype = "dashed", alpha = 0.5, color = "red")
        ggsave(paste0(output_dirs$clustering, "mitochondrial_qc_by_batch.pdf"), p_mito, width = 12, height = 8)
        return(mito_qc_df)
    } else {
        return(NULL)
    }
}
mito_qc_data <- plot_mitochondrial_qc(batch_names, rds_data_dir, output_dirs)

integrated_merged_obj_path <- paste0(rds_data_dir, output.file.prefix, "_merged.rds")
if (file.exists(integrated_merged_obj_path)) {
    message("\nMerged object already exists, load the object: ", integrated_merged_obj_path)
    merged_obj <- readRDS(integrated_merged_obj_path)
    integration_success <- TRUE
} else {
    message("\nMerged object does not exist, proceeding with sketch sampling, integration, and merging.")

    # ===========================================
    #   NormalizeData, FindVariableFeatures
    #   and    Sketch Sampling (On-Disk)
    # ===========================================
    # https://satijalab.org/seurat/articles/seurat5_sketch_analysis

    message("=== Sketch Sampling ===")
    sketched_paths <- list()
    for (i in seq_along(batch_file_names)) {
        prefix <- batch_names[i]
        sketched_rds_path <- paste0(rds_data_dir, prefix, ifelse(VisiumHD, "_sketched.rds", "_raw.rds"))
        if (file.exists(sketched_rds_path)) {
            sketched_paths[[prefix]] <- sketched_rds_path
            next
        }
        gbm_subset <- readRDS(paste0(rds_data_dir, prefix, "_raw.rds"))
        if (VisiumHD) {
            Seurat::DefaultAssay(gbm_subset) <- "Spatial.008um"
            gbm_subset <- Seurat::NormalizeData(gbm_subset) |> Seurat::FindVariableFeatures()
            gbm_subset <- Seurat::SketchData(
                object = gbm_subset,
                ncells = 12000,
                method = "Uniform",
                sketched.assay = "sketch"
            )
            Seurat::DefaultAssay(gbm_subset) <- "sketch"
            saveRDS(gbm_subset, sketched_rds_path)
            sketched_paths[[prefix]] <- sketched_rds_path
        } else {
            gbm_subset <- Seurat::NormalizeData(gbm_subset) |> Seurat::FindVariableFeatures()
            saveRDS(gbm_subset, sketched_rds_path)
            sketched_paths[[prefix]] <- sketched_rds_path
        }
        gc()
    }

    # ===========================================
    #           Merge (Sketched) Objects (On-Disk)
    # ===========================================
    message("=== Merging (Sketched) Objects (On-Disk) ===")
    merged_obj <- NULL
    for (i in seq_along(batch_names)) {
        obj <- readRDS(sketched_paths[[batch_names[i]]])
        obj@tools <- list()
        obj@misc <- list()
        if (is.null(merged_obj)) {
            merged_obj <- obj
        } else {
            merged_obj <- merge(merged_obj, y = obj)
        }
        rm(obj)
        gc()
    }

    # ===========================================
    #           Integration
    # ===========================================
    message("=== Integration (On-Disk) ===")
    merged_obj <- NormalizeData(merged_obj) |>
        FindVariableFeatures() |>
        ScaleData() |>
        RunPCA(reduction.name = "PCA")

    integration_success <- FALSE
    # Validate assays before integration
    if (target_assay %in% names(merged_obj@assays)) {
        message("Target assay '", target_assay, "' found. Attempting integration...")

        tryCatch(
            {
                DefaultAssay(merged_obj) <- target_assay

                # Ensure the object has the required structure for integration
                if (length(unique(merged_obj$batch)) > 1) {
                    merged_obj <- IntegrateLayers(
                        merged_obj,
                        method = RPCAIntegration,
                        orig.reduction = "PCA",
                        new.reduction = ifelse(VisiumHD, "sketch_integrated.rpca", "integrated.rpca"),
                        verbose = TRUE
                    )
                    integration_success <- TRUE
                    message("Integration completed successfully")
                } else {
                    message("Only one batch detected, skipping integration")
                    integration_success <- FALSE
                }
            },
            error = function(e) {
                message("Error during integration: ", e$message)
                message("Proceeding without integration...")
                integration_success <- FALSE
            }
        )
    } else {
        message("Target assay '", target_assay, "' not found in merged object.")
        message("Available assays: ", paste(names(merged_obj@assays), collapse = ", "))
        integration_success <- FALSE
    }


    # ===========================================
    #           Join Layers
    # ===========================================
    merged_obj <- JoinLayers(merged_obj)

    # --- Convert merged counts layer to on-disk matrix ---
    merged_counts_dir <- paste0(rds_data_dir, output.file.prefix, "_merged_counts_bpcells")
    merged_counts_matrix <- Seurat::GetAssayData(merged_obj, assay = target_assay, layer = "counts")
    BPCells::write_matrix_dir(mat = merged_counts_matrix, dir = merged_counts_dir)
    merged_obj[[target_assay]]$counts <- BPCells::open_matrix_dir(dir = merged_counts_dir)
    merged_data_dir <- paste0(rds_data_dir, output.file.prefix, "_merged_data_bpcells")
    merged_data_matrix <- Seurat::GetAssayData(merged_obj, assay = target_assay, layer = "data")
    BPCells::write_matrix_dir(mat = merged_data_matrix, dir = merged_data_dir)
    merged_obj[[target_assay]]$data <- BPCells::open_matrix_dir(dir = merged_data_dir)
    saveRDS(merged_obj, paste0(rds_data_dir, output.file.prefix, "_merged.rds"))
}

# ===========================================
#           Clustering
# ===========================================
message("=== Clustering ===")

if (integration_success && ifelse(VisiumHD, "sketch_integrated.rpca", "integrated.rpca") %in% names(merged_obj@reductions)) {
    target_reduction <- ifelse(VisiumHD, "sketch_integrated.rpca", "integrated.rpca")
} else {
    target_reduction <- "PCA"
}

merged_obj <- FindNeighbors(
    merged_obj,
    assay = target_assay,
    reduction = target_reduction,
    dims = 1:30,
    k.param = Cluster_n_neighbors
) |>
    FindClusters(cluster.name = "seurat_clusters", resolution = Cluster_resolution)

merged_obj <- RunUMAP(
    merged_obj,
    reduction = target_reduction,
    reduction.name = "umap",
    return.model = TRUE,
    dims = 1:30,
    n.components = 3
)

# integration_results <- list(
#     reductions = merged_obj@reductions,
#     clusters = merged_obj$seurat_clusters,
#     cell_names = colnames(merged_obj)
# )
# saveRDS(integration_results, paste0(rds_data_dir, output.file.prefix, "_umap_clusters.rds"))
# DimPlot(merged_obj, label = T, label.size = 3, reduction = "umap", group.by = "seurat_clusters", alpha = 0.1) + NoLegend()

message("=== Projecting Data (On-Disk) ===")
merged_obj <- ProjectData(
    object = merged_obj,
    assay = assay_name,
    full.reduction = "pca.full",
    sketched.assay = "sketch",
    sketched.reduction = ifelse(VisiumHD, "sketch_integrated.rpca", "integrated.rpca"),
    umap.model = "umap",
    dims = 1:30,
    refdata = list(seurat_cluster_full = "seurat_clusters")
)
# DimPlot(merged_obj, label = T, label.size = 3, reduction = "full.umap", group.by = "seurat_cluster_full", alpha = 0.1) + NoLegend()

saveRDS(merged_obj, paste0(rds_data_dir, output.file.prefix, "_merged_clustered.rds"))

message("=== Cleaning Up ===")
rm(merged_obj)
gc()
message("===========================================")
message("Done preprocessing and clustering (on-disk).")
