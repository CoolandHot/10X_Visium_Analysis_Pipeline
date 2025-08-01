# ===========================================
#           Setup and Configuration
# ===========================================
message("=== Setup and Configuration ===")

source("util_headers.r") # Assuming contents are integrated or loaded

# Create output directories
lapply(output_dirs, function(dir) dir.create(dir, recursive = TRUE, showWarnings = FALSE))
dir.create(rds_data_dir, recursive = TRUE)

# Set memory options
options(future.globals.maxSize = 45 * 1024^3) # 45GB

# ===========================================
#           Read and Preprocess Batches (On-Disk)
# ===========================================
message("=== Reading and Preprocessing Batches ===")

low_count_threshold <- 50000
batch_counts_df <- data.frame(
    batch = character(),
    total_counts = numeric(),
    n_cells = numeric(),
    low_counts = logical(),
    stringsAsFactors = FALSE
)

# --- Function for QC filtering ---
filter.genes.cells <- function(obj, min.value, min.cells, min.genes) {
    # filter genes
    data.slot <- Seurat::GetAssayData(obj, layer = "counts")
    # genes_per_cell <- Matrix::colSums(data.slot>0)
    num.cells <- Matrix::rowSums(data.slot >= min.value)
    genes.use <- names(num.cells[which(num.cells >= min.cells)])
    obj <- subset(obj, features = genes.use)

    # filter cells
    num.genes <- Matrix::colSums(data.slot >= min.value)
    obj$num.genes <- num.genes
    obj <- subset(obj, subset = num.genes >= min.genes)
    obj$num.genes <- NULL
    return(obj)
}

# --- Function for batch processing ---
process_batch <- function(i, batch_file_names, batch_names, raw_data_dir, rds_data_dir, low_count_threshold, VisiumHD) {
    prefix <- batch_names[i]
    message("Processing batch: ", prefix)

    # Check if raw RDS file already exists
    rds_path <- paste0(rds_data_dir, prefix, "_raw.rds")

    if (file.exists(rds_path)) {
        message("Raw RDS file already exists for batch ", prefix, ". Loading existing file...")
        gbm_subset <- readRDS(rds_path)
    } else {
        # 1. Load Data
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

        # 2. Add metadata
        gbm_subset[["batch"]] <- prefix
        colnames(gbm_subset) <- paste0(prefix, "_", colnames(gbm_subset)) # Prefix cell names
        gbm_subset@meta.data$orig.ident <- gbm_subset@meta.data$batch

        # 3. Apply QC filtering
        message("Applying QC filtering: min_genes=", on_disk$qc_min_genes, ", min_cells=", on_disk$qc_min_cells)
        gbm_subset <- filter.genes.cells(gbm_subset, min.value = 0, min.cells = on_disk$qc_min_cells, min.genes = on_disk$qc_min_genes)

        # 4. Save raw object (RDS)
        saveRDS(gbm_subset, rds_path)
    }

    # 3. Inspect counts
    if (VisiumHD) {
        total_counts <- sum(Matrix::colSums(GetAssayData(gbm_subset, layer = "counts")))
    } else {
        total_counts <- sum(Matrix::colSums(GetAssayData(gbm_subset, assay = "Spatial", layer = "counts")))
    }
    n_cells <- ncol(gbm_subset)
    is_low_count <- total_counts < low_count_threshold
    message(sprintf("Batch %s: total counts = %d, cells = %d", prefix, total_counts, n_cells))
    if (is_low_count) {
        warning(sprintf("Batch %s has low total counts (%d)", prefix, total_counts))
    }

    # 5. Return counts info
    return(data.frame(
        batch = prefix,
        total_counts = total_counts,
        n_cells = n_cells,
        low_counts = is_low_count
    ))
}

# --- Process all batches ---
# Process batches sequentially to avoid memory explosion
batch_info_list <- list()
for (i in seq_along(batch_file_names)) {
    batch_info <- process_batch(i,
        batch_file_names = batch_file_names,
        batch_names = batch_names,
        raw_data_dir = raw_data_dir,
        rds_data_dir = rds_data_dir,
        low_count_threshold = low_count_threshold,
        VisiumHD = VisiumHD
    )
    batch_info_list[[i]] <- batch_info

    # Force garbage collection between batches
    gc()
}
batch_counts_df <- do.call(rbind, batch_info_list)

# --- Save batch counts summary ---
write.csv(batch_counts_df, paste0(rds_data_dir, "batch_counts_summary.csv"), row.names = FALSE)

# ===========================================
#           Mitochondrial QC Analysis
# ===========================================
message("=== Mitochondrial QC Analysis ===")

# Function to calculate and plot mitochondrial percentages
plot_mitochondrial_qc <- function(batch_names, rds_data_dir, output_dirs) {
    mito_qc_list <- list()

    for (batch_name in batch_names) {
        message("Calculating mitochondrial QC for batch: ", batch_name)

        # Load batch object
        gbm_subset <- readRDS(paste0(rds_data_dir, batch_name, "_raw.rds"))

        # Calculate mitochondrial percentages
        # For mouse: genes starting with "mt-", for human: genes starting with "MT-"
        mito_genes <- grep("^[Mm][Tt]-", rownames(gbm_subset), value = TRUE)
        if (length(mito_genes) > 0) {
            # Handle assay selection based on VisiumHD flag
            gbm_subset <- PercentageFeatureSet(gbm_subset, pattern = "^[Mm][Tt]-|^MT-", col.name = "percent.mt", assay = ifelse(VisiumHD, "Spatial.008um", "Spatial"))

            # Store QC data
            mito_qc_list[[batch_name]] <- data.frame(
                batch = batch_name,
                cell_id = colnames(gbm_subset),
                percent_mt = gbm_subset$percent.mt,
                n_genes = ifelse(VisiumHD, gbm_subset$nFeature_Spatial.008um, gbm_subset$nFeature_Spatial),
                n_counts = ifelse(VisiumHD, gbm_subset$nCount_Spatial.008um, gbm_subset$nCount_Spatial)
            )
        } else {
            message("No mitochondrial genes found for batch: ", batch_name)
            mito_qc_list[[batch_name]] <- data.frame(
                batch = batch_name,
                cell_id = character(0),
                percent_mt = numeric(0),
                n_genes = numeric(0),
                n_counts = numeric(0)
            )
        }
    }

    # Combine all batches
    mito_qc_df <- do.call(rbind, mito_qc_list)

    if (nrow(mito_qc_df) > 0) {
        # Create mitochondrial percentage plot
        p_mito <- ggplot(mito_qc_df, aes(x = batch, y = percent_mt, fill = batch)) +
            geom_violin(alpha = 0.7) +
            geom_boxplot(width = 0.1, fill = "white", alpha = 0.8) +
            labs(
                title = "Mitochondrial Gene Expression Percentage by Batch",
                subtitle = "Use this plot to determine mitochondrial cutoff threshold",
                x = "Batch",
                y = "Mitochondrial Gene %"
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

        # Create summary statistics
        mito_summary <- mito_qc_df |>
            dplyr::group_by(batch) |>
            dplyr::summarise(
                mean_mt = round(mean(percent_mt, na.rm = TRUE), 2),
                median_mt = round(median(percent_mt, na.rm = TRUE), 2),
                q75_mt = round(quantile(percent_mt, 0.75, na.rm = TRUE), 2),
                q95_mt = round(quantile(percent_mt, 0.95, na.rm = TRUE), 2),
                max_mt = round(max(percent_mt, na.rm = TRUE), 2),
                .groups = "drop"
            )

        write.csv(mito_summary, paste0(rds_data_dir, "mitochondrial_summary.csv"), row.names = FALSE)

        message("Mitochondrial QC plots saved to: ", paste0(output_dirs$clustering, "mitochondrial_qc_by_batch.pdf"))
        message("Mitochondrial summary saved to: ", paste0(rds_data_dir, "mitochondrial_summary.csv"))

        return(mito_qc_df)
    } else {
        message("No mitochondrial data found across all batches")
        return(NULL)
    }
}

# Generate mitochondrial QC plots
mito_qc_data <- plot_mitochondrial_qc(batch_names, rds_data_dir, output_dirs)

# --- Optional: Plot batch counts summary ---
# Create visualization of total counts with cell numbers
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


# --- Optional: Plot replicate PCA (code from original) ---
plot_replicate_pca <- function(batch_names, inspect_replicates, rds_data_dir, max_genes = 2000) {
    # Calculate per-batch average expression (pseudo-bulk)
    batch_expr <- sapply(batch_names, function(b) {
        gbm_subset <- readRDS(paste0(rds_data_dir, b, "_raw.rds"))
        # Use counts layer since data layer is empty at this stage
        counts_data <- GetAssayData(gbm_subset, layer = "counts")
        if (ncol(counts_data) == 0) {
            return(rep(NA, nrow(counts_data)))
        }
        Matrix::rowMeans(counts_data)
    })

    # Convert to matrix and transpose
    batch_expr <- as.matrix(batch_expr)
    batch_expr <- t(batch_expr)

    # Remove batches with NA or zero variance genes
    valid_batches <- complete.cases(batch_expr)
    if (sum(valid_batches) < 2) {
        cat("Warning: Not enough valid batches for PCA analysis\n")
        return(NULL)
    }

    batch_expr <- batch_expr[valid_batches, ]

    # Remove genes with zero variance across batches
    gene_vars <- apply(batch_expr, 2, var, na.rm = TRUE)
    valid_genes <- !is.na(gene_vars) & gene_vars > 0

    if (sum(valid_genes) < 10) {
        cat("Warning: Not enough variable genes for PCA analysis\n")
        return(NULL)
    }

    batch_expr <- batch_expr[, valid_genes]

    # Downsample genes if too many to avoid memory issues
    if (ncol(batch_expr) > max_genes) {
        message("Downsampling from ", ncol(batch_expr), " to ", max_genes, " genes for PCA")
        set.seed(42) # For reproducibility
        selected_genes <- sample(ncol(batch_expr), max_genes)
        batch_expr <- batch_expr[, selected_genes]
    }

    # PCA
    pca <- prcomp(batch_expr, scale. = TRUE)
    pca_df <- data.frame(
        PC1 = pca$x[, 1],
        PC2 = pca$x[, 2],
        batch = rownames(pca$x)
    )

    # Assign group for coloring
    group_map <- setNames(rep(NA, length(batch_names)), batch_names)
    for (i in seq_along(inspect_replicates)) {
        group_map[inspect_replicates[[i]]] <- paste0("Group", i)
    }
    pca_df$group <- group_map[pca_df$batch]

    # Plot
    ggplot2::ggplot(pca_df, aes(x = PC1, y = PC2, color = group, label = batch)) +
        geom_point(size = 4) +
        geom_text(vjust = -0.7, size = 3) +
        labs(
            title = paste("PCA of batch replicates (", ncol(batch_expr), "genes )"),
            x = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 1), "%)"),
            y = paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 1), "%)")
        ) +
        theme_minimal() +
        theme(legend.position = "bottom")
}

if (!is.null(inspect_replicates) && length(inspect_replicates) > 0) {
    p_replicate <- plot_replicate_pca(batch_names, inspect_replicates, rds_data_dir)
    if (!is.null(p_replicate)) {
        ggsave(paste0(output_dirs$clustering, output.file.prefix, "_replicate_pca.pdf"), plot = p_replicate)
    }
    message("Replicate PCA plotting logic would go here if function is available.")
}


# ===========================================
#           Sketch Sampling (On-Disk)
# ===========================================
message("=== Sketch Sampling ===")

brain_data_list <- list()
for (i in seq_along(batch_file_names)) { # Keep loop for sequential sketching
    prefix <- batch_names[i]
    message("Sketching batch: ", prefix)

    # Check if sketched RDS file already exists
    sketched_rds_path <- paste0(rds_data_dir, prefix, ifelse(VisiumHD, "_sketched.rds", "_raw.rds"))

    if (file.exists(sketched_rds_path)) {
        message("Sketched RDS file already exists for batch ", prefix, ". Skipping sketching...")
        brain_data_list[[prefix]] <- sketched_rds_path
        next
    }

    # 1. Load raw object
    gbm_subset <- readRDS(paste0(rds_data_dir, prefix, "_raw.rds"))

    # 2. Set default assay for VisiumHD
    if (VisiumHD) {
        Seurat::DefaultAssay(gbm_subset) <- "Spatial.008um"
    }

    # 3. Normalize and find variable features (in-memory steps)
    gbm_subset <- Seurat::NormalizeData(gbm_subset) |> Seurat::FindVariableFeatures()

    # 4. Sketch Data (On-Disk)
    if (VisiumHD) {
        # Perform sketching
        tryCatch(
            {
                gbm_subset <- Seurat::SketchData(
                    object = gbm_subset,
                    ncells = 12000,
                    method = "Uniform", # Or "LeverageScore"
                    sketched.assay = "sketch"
                )

                # Switch analysis to sketched cells
                Seurat::DefaultAssay(gbm_subset) <- "sketch"

                # Save to RDS format only
                saveRDS(gbm_subset, sketched_rds_path)
                message("Saved sketched object to RDS: ", sketched_rds_path)
                brain_data_list[[prefix]] <- sketched_rds_path
            },
            error = function(e) {
                stop("Failed to sketch batch ", prefix, ": ", e$message)
            }
        )
    } else {
        # For non-VisiumHD, save the processed object without sketching
        saveRDS(gbm_subset, sketched_rds_path)
        brain_data_list[[prefix]] <- sketched_rds_path
    }
}


# ===========================================
#           Merge Sketched Objects in Chunks (Memory-Efficient)
# ===========================================
message("=== Merging Sketched Objects in Chunks ===")

# Define chunk size and create chunks
batch_chunks <- split(batch_names, ceiling(seq_along(batch_names) / on_disk$merge_n_batches))
message("Processing ", length(batch_chunks), " chunks with up to ", on_disk$merge_n_batches, " batches each")

chunk_paths <- list()
for (i in seq_along(batch_chunks)) {
    chunk_name <- paste0("chunk_", i)
    chunk_rds_path <- paste0(rds_data_dir, output.file.prefix, "_", chunk_name, "_merged.rds")

    if (file.exists(chunk_rds_path)) {
        message("Chunk ", i, " already exists. Skipping merge...")
        chunk_paths[[chunk_name]] <- chunk_rds_path
        next
    }

    message("=== Processing Chunk ", i, " (", length(batch_chunks[[i]]), " batches) ===")

    # Load first batch in chunk
    first_batch <- batch_chunks[[i]][1]
    merged_chunk <- readRDS(brain_data_list[[first_batch]])
    message("Loaded first batch in chunk: ", first_batch)

    # Clean up unnecessary slots
    merged_chunk@tools <- list()
    merged_chunk@misc <- list()
    gc()

    # Merge remaining batches in chunk
    if (length(batch_chunks[[i]]) > 1) {
        for (batch_name in batch_chunks[[i]][-1]) {
            message("Merging batch: ", batch_name)
            obj_to_merge <- readRDS(brain_data_list[[batch_name]])

            # Clean up unnecessary slots
            obj_to_merge@tools <- list()
            obj_to_merge@misc <- list()

            # Merge
            merged_chunk <- merge(merged_chunk, y = obj_to_merge)

            # Immediate cleanup
            rm(obj_to_merge)
            gc()
        }
    }

    # Save chunk
    saveRDS(merged_chunk, chunk_rds_path)
    chunk_paths[[chunk_name]] <- chunk_rds_path
    message("Saved chunk ", i, " to: ", chunk_rds_path)

    # Clean up chunk from memory
    rm(merged_chunk)
    gc()
}

# ===========================================
#           Extract Data for Integration (Memory-Efficient)
# ===========================================
message("=== Extracting Data for Integration ===")

# Function to extract essential data from a chunk
extract_essential_data <- function(chunk_path, chunk_name) {
    message("Extracting data from: ", chunk_name)

    # Load chunk
    chunk_obj <- readRDS(chunk_path)

    # Set correct assay
    if (VisiumHD) {
        DefaultAssay(chunk_obj) <- "sketch"
    } else {
        DefaultAssay(chunk_obj) <- "Spatial"
    }

    # Join layers before normalization to avoid v5 issues
    chunk_obj <- JoinLayers(chunk_obj)

    # Normalize and find variable features
    chunk_obj <- chunk_obj |>
        NormalizeData() |>
        FindVariableFeatures()

    # Extract essential components
    essential_data <- list(
        data = GetAssayData(chunk_obj, layer = "data"),
        var_features = VariableFeatures(chunk_obj),
        meta_data = chunk_obj@meta.data,
        cell_names = colnames(chunk_obj),
        chunk_name = chunk_name
    )

    # Clean up
    rm(chunk_obj)
    gc()

    return(essential_data)
}

# Extract data from all chunks
essential_data_list <- list()
for (i in seq_along(chunk_paths)) {
    chunk_name <- names(chunk_paths)[i]
    essential_data_list[[chunk_name]] <- extract_essential_data(chunk_paths[[i]], chunk_name)
}

# ===========================================
#           Create Combined Object for Integration
# ===========================================
message("=== Creating Combined Object for Integration ===")

# Combine essential data
all_var_features <- unique(unlist(lapply(essential_data_list, function(x) x$var_features)))
all_meta_data <- do.call(rbind, lapply(essential_data_list, function(x) x$meta_data))
all_cell_names <- unlist(lapply(essential_data_list, function(x) x$cell_names))

# Create a minimal Seurat object for integration
message("Creating minimal combined object...")

# Combine normalized data matrices (these are already normalized)
combined_data <- do.call(cbind, lapply(essential_data_list, function(x) x$data))

# Create minimal Seurat object with normalized data
minimal_obj <- CreateSeuratObject(
    counts = combined_data, # Using normalized data as counts
    meta.data = all_meta_data,
    assay = ifelse(VisiumHD, "sketch", "Spatial")
)

# Since we're using normalized data as counts, we need to copy it to the data layer
minimal_obj <- SetAssayData(minimal_obj, layer = "data", new.data = combined_data)

# Set variable features
VariableFeatures(minimal_obj) <- all_var_features

# Clean up intermediate data
rm(essential_data_list, combined_data)
gc()

# ===========================================
#           Integration & Clustering (Minimal Object)
# ===========================================
message("=== Integration & Clustering (Minimal Object) ===")

# Since data is already normalized, skip normalization and go directly to scaling
minimal_obj <- minimal_obj |>
    ScaleData() |>
    RunPCA(reduction.name = "PCA")

# Join layers before integration
minimal_obj <- JoinLayers(minimal_obj)

# Check if the object has the required assay and data for integration

# Ensure we're using the correct assay for integration

# Modified integration with error checking
integration_success <- FALSE # Initialize variable
tryCatch(
    {
        # Ensure the object has valid assays before integration
        if (length(minimal_obj@assays) > 0) {
            # Get the first (and likely only) assay name if "sketch" doesn't exist
            assay_names <- names(minimal_obj@assays)
            target_assay <- ifelse(VisiumHD, "sketch", "Spatial")
            if (target_assay %in% assay_names) {
                DefaultAssay(minimal_obj) <- target_assay
            } else {
                # Use the first available assay
                DefaultAssay(minimal_obj) <- assay_names[1]
                message("Using '", assay_names[1], "' assay for integration as '", target_assay, "' not found")
            }

            # Check if the assay data exists before integration
            if (!is.null(minimal_obj[[DefaultAssay(minimal_obj)]])) {
                minimal_obj <- IntegrateLayers(
                    minimal_obj,
                    method = RPCAIntegration,
                    orig.reduction = "PCA",
                    new.reduction = ifelse(VisiumHD, "sketch_integrated.rpca", "integrated.rpca"),
                    verbose = TRUE
                )
                integration_success <- TRUE
            } else {
                message("Assay data is NULL, skipping integration")
            }
        } else {
            stop("No valid assays found in minimal object for integration")
        }
    },
    error = function(e) {
        message("Error during integration: ", e$message)
        message("Attempting alternative integration approach")
        integration_success <- FALSE
    }
)

# Clustering
if (integration_success && ifelse(VisiumHD, "sketch_integrated.rpca", "integrated.rpca") %in% names(minimal_obj@reductions)) {
    target_reduction <- ifelse(VisiumHD, "sketch_integrated.rpca", "integrated.rpca")
    target_assay <- ifelse(VisiumHD, "sketch", "Spatial")
} else {
    message("Using sketch_PCA for downstream analysis due to integration failure")
    target_reduction <- "PCA"
    target_assay <- ifelse(VisiumHD, "sketch", "Spatial")
}

minimal_obj <- FindNeighbors(
    minimal_obj,
    assay = target_assay,
    reduction = target_reduction,
    dims = 1:30,
    k.param = Cluster_n_neighbors
) |>
    FindClusters(cluster.name = "seurat_cluster", resolution = Cluster_resolution)

# UMAP
minimal_obj <- RunUMAP(
    minimal_obj,
    reduction = target_reduction,
    reduction.name = "umap",
    return.model = TRUE,
    dims = 1:30,
    n.components = 2
)

# 3D UMAP
minimal_obj <- RunUMAP(
    minimal_obj,
    reduction = target_reduction,
    reduction.name = "umap_3D",
    return.model = TRUE,
    dims = 1:30,
    n.components = 3
)

# Create a copy of umap as umap_integrated to maintain consistency if needed
if (!integration_success && !"umap_integrated" %in% names(minimal_obj@reductions)) {
    minimal_obj[["umap_integrated"]] <- minimal_obj[["umap"]]
}

# Extract integration results
integration_results <- list(
    reductions = minimal_obj@reductions,
    clusters = minimal_obj$seurat_cluster,
    cell_names = colnames(minimal_obj)
)

# Save integration_results to RDS
saveRDS(integration_results, paste0(rds_data_dir, output.file.prefix, "_umap_clusters.rds"))

message("=== Cleaning Up ===")
# Clean up minimal object
rm(minimal_obj, integration_results)
gc()


# ================================================
#     Method to Transfer Results Back to Chunks
# ================================================
add_results_to_chunk <- function(chunk_obj, integration_results) {
    # Join layers to avoid v5 issues
    chunk_obj <- JoinLayers(chunk_obj)

    # Set correct assay based on VisiumHD flag
    if (VisiumHD) {
        DefaultAssay(chunk_obj) <- "sketch"
    } else {
        DefaultAssay(chunk_obj) <- "Spatial"
    }

    # Get cell indices for this chunk
    chunk_cells <- colnames(chunk_obj)
    cell_indices <- match(chunk_cells, integration_results$cell_names)

    # Add reductions
    for (reduction_name in names(integration_results$reductions)) {
        chunk_obj@reductions[[reduction_name]] <- integration_results$reductions[[reduction_name]]
        # Subset embeddings to this chunk's cells
        chunk_obj@reductions[[reduction_name]]@cell.embeddings <-
            integration_results$reductions[[reduction_name]]@cell.embeddings[cell_indices, , drop = FALSE]
    }

    # Add clusters
    chunk_obj$seurat_cluster <- integration_results$clusters[cell_indices]

    # Return updated chunk object
    return(chunk_obj)
}

# Load chunk
# output_path <- paste0(rds_data_dir, output.file.prefix, "_chunk_", 1, "_merged.rds")
# merged_obj <- readRDS(output_path)
# integration_results <- readRDS(paste0(rds_data_dir, output.file.prefix, "_umap_clusters.rds"))
# merged_obj <- add_results_to_chunk(merged_obj, integration_results)
# =============================================




# ===========================================
#       Project Data Back to Full Resolution (Optional/On-Demand)
# ===========================================
message("=== Projecting Data (Optional) ===")
# This step is memory intensive. It's often done separately or on-demand.
# It requires the original full-resolution data and the sketched model.

# --- Function for projecting a single batch ---
project_single_batch <- function(batch_name, rds_data_dir, clustered_sketched_rds_path, output.file.prefix) {
    message("Projecting batch: ", batch_name)

    # 1. Load the clustered sketched object (contains the UMAP model)
    clustered_sketched_obj <- readRDS(clustered_sketched_rds_path)

    # 2. Load the original raw full-resolution object
    full_res_obj <- readRDS(paste0(rds_data_dir, batch_name, "_raw.rds"))
    if (VisiumHD) {
        Seurat::DefaultAssay(full_res_obj) <- "Spatial.008um"
    } else {
        Seurat::DefaultAssay(full_res_obj) <- "Spatial"
    }
    # Normalize full res data (ensure consistency)
    full_res_obj <- NormalizeData(full_res_obj)

    # 3. Project Data from sketched to full resolution
    # This requires the sketched object with the model and the full-res object
    projected_obj <- Seurat::ProjectData(
        object = full_res_obj, # Target: Full resolution
        assay = "Spatial.008um", # Assay in full_res_obj
        # Provide the reductions and model from the clustered sketched object
        full.reduction = "full.sketch_integrated.rpca", # Name for the projected PCA
        sketched.assay = clustered_sketched_obj[["sketch"]], # Sketched assay data
        sketched.reduction = clustered_sketched_obj[["sketch_PCA"]], # Sketched PCA
        umap.model = clustered_sketched_obj[["sketch_umap"]], # 2D UMAP model
        dims = 1:50, # Use dims from integration
        refdata = list(projected_seurat_cluster = "seurat_cluster") # Transfer clusters
    )

    # 4. Save the projected object for this batch
    projected_batch_rds_path <- paste0(rds_data_dir, batch_name, "_projected.rds")
    saveRDS(projected_obj, projected_batch_rds_path)
    message("Saved projected object for batch ", batch_name, " to RDS: ", projected_batch_rds_path)

    # Clean up
    rm(full_res_obj, projected_obj)
    gc()
}

# --- Execute Projection for All Batches (Memory Intensive) ---
# Uncomment the lines below if you want to run projection immediately.
# It's recommended to run this in a separate script or session due to memory demands.

# lapply(batch_names, project_single_batch,
#        rds_data_dir = rds_data_dir,
#        clustered_sketched_rds_path = clustered_sketched_rds_path,
#        output.file.prefix = output.file.prefix)

message("Projection step completed (or skipped).")
message("===========================================")
message("Done preprocessing and clustering (sketched).")
message("Projected data step is optional and memory-intensive.")
message("===========================================")
