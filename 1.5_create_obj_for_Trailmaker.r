source("util_headers.r")

merged_obj <- readRDS(paste0(rds_data_dir, output.file.prefix, "_merged.rds"))
# Validate cluster method and load if necessary
if (cluster_method %in% colnames(merged_obj@meta.data)) {
    valid_cluster_method <- cluster_method
} else {
    # Load cluster assignments from CSV file
    cluster_csv_path <- paste0(output_dirs$clustering, cluster_method, ".csv")
    if (file.exists(cluster_csv_path)) {
        cluster_data <- read.csv(cluster_csv_path, row.names = 1)
        cluster_vec <- setNames(cluster_data[[1]], rownames(cluster_data))
        merged_obj@meta.data[[cluster_method]] <- cluster_vec[rownames(merged_obj@meta.data)]
        valid_cluster_method <- cluster_method
    } else {
        stop(paste("Cluster method", cluster_method, "not found in metadata or CSV."))
    }
}

Idents(merged_obj) <- valid_cluster_method

# Control variable: TRUE for spatial coordinates, FALSE for existing reductions
use_spatial_coords <- TRUE

# 1. Only keep the RNA assay (with raw counts)
if ("RNA" %in% names(merged_obj@assays)) {
    raw_counts <- Seurat::GetAssayData(merged_obj[["RNA"]], layer = "counts")
} else if ("Spatial" %in% names(merged_obj@assays)) {
    warning("RNA assay not found, using raw counts from Spatial assay")
    raw_counts <- Seurat::GetAssayData(merged_obj[["Spatial"]], layer = "counts")
} else {
    stop("Neither RNA nor Spatial assay found in merged_obj")
}

scdata <- Seurat::CreateSeuratObject(
    counts = raw_counts,
    meta.data = merged_obj@meta.data
)
# Convert to v4 format immediately after creation
scdata@assays$RNA <- as(object = scdata@assays$RNA, Class = "Assay")

# 2. Add sample assignment (samples slot)
if ("batch" %in% colnames(merged_obj@meta.data)) {
    scdata$samples <- merged_obj@meta.data$batch
} else {
    scdata$samples <- rep("sample1", ncol(scdata))
}

# 3. Copy cluster metadata (already in meta.data)

# 4. Copy PCA reduction
if ("PCA" %in% names(merged_obj@reductions)) {
    pca_red <- merged_obj[["PCA"]]
    pca_red@assay.used <- "RNA"
    # Ensure v4 compatibility for reduction
    # pca_red <- as(object = pca_red, Class = "DimReduc")
    scdata[["pca"]] <- pca_red
}

# 5. Copy UMAP reduction or create from spatial coordinates
if (use_spatial_coords && ("Spatial" %in% names(merged_obj@assays) || !is.null(merged_obj@images))) {
    # Handle multiple spatial fields of view with coordinate shifting
    if ("batch" %in% colnames(merged_obj@meta.data) && length(unique(merged_obj$batch)) > 1) {
        # Get all unique batch names
        batch_names <- unique(merged_obj$batch)
        n_batches <- length(batch_names)

        if (n_batches >= 2 && n_batches <= 8) {
            # Extract coordinates for each batch
            batch_coords <- list()
            # Initialize variables to track max shifts
            max_x_shift <- 0
            max_y_shift <- 0

            for (batch_name in batch_names) {
                batch_mask <- merged_obj$batch == batch_name
                batch_obj <- subset(merged_obj, cells = colnames(merged_obj)[batch_mask])
                batch_coord <- Seurat::GetTissueCoordinates(batch_obj)
                batch_coords[[batch_name]] <- batch_coord

                # Calculate max shifts during the loop for efficiency
                batch_coord_matrix <- as.matrix(batch_coord[, 1:2])
                max_x_shift <- max(max_x_shift, max(batch_coord_matrix[, 1]))
                max_y_shift <- max(max_y_shift, max(batch_coord_matrix[, 2]))
            }

            # Apply shifts to coordinates using batch_coords
            spatial_coords <- do.call(rbind, batch_coords)
            rownames(spatial_coords) <- spatial_coords$cell
            # Ensure the order matches the original data
            spatial_coords <- spatial_coords[match(colnames(merged_obj), rownames(spatial_coords)), ]
            spatial_coords <- as.matrix(spatial_coords[, 1:2])

            # Define grid size for up to 8 batches
            grid_sizes <- list(
                `2` = c(2, 1),
                `3` = c(2, 2),
                `4` = c(2, 2),
                `5` = c(3, 2),
                `6` = c(3, 2),
                `7` = c(4, 2),
                `8` = c(4, 2)
            )

            grid_dims <- grid_sizes[[as.character(n_batches)]]
            n_cols <- grid_dims[1]
            n_rows <- grid_dims[2]

            # Calculate shifts for each batch
            shifts <- list()
            for (i in seq_along(batch_names)) {
                batch_name <- batch_names[i]
                row <- (i - 1) %/% n_cols
                col <- (i - 1) %% n_cols
                shifts[[batch_name]] <- c(col * max_x_shift * 1.05, row * max_y_shift * -1.05)
            }

            for (batch_name in batch_names) {
                batch_mask <- merged_obj$batch == batch_name
                shift <- shifts[[batch_name]]
                spatial_coords[batch_mask, 1] <- spatial_coords[batch_mask, 1] + shift[1]
                spatial_coords[batch_mask, 2] <- spatial_coords[batch_mask, 2] + shift[2]
            }
        } else {
            stop(paste("Only 2-8 batches are supported, got", n_batches))
        }
    } else {
        # Single batch - extract spatial coordinates normally
        spatial_coords <- Seurat::GetTissueCoordinates(merged_obj)
        spatial_coords <- as.matrix(spatial_coords[, 1:2])
    }

    # Convert to matrix and standardize column names
    colnames(spatial_coords) <- c("umap_1", "umap_2")

    # Create a proper DimReduc object matching the structure
    spatial_red <- new("DimReduc",
        cell.embeddings = spatial_coords,
        assay.used = "RNA",
        global = TRUE,
        key = "umap_"
    )

    # Convert to v4 compatibility
    spatial_red <- as(object = spatial_red, Class = "DimReduc")
    scdata[["umap"]] <- spatial_red
}

if (!use_spatial_coords) {
    # Use existing reductions
    umap_name <- grep("umap", names(merged_obj@reductions), value = TRUE)
    if (length(umap_name) > 0) {
        umap_red <- merged_obj[[umap_name[1]]]
        umap_red@assay.used <- "RNA"
        # Ensure v4 compatibility for reduction
        umap_red <- as(object = umap_red, Class = "DimReduc")
        scdata[["umap"]] <- umap_red
    } else if ("tsne" %in% names(merged_obj@reductions)) {
        tsne_red <- merged_obj[["tsne"]]
        tsne_red@assay.used <- "RNA"
        # Ensure v4 compatibility for reduction
        tsne_red <- as(object = tsne_red, Class = "DimReduc")
        scdata[["tsne"]] <- tsne_red
    } else {
        warning("No UMAP or tSNE reduction found in merged_obj")
    }
}

# 6. Remove all other assays and reductions (keep only RNA, pca, umap/tsne)
scdata@assays <- scdata@assays["RNA"]
scdata@reductions <- scdata@reductions[intersect(names(scdata@reductions), c("pca", "umap", "tsne"))]

# Seurat::DimPlot(scdata, reduction = "umap", group.by = "samples")

# 7. Save the object for Trailmaker (object is already in v4 format)
saveRDS(scdata, file = paste0(
    project_dir, "rds_data/", output.file.prefix,
    ifelse(use_spatial_coords, "_spatial", ""), "_forTrailmaker.rds"
))
