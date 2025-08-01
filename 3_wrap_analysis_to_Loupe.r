# Run this script only after cell type deconvolution and `predict_GBM_scSeq_cellType_X` are written to the object


# remotes::install_github("10xGenomics/loupeR")
loupeR::setup()

# ===========================================
#           Define Tumour Clusters
# ***********************************************
# Please identify the cluster numbers
# in the Plotly UMAP html
# ***********************************************

# ===========================================
#           Setup and Configuration
# ===========================================
source("util_headers.r")

merged_obj <- readRDS(paste0(project_dir, "rds_data/", output.file.prefix, "_clustered_12k.rds"))

# ===========================================
# load cell type prediction results
# ===========================================
cellType_1_cluster_df <- read.csv(paste0(output_dirs$clustering, "predict_GBM_scSeq_cellType_1.csv"), stringsAsFactors = FALSE, check.names = FALSE)
rownames(cellType_1_cluster_df) <- cellType_1_cluster_df$Cell
cellType_1_cluster_df$Cell <- NULL

cellType_2_cluster_df <- read.csv(paste0(output_dirs$clustering, "predict_GBM_scSeq_cellType_2.csv"), stringsAsFactors = FALSE, check.names = FALSE)
rownames(cellType_2_cluster_df) <- cellType_2_cluster_df$Cell
cellType_2_cluster_df$Cell <- NULL

merged_obj$predict_GBM_scSeq_cellType_1 <- cellType_1_cluster_df
merged_obj$predict_GBM_scSeq_cellType_2 <- cellType_2_cluster_df

merged_obj <- JoinLayers(merged_obj)
# ===========================================

Idents(merged_obj) <- valid_cluster_method
cell_region_cluster <- Idents(merged_obj)
merged_obj$cell_region_cluster <- dplyr::case_when(
  cell_region_cluster %in% outTumour_cluster_nums_vector ~ "outTumour",
  cell_region_cluster %in% inTumour_cluster_nums_vector ~ "inTumour",
  cell_region_cluster %in% edgeTumour_cluster_nums_vector ~ "edgeTumour",
  .default = as.character(cell_region_cluster)
) |> as.factor()

merged_obj_list <- SplitObject(merged_obj, split.by = "batch")

# ===========================
# single files
# ===========================
for (batch_id in batch_names) {
  one_subset <- merged_obj_list[[batch_id]]

  DefaultAssay(one_subset) <- ifelse(VisiumHD, "Spatial.008um", "Spatial")
  Idents(one_subset) <- valid_cluster_method

  count_mat <- one_subset[[ifelse(VisiumHD, "Spatial.008um", "Spatial")]]$counts

  # ===========================
  # extract clusters
  # ===========================
  clusters <- tibble::tibble(
    projected_seurat_cluster = Idents(one_subset)[1:ncol(count_mat)],
    cell_region = one_subset$cell_region_cluster[1:ncol(count_mat)],
    predict_GBM_scSeq_cellType1 = one_subset$predict_GBM_scSeq_cellType_1[1:ncol(count_mat)],
    predict_GBM_scSeq_cellType2 = one_subset$predict_GBM_scSeq_cellType_2[1:ncol(count_mat)]
  ) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::starts_with("predict_"),
        ~ as.factor(tidyr::replace_na(as.character(.), "null"))
      )
    ) |>
    as.list()

  # ===========================
  # extract tissue coordinates
  # ===========================
  st_location <- GetTissueCoordinates(one_subset)
  st_location$cell <- NULL
  if (VisiumHD) {
    dim_name_vector <- c("fullsketchumap_1", "fullsketchumap_2")
  } else {
    dim_name_vector <- c("umap_1", "umap_2")
  }
  projections <- list(
    full_sketch_umap = one_subset[[ifelse(VisiumHD, "full.sketch_umap", "umap")]]@cell.embeddings[, dim_name_vector] |> as.matrix(),
    spatial = st_location |> as.matrix()
  )

  loupeR::create_loupe(
    count_mat = count_mat,
    clusters = clusters,
    projections = projections,
    output_dir = output_dirs$loupe,
    output_name = paste0(batch_id, "_manual")
  )
}

# ===========================
# merged files
# ===========================
DefaultAssay(merged_obj) <- ifelse(VisiumHD, "Spatial.008um", "Spatial")
Idents(merged_obj) <- valid_cluster_method

# The Idents(merged_obj) is the concatenation of sample "SAL", "KD", "RAD", "KD_RAD",
# within each sample is also the resulting concatenation of 'Spatial.008um' and 'Spatial.016um',
# e.g. length(merged_obj[['Spatial.008um']]$SAL)+length(merged_obj[['Spatial.016um']]$SAL)
# +length(merged_obj[['Spatial.008um']]$KD)+length(merged_obj[['Spatial.016um']]$KD)+...

# ===========================
# Initialize empty factors with correct levels
initialize_factors <- function(merged_obj) {
  list(
    cell_region = factor(levels = levels(merged_obj$cell_region_cluster)),
    predict_GBM_scSeq_cellType1 = factor(levels = levels(merged_obj$predict_GBM_scSeq_cellType_1)),
    predict_GBM_scSeq_cellType2 = factor(levels = levels(merged_obj$predict_GBM_scSeq_cellType_2)),
    batch = factor(levels = levels(merged_obj$batch)),
    idents = factor(levels = levels(Idents(merged_obj)))
  )
}

# ===========================
# process batches
# ===========================
process_batches <- function(merged_obj, merged_obj_list) {
  new_factors <- initialize_factors(merged_obj)
  current_pos <- 1
  counts_list <- list()

  # ===========================
  # interpolate counts
  # ===========================
  for (batch_id in batch_names) {
    one_subset <- merged_obj_list[[batch_id]]

    counts_008um <- one_subset[[ifelse(VisiumHD, "Spatial.008um", "Spatial")]]$counts
    len_008 <- ncol(counts_008um)
    if (VisiumHD) {
      len_016 <- ncol(one_subset[["Spatial.016um"]]$counts)
    }

    counts_list[[batch_id]] <- counts_008um

    idents_range <- current_pos:(current_pos + len_008 - 1)

    # Update factors
    new_factors$idents <- c(new_factors$idents, Idents(merged_obj)[idents_range])
    new_factors$batch <- c(new_factors$batch, merged_obj$batch[idents_range])
    new_factors$cell_region <- c(
      new_factors$cell_region,
      merged_obj$cell_region_cluster[idents_range]
    )
    new_factors$predict_GBM_scSeq_cellType1 <- c(
      new_factors$predict_GBM_scSeq_cellType1,
      merged_obj$predict_GBM_scSeq_cellType_1[idents_range]
    )
    new_factors$predict_GBM_scSeq_cellType2 <- c(
      new_factors$predict_GBM_scSeq_cellType2,
      merged_obj$predict_GBM_scSeq_cellType_2[idents_range]
    )

    # Update position counter
    if (VisiumHD) {
      current_pos <- current_pos + len_008 + len_016
    } else {
      current_pos <- current_pos + len_008
    }
  }

  # Combine counts
  combined_counts <- do.call(cbind, counts_list)

  # Create final clusters list with proper factor conversion and NA handling
  clusters <- as.data.frame(new_factors) |>
    dplyr::mutate(
      projected_seurat_cluster = as.factor(idents),
      sample_groups = as.factor(batch),
      cell_region = as.factor(cell_region),
      dplyr::across(
        dplyr::starts_with("predict_"),
        ~ as.factor(tidyr::replace_na(as.character(.), "null"))
      )
    ) |>
    dplyr::select(-idents, -batch) |>
    as.list()

  # Return both counts and clusters
  list(
    combined_counts = combined_counts,
    clusters = clusters
  )
}

result <- process_batches(merged_obj, merged_obj_list)
combined_counts <- result$combined_counts
clusters <- result$clusters

# ===========================
# form square tile coordinates
# ===========================
form_square_tile_coor <- function(merged_obj_list) {
  st_location_list <- list()
  for (batch_id in batch_names) {
    one_subset <- merged_obj_list[[batch_id]]
    st_location <- GetTissueCoordinates(one_subset)
    st_location$cell <- NULL
    st_location_list[[batch_id]] <- st_location
  }

  # ===========================
  # shift coordinates
  # ===========================
  # Calculate dynamic shift values based on the minimum and maximum coordinates
  max_x_shift <- max(st_location_list[[1]][, 1])
  max_y_shift <- max(st_location_list[[1]][, 2])

  shifts <- list()

  if (length(st_location_list) == 2) {
    shifts <- list(
      a = c(-max_x_shift, 0), # Shift x left by dynamic value, y stays
      b = c(0, 0) # Shift x right, y stays
    )
  } else if (length(st_location_list) == 3) {
    shifts <- list(
      a = c(-max_x_shift, -max_y_shift), # Top-left
      b = c(0, -max_y_shift), # Top-right
      c = c(-max_x_shift, 0) # Bottom-left
    )
  } else if (length(st_location_list) == 4) {
    shifts <- list(
      a = c(-max_x_shift, -max_y_shift), # Top-left
      b = c(0, -max_y_shift), # Top-right
      c = c(-max_x_shift, 0), # Bottom-left
      d = c(0, 0) # Bottom-right
    )
  } else {
    stop("Only 2, 3, or 4 batches are supported")
  }
  names(shifts) <- names(st_location_list)

  for (plot_name in names(st_location_list)) {
    shift <- shifts[[plot_name]]
    st_location_list[[plot_name]][, 1] <- st_location_list[[plot_name]][, 1] + shift[1] # Shift x
    st_location_list[[plot_name]][, 2] <- st_location_list[[plot_name]][, 2] + shift[2] # Shift y
  }
  names(st_location_list) <- NULL # prevent the `do.call(rbind, )` adding the names to rownames
  combined_location <- do.call(rbind, st_location_list) |> as.matrix()
  return(combined_location)
}

st_location_merged <- form_square_tile_coor(merged_obj_list)

if (VisiumHD) {
  dim_name_vector <- c("fullsketchumap_1", "fullsketchumap_2")
} else {
  dim_name_vector <- c("umap_1", "umap_2")
}
projections <- list(
  full_sketch_umap = merged_obj[[ifelse(VisiumHD, "full.sketch_umap", "umap")]]@cell.embeddings[, dim_name_vector] |> as.matrix(),
  spatial = st_location_merged |> as.matrix()
)


loupeR::create_loupe(
  count_mat = combined_counts,
  clusters = clusters,
  projections = projections,
  output_dir = output_dirs$loupe,
  output_name = paste0(output.file.prefix, "_manual")
)

cat("=====================================\n")
cat("Done exporting to Loupe files\n")
cat("=====================================\n")

# ===========================
# feature list
# https://www.cellsignal.com/pathways/immune-cell-markers-mouse
# ===========================
Macrophages_genes <- grep("^NOS2$|^cd86$|^cd80$|^ADGRE1$|^CD163$|^CD206$|^ARG", rownames(merged_obj), ignore.case = T, value = T)

Tcell_genes <- grep("^cd[43]$|^cd8$|^cd62L$|^CD25$|^CD69$|^IL2$|^AHR|^CXCR5|^IL9|^SPI1|^BCL6|^TBX21|^IFNG|^IL4|^IL17|^FOXP3|^RORC|^GATA3", rownames(merged_obj), ignore.case = T, value = T)

data.frame(
  Name = Macrophages_genes,
  List = "Macrophages_markers"
) |> write.csv(paste0(output_dirs$loupe, "/Macrophages_markers.csv"), row.names = FALSE)

data.frame(
  Name = Tcell_genes,
  List = "Tcell_markers"
) |> write.csv(paste0(output_dirs$loupe, "/Tcell_markers.csv"), row.names = FALSE)

cat("=====================================\n")
cat("Done marker genes export\n")
cat("=====================================\n")
quit("no")
