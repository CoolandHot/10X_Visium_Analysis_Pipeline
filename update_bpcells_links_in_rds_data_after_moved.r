old_dir <- "/scratch/hh01116/10X_KD/rds_data/"
new_dir <- "/scratch/hh01116/10X_KD/KD_rds_data/"

library(BPCells)

update_bpcells_slots <- function(obj, old_dir, new_dir) {
    for (assay_name in names(obj@assays)) {
        assay <- obj@assays[[assay_name]]
        # For Assay5 objects, use @layers slot
        if ("layers" %in% slotNames(assay)) {
            # Update counts layer if it's a BPCells matrix and path starts with old_dir
            if ("counts" %in% names(assay@layers) && inherits(assay@layers[["counts"]]@matrix, "MatrixDir")) {
                bp_dir <- assay@layers[["counts"]]@matrix@dir
                if (startsWith(bp_dir, old_dir)) {
                    new_bp_dir <- sub(old_dir, new_dir, bp_dir, fixed = TRUE)
                    assay@layers[["counts"]]@matrix@dir <- new_bp_dir
                }
            }
            # Update data layer if it's a BPCells matrix and path starts with old_dir
            if ("data" %in% names(assay@layers) && inherits(assay@layers[["data"]]@matrix, "MatrixDir")) {
                bp_dir <- assay@layers[["data"]]@matrix@dir
                if (startsWith(bp_dir, old_dir)) {
                    new_bp_dir <- sub(old_dir, new_dir, bp_dir, fixed = TRUE)
                    assay@layers[["data"]]@matrix@dir <- new_bp_dir
                }
            }
        }
        obj@assays[[assay_name]] <- assay
    }
    return(obj)
}

rds_files <- list.files(new_dir, pattern = "\\.rds$", full.names = TRUE)
for (rds_path in rds_files) {
    cat("Processing:", rds_path, "\n")
    obj <- readRDS(rds_path)
    obj2 <- update_bpcells_slots(obj, old_dir, new_dir)
    saveRDS(obj2, rds_path)
}
cat("BPCells links updated in all RDS files.\n")
