old_dir <- "/scratch/hh01116/10X_KD/rds_data/"
new_dir <- "/scratch/hh01116/10X_KD/KD_rds_data/"

library(BPCells)

update_bpcells_slots <- function(obj, old_dir, new_dir) {
    for (assay_name in names(obj@assays)) {
        assay <- obj@assays[[assay_name]]
        # For Assay5 objects, use @layers slot
        if ("layers" %in% slotNames(assay)) {
            for (key in names(assay@layers)) {
                # Check if the layer is a RenameDims object (BPCells)
                if (inherits(assay@layers[[key]], "RenameDims")) {
                    if ("dir" %in% slotNames(assay@layers[[key]]@matrix)) {
                        bp_dir <- assay@layers[[key]]@matrix@dir
                        if (startsWith(bp_dir, old_dir)) {
                            new_bp_dir <- sub(old_dir, new_dir, bp_dir, fixed = TRUE)
                            assay@layers[[key]]@matrix@dir <- new_bp_dir
                        }
                    }
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
