# data `GSE195848_Seurat_object.RDS` obtained from
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE195848
source("util_headers.r")

library(Seurat)
scRNA <- readRDS("GSE195848_Seurat_object.RDS")
# only keep the early stage cells
# early-stage (10 days) and late-stage (>30 days)
scRNA_early <- subset(scRNA, stage == "Early")


head(scRNA_early@meta.data)

cols <- c("batch", "sample", "group", "replicate", "cell", "cell_type")
for (col in cols) {
    cat("Values in", col, ":\n")
    print(unique(scRNA_early@meta.data[[col]]))
    cat("\n")
}

# ===================================
#  spacexr
# ===================================
# Create a reference object from scRNA-seq data for spacexr
cell_type <- scRNA_early@meta.data$cell_type
sc_count <- scRNA_early@assays$RNA@counts
names(cell_type) <- colnames(sc_count)
nUMI <- colSums(sc_count)
names(nUMI) <- colnames(sc_count)

reference <- spacexr::Reference(sc_count, cell_type, nUMI)
saveRDS(reference, paste0(project_dir, "rds_data/", "ref_GBM_scSeq.rds"))


# ===================================
# cell2location
# ===================================
# Create a reference object for cell2location

# BiocManager::install("zellkonverter", dependencies = TRUE)
sce_obj <- Seurat::as.SingleCellExperiment(scRNA_early)
zellkonverter::writeH5AD(sce_obj, paste0(project_dir, "rds_data/", "ref_GBM_scSeq.h5ad"))


# # Manual export approach if the above doesn't work
# counts_matrix <- GetAssayData(scRNA_early, layer = "counts", assay = "RNA")
# metadata <- scRNA_early@meta.data
# # Export to files
# write.csv(as.matrix(counts_matrix), paste0(project_dir, "rds_data/", "GBM_counts.csv"))
# write.csv(metadata, paste0(project_dir, "rds_data/", "GBM_metadata.csv"), row.names = TRUE)
# # Create gene information
# genes_df <- data.frame(
#     gene_ids = rownames(counts_matrix),
#     gene_symbols = rownames(counts_matrix),
#     row.names = rownames(counts_matrix)
# )
# write.csv(genes_df, paste0(project_dir, "rds_data/", "GBM_genes.csv"), row.names = TRUE)
# cat("Files exported successfully:\n")
# cat("- GBM_counts.csv\n")
# cat("- GBM_metadata.csv\n")
# cat("- GBM_genes.csv\n")
# cat("Now run the Python script to create h5ad file.\n")
