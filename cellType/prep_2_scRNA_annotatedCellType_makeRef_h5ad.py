import scanpy as sc
import pandas as pd
import numpy as np
import yaml
import os

config = yaml.load(open("config/batch_config.yaml"), Loader=yaml.FullLoader)
project_dir = config["project_dir"]

print("Reading exported data...")

# Read the exported data
counts = pd.read_csv(f"{project_dir}rds_data/GBM_counts.csv", index_col=0)
metadata = pd.read_csv(f"{project_dir}rds_data/GBM_metadata.csv", index_col=0)
genes = pd.read_csv(f"{project_dir}rds_data/GBM_genes.csv", index_col=0)

print(f"Counts matrix shape: {counts.shape}")
print(f"Metadata shape: {metadata.shape}")
print(f"Genes shape: {genes.shape}")

# Create AnnData object (transpose counts for cells x genes format)
print("Creating AnnData object...")
adata = sc.AnnData(X=counts.T.values.astype(np.float32))

# Add metadata
adata.obs = metadata
adata.var = genes
adata.var_names = counts.index
adata.obs_names = counts.columns

# Verify required columns exist
required_cols = ["sample", "group", "replicate", "cell", "cell_type"]
for col in required_cols:
    if col in adata.obs.columns:
        print(f"✓ Found required column: {col}")
        print(f"  Unique values: {adata.obs[col].unique()}")
    else:
        print(f"✗ Missing required column: {col}")

# Save the h5ad file
print("Saving h5ad file...")
adata.write(f"{project_dir}rds_data/GBM_earlyStage_scSeq_annotated_cell2location.h5ad")
print("Successfully created GBM_earlyStage_scSeq_annotated_cell2location.h5ad")

print("removing GBM_counts.csv, GBM_metadata.csv, GBM_genes.csv")
os.remove(f"{project_dir}rds_data/GBM_counts.csv")
os.remove(f"{project_dir}rds_data/GBM_metadata.csv")
os.remove(f"{project_dir}rds_data/GBM_genes.csv")

# Print summary
print(f"\nSummary:")
print(f"- Cells: {adata.n_obs}")
print(f"- Genes: {adata.n_vars}")
print(f"- Cell types: {len(adata.obs['cell_type'].unique())}")
print(f"- Samples: {len(adata.obs['sample'].unique())}")
