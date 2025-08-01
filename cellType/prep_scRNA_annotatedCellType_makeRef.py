#%%
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245716
import os
import glob
import pandas as pd
import scanpy as sc
from scipy import io
import scipy.sparse
import re
import gzip

# Configuration paths - move these outside functions
DEFAULT_GTF_FILENAME = "Mus_musculus.GRCm39.114.gtf.gz"
GENE_METADATA_NAME = "GSE245716_Metadata.celltype.csv"
GENE_DATA_FOLDER = "GSE245716"
OUTPUT_FILENAME = "refer_adata.h5ad"

def parse_gtf_attributes(attr_string):
    """Parse GTF attributes string into dictionary"""
    attributes = {}
    for item in attr_string.split(';'):
        if item.strip():
            parts = item.strip().split(' ', 1)
            if len(parts) == 2:
                key, value = parts
                attributes[key] = value.strip('"')
    return attributes

def load_gene_mapping_from_gtf(gtf_file):
    """Create gene ID to gene name mapping from GTF file"""
    gene_mapping = {}
    
    with gzip.open(gtf_file, 'rt') if gtf_file.endswith('.gz') else open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) >= 9 and parts[2] == 'gene':
                attributes = parse_gtf_attributes(parts[8])
                gene_id = attributes.get('gene_id')
                gene_name = attributes.get('gene_name')
                
                if gene_id and gene_name:
                    gene_mapping[gene_id] = gene_name
    
    return gene_mapping

def convert_ensembl_to_symbol(adata, gtf_file):
    """
    Replace adata.var_names (Ensembl IDs) with gene symbols using local GTF file.
    Removes duplicate gene symbols by keeping only the first occurrence.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object with Ensembl IDs as var_names
    gtf_file : str
        Path to the GTF file
    
    Returns:
    --------
    adata : AnnData
        Modified AnnData object with gene symbols as var_names (no duplicates)
    """
    print("Loading gene mapping from GTF file...")
    gene_mapping = load_gene_mapping_from_gtf(gtf_file)
    print(f"Loaded {len(gene_mapping)} gene mappings")
    
    ensembl_ids = adata.var_names.tolist()
    
    # Map Ensembl IDs to gene symbols
    mapped_symbols = []
    unmapped_count = 0
    
    for eid in ensembl_ids:
        symbol = gene_mapping.get(eid)
        if symbol:
            mapped_symbols.append(symbol)
        else:
            mapped_symbols.append(eid)  # Keep original ID if no mapping found
            unmapped_count += 1
    
    print(f"Successfully mapped {len(ensembl_ids) - unmapped_count}/{len(ensembl_ids)} genes")
    if unmapped_count > 0:
        print(f"Warning: {unmapped_count} genes could not be mapped and kept their original IDs")
    
    # Handle duplicate gene symbols by keeping only first occurrence
    print("\nChecking for duplicate gene symbols...")
    symbol_counts = pd.Series(mapped_symbols).value_counts()
    duplicated_symbols = symbol_counts[symbol_counts > 1]
    
    if len(duplicated_symbols) > 0:
        print(f"Found {len(duplicated_symbols)} gene symbols with duplicates:")
        print(f"Total duplicate genes to remove: {(duplicated_symbols - 1).sum()}")
        print("Top 10 most duplicated symbols:")
        print(duplicated_symbols.head(10))
        
        # Find indices to keep (first occurrence of each symbol)
        seen_symbols = set()
        indices_to_keep = []
        
        for i, symbol in enumerate(mapped_symbols):
            if symbol not in seen_symbols:
                seen_symbols.add(symbol)
                indices_to_keep.append(i)
        
        print(f"Keeping {len(indices_to_keep)} out of {len(mapped_symbols)} genes")
        
        # Filter adata to keep only first occurrence of each gene symbol
        adata_filtered = adata[:, indices_to_keep].copy()
        adata_filtered.var_names = [mapped_symbols[i] for i in indices_to_keep]
        
        print(f"Final shape: {adata_filtered.shape} (was {adata.shape})")
        return adata_filtered
    
    else:
        print("No duplicate gene symbols found - using original mapping")
        if len(mapped_symbols) == adata.n_vars:
            adata.var_names = mapped_symbols
        else:
            print(f"Warning: Length of mapped_symbols ({len(mapped_symbols)}) does not match adata.n_vars ({adata.n_vars}).")
            print("Gene names (adata.var_names) will not be updated to prevent errors.")
        return adata

#%%
def main(root_path="/scratch/hh01116/download"):
    """
    Main function to process scRNA-seq data and create reference dataset.
    
    Parameters:
    -----------
    root_path : str
        Root directory path where data files are located and output will be saved
    """
    print(f"Using root path: {root_path}")
    
    # Construct full paths
    gtf_path = os.path.join(root_path, DEFAULT_GTF_FILENAME)
    gene_folder = os.path.join(root_path, GENE_DATA_FOLDER)
    metadata_path = os.path.join(gene_folder, GENE_METADATA_NAME)
    output_path = os.path.join(root_path, OUTPUT_FILENAME)
    
    # Verify GTF file exists
    if not os.path.exists(gtf_path):
        raise FileNotFoundError(f"GTF file not found at: {gtf_path}")
    
    # 1) find all sample prefixes by locating barcodes files
    bc_paths = glob.glob(f"{gene_folder}/*_barcodes.tsv.gz", recursive=True)

    if not bc_paths:
        raise FileNotFoundError("No barcodes files found. Make sure you're in the correct directory.")
    
    print(f"Found {len(bc_paths)} sample directories")

    adatas = []
    for bc in bc_paths:
        prefix = bc[:-len("_barcodes.tsv.gz")]
        # original sample string, e.g. "GSM7847300_BS-5307-06162021"
        sample = os.path.basename(prefix)
        # take only the part after the first underscore: "BS-5307-06162021"
        partial = re.sub(r'^.*?_', '', sample)

        mtx = prefix + "_matrix.mtx.gz"
        feat = prefix + "_features.tsv.gz"

        print(f"Processing sample: {partial}")

        # 2) read raw data
        mat = io.mmread(mtx).tocsr().T
        barcodes = pd.read_csv(bc, header=None)[0].astype(str)
        # prefix each barcode with the partial sample ID
        barcodes = partial + "_" + barcodes

        features = pd.read_csv(feat, header=None, sep="\t")[0].astype(str)

        # 3) build AnnData - explicitly ensure X is sparse
        ad = sc.AnnData(X=mat)  # mat is already sparse CSR
        ad.obs_names = barcodes
        ad.var_names = features
        
        # Verify matrix is sparse
        if not scipy.sparse.issparse(ad.X):
            print(f"Warning: Converting dense matrix to sparse for sample {partial}")
            ad.X = scipy.sparse.csr_matrix(ad.X)
        
        adatas.append(ad)

    # 4) concatenate
    print("Concatenating all samples...")
    adata = sc.concat(adatas, join="outer", label="sample", index_unique=None)
    print(f"Combined data shape: {adata.shape}")
    
    # Ensure concatenated matrix is sparse
    if not scipy.sparse.issparse(adata.X):
        print("Converting concatenated matrix to sparse format...")
        adata.X = scipy.sparse.csr_matrix(adata.X)
    
    print(f"Matrix format: {type(adata.X)}")
    print(f"Matrix sparsity: {1 - adata.X.nnz / (adata.X.shape[0] * adata.X.shape[1]):.4f}")
    
    # 5) read metadata and join all available annotations
    print(f"Loading metadata from: {metadata_path}")
    if not os.path.exists(metadata_path):
        raise FileNotFoundError(f"Metadata file not found at: {metadata_path}")
        
    meta = pd.read_csv(metadata_path, index_col=0)
    print(f"Metadata shape: {meta.shape}")
    
    adata.obs = adata.obs.join(meta, how="left")

    # 5b) filter to retain only spots present in the CSV
    adata = adata[adata.obs.index.isin(meta.index)].copy()
    print(f"After filtering to metadata: {adata.shape}")
    
    return adata, gtf_path, output_path


if __name__ == "__main__":
    # You can change the root path here or pass it as a command line argument
    ROOT_PATH = "/scratch/hh01116/download"
    
    # Uncomment below if you want to use command line arguments
    # import sys
    # if len(sys.argv) > 1:
    #     ROOT_PATH = sys.argv[1]
    
    adata, gtf_file_path, output_file_path = main(root_path=ROOT_PATH)
    
    # Convert Ensembl IDs to gene symbols (with duplicate handling)
    print(f"\nConverting Ensembl IDs to gene symbols using: {gtf_file_path}")
    adata = convert_ensembl_to_symbol(adata, gtf_file=gtf_file_path)
    
    # Additional validation: ensure no duplicate gene names
    print("\nFinal validation:")
    print(f"Total genes: {adata.n_vars}")
    print(f"Unique gene names: {len(adata.var_names.unique())}")
    
    duplicate_check = adata.var_names.duplicated()
    if duplicate_check.any():
        print(f"ERROR: Still found {duplicate_check.sum()} duplicate gene names!")
        duplicated_genes = adata.var_names[duplicate_check]
        print(f"Duplicated genes: {duplicated_genes.tolist()}")
        raise ValueError("Duplicate gene names still present after processing")
    else:
        print("✓ No duplicate gene names found")
    
    print(f"\nSaving processed data to: {output_file_path}")
    adata.write_h5ad(output_file_path)
    print("Processing completed successfully!")