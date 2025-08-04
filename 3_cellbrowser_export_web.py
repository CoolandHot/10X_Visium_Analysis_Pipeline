#!/usr/bin/env python3
"""
Scanpy-based Cell Browser exporter for spatial transcriptomics data.
Reads h5ad files and exports to UCSC Cell Browser format with cluster data.
https://cellbrowser.readthedocs.io/en/master/basic_usage.html

cd <path_to_the_exported_folder>
# build & start the server. But this script has included the build step.
cbBuild -o $VOL_NOBACKUP/10X_DIPG/output/html_reports/cellbrowser -p 8899


<at any folder>
# start the server withouth building
cbUpgrade -o $VOL_NOBACKUP/10X_DIPG/output/html_reports/cellbrowser -p 8899
"""

import pandas as pd
from pathlib import Path
from typing import Dict, List, Any, Optional
import yaml
import numpy as np
import copy
import subprocess
import shutil
import re

# Core scverse ecosystem
import scanpy as sc
import anndata as ad


class CellBrowserExporter:
    """
    Scanpy-based exporter for UCSC Cell Browser from h5ad files.
    Creates all necessary files for Cell Browser including cluster data from CSV files.
    """

    def __init__(self, config: Dict[str, Any]):
        """
        Initialize Cell Browser exporter.

        Parameters:
        -----------
        config : Dict[str, Any]
            Configuration dictionary with analysis parameters
        """
        self.config = config
        self.cellbrowser_config = config.get("cellbrowser", {})

        # Create Cell Browser output directory
        self.cellbrowser_dir = Path(self.cellbrowser_config.get(
            "cellbrowser_dir", "output/cellbrowser_geneExpr"
        ))
        self.cellbrowser_dir.mkdir(exist_ok=True)

        # Cluster file paths from config
        self.cluster_files = {}
        if "cluster_files" in self.cellbrowser_config:
            for cluster_name, file_path in self.cellbrowser_config[
                "cluster_files"
            ].items():
                self.cluster_files[cluster_name] = Path(file_path)

    def load_h5ad_data(self, h5ad_path: Path) -> ad.AnnData:
        """Load h5ad file using scanpy"""
        print(f"Loading h5ad data from: {h5ad_path}")
        adata = sc.read_h5ad(h5ad_path)
        print(f"Loaded data with {adata.n_obs} cells and {adata.n_vars} genes")
        return adata

    def load_cluster_data(self) -> Dict[str, pd.Series]:
        """Load cluster data from CSV files and handle missing cells"""
        cluster_data = {}

        cluster_method = self.config.get("cluster_method", None)
        dga = self.config.get("Differential_Gene_Analysis", {})
        outTumour = set(dga.get("outTumour_cluster_nums", []))
        inTumour = set(dga.get("inTumour_cluster_nums", []))
        edgeTumour = set(dga.get("edgeTumour_cluster_nums", []))

        for cluster_name, file_path in self.cluster_files.items():
            if file_path.exists():
                print(f"Loading cluster data from: {file_path}")
                df = pd.read_csv(file_path, header=0)

                # Assume first column is cell names, second is cluster assignments
                df.columns = ["cell_id", cluster_name]

                # If filename includes cluster_method, apply mapping logic
                if cluster_method and cluster_method in str(file_path):
                    def map_cluster(x):
                        try:
                            val = int(x)
                        except Exception:
                            return str(x)
                        if val in outTumour:
                            return "outTumour"
                        elif val in inTumour:
                            return "inTumour"
                        elif val in edgeTumour:
                            return "edgeTumour"
                        else:
                            return str(x)
                    df[cluster_name] = df[cluster_name].apply(map_cluster)

                cluster_series = df.set_index("cell_id")[cluster_name]
                cluster_data[cluster_name] = cluster_series
                print(
                    f"Loaded {len(cluster_series)} cluster assignments for {cluster_name}"
                )
            else:
                print(f"Warning: Cluster file not found: {file_path}")

        return cluster_data

    def merge_cluster_data(
        self, adata: ad.AnnData, cluster_data: Dict[str, pd.Series]
    ) -> ad.AnnData:
        """Merge cluster data into adata.obs, adding NA for missing cells"""
        adata_copy = adata.copy()

        for cluster_name, cluster_series in cluster_data.items():
            print(f"Merging {cluster_name} data...")

            # Create a series with all cell IDs from adata, initially as object dtype
            merged_cluster = pd.Series(
                index=adata_copy.obs_names,
                dtype="object",  # Changed from 'category' to 'object'
                name=cluster_name,
            )

            # Fill missing cells with 'NA' first
            merged_cluster[:] = "NA"

            # Fill in available cluster assignments
            available_cells = merged_cluster.index.intersection(cluster_series.index)
            merged_cluster.loc[available_cells] = cluster_series.loc[
                available_cells
            ].astype(str)

            # Convert to categorical after all values are set
            merged_cluster = merged_cluster.astype("category")

            # Add to adata.obs
            adata_copy.obs[cluster_name] = merged_cluster

            print(
                f"Added {cluster_name} to metadata: {len(available_cells)} matched, "
                f"{len(merged_cluster) - len(available_cells)} set to NA"
            )

        return adata_copy

    def _shift_spatial_coordinates(self, adata: ad.AnnData) -> ad.AnnData:
        """
        Shift spatial coordinates to create a tiled layout for multiple batches.
        Similar to the R function form_square_tile_coor.
        """
        if "spatial" not in adata.obsm:
            print("No spatial coordinates found, skipping coordinate shifting")
            return adata

        adata_copy = adata.copy()
        batch_names = adata_copy.obs["batch"].unique()
        n_batches = len(batch_names)

        # Extract coordinates for each batch
        batch_coords = {}
        for batch_name in batch_names:
            batch_mask = adata_copy.obs["batch"] == batch_name
            batch_coords[batch_name] = adata_copy.obsm["spatial"][batch_mask]

        # Calculate shift values based on all batches
        all_coords = np.vstack(list(batch_coords.values()))
        max_x_shift = all_coords[:, 0].max()
        max_y_shift = all_coords[:, 1].max()

        # Define shifts based on number of batches
        if 2 <= n_batches <= 8:
            # Determine grid size (cols, rows) for up to 8 batches
            grid_sizes = {
                2: (2, 1),
                3: (2, 2),
                4: (2, 2),
                5: (3, 2),
                6: (3, 2),
                7: (4, 2),
                8: (4, 2),
            }
            n_cols, n_rows = grid_sizes[n_batches]
            shifts = {}
            for i, batch_name in enumerate(batch_names):
                row = i // n_cols
                col = i % n_cols
                shifts[batch_name] = [col * max_x_shift * 1.05, row * max_y_shift * -1.05]
        else:
            raise ValueError(f"Only 2-8 batches are supported, got {n_batches}")

        # Apply shifts to coordinates
        shifted_coords = adata_copy.obsm["spatial"].copy()
        for batch_name in batch_names:
            batch_mask = adata_copy.obs["batch"] == batch_name
            shift = shifts[batch_name]
            shifted_coords[batch_mask, 0] += shift[0]  # x shift
            shifted_coords[batch_mask, 1] += shift[1]  # y shift

        adata_copy.obsm["spatial"] = shifted_coords
        print(f"Applied coordinate shifts for {n_batches} batches")

        return adata_copy

    def export_for_cellbrowser(self, h5ad_path: Optional[Path] = None):
        """
        Export h5ad data for UCSC Cell Browser.

        Parameters:
        -----------
        h5ad_path : Optional[Path]
            Path to h5ad file. If None, will look for default locations.
        """
        # Find h5ad file if not provided
        if h5ad_path is None:
            h5ad_path = self._find_h5ad_file()

        # Load data
        adata = self.load_h5ad_data(h5ad_path)

        # Load and merge cluster data
        cluster_data = self.load_cluster_data()
        if cluster_data:
            adata = self.merge_cluster_data(adata, cluster_data)

        # Get configuration from cellbrowser section
        color_field = self.cellbrowser_config.get("color_field", "leiden")
        label_field = self.cellbrowser_config.get("label_field", "leiden")
        enum_fields = self.cellbrowser_config.get("enum_fields", ["leiden", "batch"])

        # Add cluster fields to enum_fields
        cluster_fields = list(cluster_data.keys()) if cluster_data else []
        enum_fields.extend(cluster_fields)
        enum_fields = list(set(enum_fields))  # Remove duplicates

        # Process data based on batch presence
        if "batch" in adata.obs.columns:
            print("Found batch column - combining all batches with shifted coordinates")
            # Shift spatial coordinates for tiled layout
            adata = self._shift_spatial_coordinates(adata)
            # Export all data as single combined batch
            self._export_batch(adata, color_field, label_field, enum_fields)
        else:
            # Export all data as single batch
            print("No batch column found - exporting all data")
            self._export_batch(adata, color_field, label_field, enum_fields)

    def _find_h5ad_file(self) -> Path:
        """Find h5ad file in expected locations"""
        # Look for common h5ad file patterns
        possible_paths = [
            Path(self.config['rds_data_dir'])
            / f"{self.config['output_file_prefix']}_merged.h5ad",
            Path(self.config['rds_data_dir'])
            / f"{self.config['output_file_prefix']}_banksy_merged.h5ad",
            Path(self.config['rds_data_dir'])
            / f"{self.config['output_file_prefix']}_BayesSpace_merged.h5ad",
        ]

        for path in possible_paths:
            if path.exists():
                return path

        # Look for any h5ad files in output directory
        h5ad_files = list(Path(self.config['rds_data_dir']).glob("*.h5ad"))
        if h5ad_files:
            return h5ad_files[0]

        raise FileNotFoundError("No h5ad file found in expected locations")

    def _compute_marker_genes(self, adata: ad.AnnData, enum_fields: List[str]):
        """
        Compute marker genes for each cluster key if not already computed.

        Parameters:
        -----------
        adata : ad.AnnData
            Annotated data object
        enum_fields : List[str]
            List of cluster fields to compute markers for
        """
        cluster_fields = [
            field
            for field in enum_fields
            if field != "batch" and field in adata.obs.columns
        ]

        for cluster_key in cluster_fields:
            marker_key = f"rank_genes_groups_{cluster_key}"
            marker_file = self.cellbrowser_dir / f"markers_{cluster_key}.tsv"

            # Check if marker file already exists
            if marker_file.exists():
                print(
                    f"Marker file already exists for {cluster_key}, skipping computation"
                )
                continue

            # Check if markers already computed
            if marker_key not in adata.uns:
                print(f"Computing marker genes for {cluster_key}...")
                try:
                    # Exclude 'NA' clusters from marker computation
                    adata_subset = adata[adata.obs[cluster_key] != "NA"].copy()
                    if adata_subset.n_obs > 0:
                        sc.tl.rank_genes_groups(
                            adata_subset, cluster_key, method="wilcoxon", key_added=marker_key
                        )
                        # Copy computed markers back to the original adata object
                        adata.uns[marker_key] = adata_subset.uns[marker_key]
                        print(f"Successfully computed markers for {cluster_key}")
                    else:
                        print(f"No cells with valid clusters for {cluster_key}, skipping marker computation.")

                except Exception as e:
                    print(f"Error computing markers for {cluster_key}: {e}")
            else:
                print(f"Markers already computed for {cluster_key}")

    def _export_batch(
        self,
        adata: ad.AnnData,
        color_field: str,
        label_field: str,
        enum_fields: List[str],
    ):
        """Export data for a specific batch"""
        # Create batch-specific directory
        batch_dir = self.cellbrowser_dir

        # Use all data (no filtering since we're combining batches)
        batch_adata = adata.copy()
        print(f"Exporting data for batch: {self.config['output_file_prefix']}")

        try:
            # Remove duplicated gene names before exporting
            if batch_adata.var_names.str.lower().duplicated().any():
                print("Removing duplicated gene names (case-insensitive)...")
                batch_adata = batch_adata[:, ~batch_adata.var_names.str.lower().duplicated(keep='first')].copy()

            # Compute marker genes for all cluster fields
            self._compute_marker_genes(batch_adata, enum_fields)

            # Export marker genes for each cluster field
            marker_files = self._export_markers_scanpy(
                batch_adata, batch_dir, enum_fields
            )

            # Export expression matrix using scanpy
            self._export_expression_matrix(batch_adata, batch_dir)

            # Export barcodes
            self._export_barcodes(batch_adata, batch_dir)

            # Export features/genes
            self._export_features(batch_adata, batch_dir)

            # Export metadata
            self._export_metadata(batch_adata, batch_dir, enum_fields)

            # Export coordinates
            coord_bounds = self._export_coordinates(batch_adata, batch_dir)

            # Export quick genes
            self._export_quick_genes(batch_adata, batch_dir)

            # Create cellbrowser.conf
            self._create_cellbrowser_config(
                batch_dir,
                color_field,
                label_field,
                enum_fields,
                coord_bounds,
                marker_files
            )

            # Create desc.conf file
            self._create_desc_config(batch_dir, "gene expression")

            print("Successfully exported Cell Browser files")

        except Exception as e:
            print(f"Error exporting Cell Browser files: {e}")
            raise

    def _export_expression_matrix(self, adata: ad.AnnData, output_dir: Path):
        """Export expression matrix using scanpy functionality"""
        matrix_file_gz = output_dir / "matrix.mtx.gz"

        # Check if matrix file already exists
        if matrix_file_gz.exists():
            print("Expression matrix file already exists, skipping export")
            return

        from scipy.io import mmwrite
        import gzip

        # Use raw data if available, otherwise use processed data
        if adata.raw is not None:
            X = adata.raw.X
        else:
            X = adata.X

        # Fill NA with 0 before exporting
        if hasattr(X, "toarray"):
            # Dense numpy array or matrix
            X = np.nan_to_num(X, nan=0)
            matrix = X.T  # Cell Browser expects genes x cells
        else:
            from scipy import sparse
            # Convert to sparse, then fillna if possible
            X = sparse.csr_matrix(X)
            # Fill NaN with 0 in sparse matrix
            if np.isnan(X.data).any():
                X.data[np.isnan(X.data)] = 0
            matrix = X.T

        # Write matrix
        matrix_file = output_dir / "matrix.mtx"
        mmwrite(matrix_file, matrix)

        # Compress
        with open(matrix_file, "rb") as f_in:
            with gzip.open(matrix_file_gz, "wb") as f_out:
                f_out.writelines(f_in)

        matrix_file.unlink()

    def _export_barcodes(self, adata: ad.AnnData, output_dir: Path):
        """Export cell barcodes"""
        barcodes_file = output_dir / "barcodes.tsv.gz"

        # Check if barcodes file already exists
        if barcodes_file.exists():
            print("Barcodes file already exists, skipping export")
            return

        import gzip

        with gzip.open(barcodes_file, "wt") as f:
            for barcode in adata.obs_names:
                f.write(f"{barcode}\n")

    def _export_features(self, adata: ad.AnnData, output_dir: Path):
        """Export gene features"""
        features_file = output_dir / "features.tsv.gz"

        # Check if features file already exists
        if features_file.exists():
            print("Features file already exists, skipping export")
            return

        import gzip

        if adata.raw is not None:
            gene_names = adata.raw.var_names
            gene_ids = (
                adata.raw.var.index if hasattr(adata.raw.var, "index") else gene_names
            )
        else:
            gene_names = adata.var_names
            gene_ids = adata.var.index if hasattr(adata.var, "index") else gene_names

        with gzip.open(features_file, "wt") as f:
            for gene_id, gene_name in zip(gene_ids, gene_names):
                f.write(f"{gene_id}\t{gene_name}\tGene Expression\n")

    def _export_metadata(
        self, adata: ad.AnnData, output_dir: Path, enum_fields: List[str]
    ):
        """Export cell metadata"""
        meta_file = output_dir / "meta.tsv"

        # Check if metadata file already exists
        if meta_file.exists():
            print("Metadata file already exists, skipping export")
            return

        meta_df = adata.obs.copy()
        meta_df.insert(0, "cell_id", adata.obs_names)

        # Drop 'library_id' column as it's same as 'batch'
        if "library_id" in meta_df.columns:
            meta_df = meta_df.drop(columns=["library_id"])

        # Ensure enum fields are strings and handle categorical data
        for field in enum_fields:
            if field in meta_df.columns:
                if isinstance(meta_df[field].dtype, pd.CategoricalDtype):
                    meta_df[field] = meta_df[field].astype(str)
                else:
                    meta_df[field] = meta_df[field].astype(str)

        meta_df.to_csv(meta_file, sep="\t", index=False)

    def _export_coordinates(
        self, adata: ad.AnnData, output_dir: Path
    ) -> Dict[str, Dict[str, float]]:
        """Export spatial and UMAP coordinates"""
        spatial_file = output_dir / "spatial_coords.tsv"
        umap_file = output_dir / "umap_coords.tsv"
        coord_bounds = {}

        # Export spatial coordinates
        if "spatial" in adata.obsm:
            spatial_coords = adata.obsm["spatial"]

            # Check if spatial coordinates file already exists
            if not spatial_file.exists():
                spatial_df = pd.DataFrame(
                    {
                        "cell_id": adata.obs_names,
                        "x": spatial_coords[:, 0],
                        "y": spatial_coords[:, 1],
                    }
                )
                spatial_df.to_csv(spatial_file, sep="\t", index=False)
            else:
                print("Spatial coordinates file already exists, skipping export")

            coord_bounds["spatial"] = {
                "minX": float(spatial_coords[:, 0].min()),
                "maxX": float(spatial_coords[:, 0].max()),
                "minY": float(spatial_coords[:, 1].min()),
                "maxY": float(spatial_coords[:, 1].max()),
            }

        # Export UMAP coordinates
        if "X_umap" in adata.obsm:
            umap_coords = adata.obsm["X_umap"]

            # Check if UMAP coordinates file already exists
            if not umap_file.exists():
                umap_df = pd.DataFrame(
                    {
                        "cell_id": adata.obs_names,
                        "x": umap_coords[:, 0],
                        "y": umap_coords[:, 1],
                    }
                )
                umap_df.to_csv(umap_file, sep="\t", index=False)
            else:
                print("UMAP coordinates file already exists, skipping export")

            coord_bounds["umap"] = {
                "minX": float(umap_coords[:, 0].min()),
                "maxX": float(umap_coords[:, 0].max()),
                "minY": float(umap_coords[:, 1].min()),
                "maxY": float(umap_coords[:, 1].max()),
            }

        return coord_bounds

    def _export_markers_scanpy(
        self, adata: ad.AnnData, output_dir: Path, enum_fields: List[str]
    ) -> List[Dict[str, str]]:
        """
        Export marker genes for each cluster key using scanpy

        Parameters:
        -----------
        adata : ad.AnnData
            Annotated data object
        output_dir : Path
            Output directory for marker files
        enum_fields : List[str]
            List of cluster fields to export markers for

        Returns:
        --------
        List[Dict[str, str]]
            List of marker file configurations for cellbrowser
        """
        marker_files = []
        cluster_fields = [
            field
            for field in enum_fields
            if field != "batch" and field in adata.obs.columns
        ]

        for cluster_key in cluster_fields:
            marker_key = f"rank_genes_groups_{cluster_key}"
            markers_file = output_dir / f"markers_{cluster_key}.tsv"

            # Check if marker file already exists
            if markers_file.exists():
                print(f"Marker file already exists for {cluster_key}, skipping export")
                # Still add to marker files list for config
                marker_files.append(
                    {
                        "file": f"markers_{cluster_key}.tsv",
                        "shortLabel": f"{cluster_key} Markers",
                    }
                )
                continue

            try:
                # Check if marker data exists
                if marker_key not in adata.uns:
                    print(f"No marker data found for {cluster_key}, skipping...")
                    continue

                print(f"Exporting marker genes for {cluster_key}...")

                # Get marker genes dataframe
                markers_df = sc.get.rank_genes_groups_df(
                    adata, group=None, key=marker_key
                )

                # Reformat for Cell Browser
                cb_markers = []
                
                # Ensure we only use markers for clusters present in the metadata
                valid_clusters = set(adata.obs[cluster_key].unique())

                for _, row in markers_df.iterrows():
                    cluster_name = str(row["group"])
                    if cluster_name not in valid_clusters:
                        continue

                    cb_markers.append(
                        {
                            "cluster": cluster_name,
                            "gene": row["names"],
                            "pval": row["pvals_adj"]
                            if "pvals_adj" in row
                            else row.get("pvals", 1.0),
                            "logfc": row["logfoldchanges"]
                            if "logfoldchanges" in row
                            else 0.0,
                        }
                    )

                markers_cb_df = pd.DataFrame(cb_markers)
                markers_cb_df.to_csv(markers_file, sep="\t", index=False)

                # Add to marker files list
                marker_files.append(
                    {
                        "file": f"markers_{cluster_key}.tsv",
                        "shortLabel": f"{cluster_key} Markers",
                    }
                )

                print(f"Exported markers for {cluster_key} to {markers_file}")

            except Exception as e:
                print(f"Could not export markers for {cluster_key}: {e}")
                # Create empty file
                pd.DataFrame(columns=["cluster", "gene", "pval", "logfc"]).to_csv(
                    markers_file, sep="\t", index=False
                )

        return marker_files

    def _export_quick_genes(self, adata: ad.AnnData, output_dir: Path):
        """Export quick genes list"""
        quick_genes_file = output_dir / "quickGenes.tsv"

        # Get all available gene names
        if adata.raw is not None:
            all_genes = adata.raw.var_names.tolist()
        else:
            all_genes = adata.var_names.tolist()

        # Default quick genes from config or use defaults
        default_genes = ["Arg1", "IDO1", "Cd45", "Cd3e", "Cd8a", "Cd4", "Foxp3", "Gfp", "Rfp", "Egfr"]
        quick_genes = self.cellbrowser_config.get("quick_genes", default_genes)

        def normalize_gene(gene):
            return gene.lower().replace("-", "").replace(".", "").replace("_", "").replace(" ", "").strip()

        all_genes_norm = [normalize_gene(gene) for gene in all_genes]
        quick_genes = [normalize_gene(gene) for gene in quick_genes]

        # Filter to genes present in data
        available_genes_indices = [all_genes_norm.index(gene) for gene in quick_genes if gene in all_genes_norm]
        available_genes = [all_genes[i] for i in available_genes_indices]

        with open(quick_genes_file, "w") as f:
            for i, gene in enumerate(available_genes):
                # Write as gene|gene to suppress Cell Browser warning
                f.write(f"{gene}|{gene}\t{i}\n")

        print(f"Exported {len(available_genes)} quick genes")

    def _create_desc_config(self, output_dir: Path, data_type: str = "gene expression"):
        """Create desc.conf configuration file"""
        desc_file = output_dir / "desc.conf"

        if data_type == "gene expression":
            title = "Jake's GBM project - data on gene expression"
            abstract = """This dataset builds upon the gene expression data for four groups of GBM samples sequenced using 10X Genomics Visium technology. The groups are ADI, RAD, COMB, and SAL."""
        else:  # cell type population density
            title = "Jake's GBM project - data on cell type population density"
            abstract = """This dataset builds upon the cell type population density data for four groups of GBM samples sequenced using 10X Genomics Visium technology. The groups are ADI, RAD, COMB, and SAL."""

        desc_content = f"""# https://github.com/maximilianh/cellBrowser/blob/master/src/cbPyLib/cellbrowser/sampleConfig/desc.conf
title = "{title}"

abstract = \"\"\"
{abstract}
\"\"\"

methods = \"\"\"
Four groups of GBM samples were sequenced using 10X Genomics Visium technology. ADI, RAD, COMB and SAL. 
\"\"\"
"""

        with open(desc_file, "w") as f:
            f.write(desc_content)

        print(f"Created desc.conf file: {desc_file}")

    def _create_cellbrowser_config(
        self,
        output_dir: Path,
        color_field: str,
        label_field: str,
        enum_fields: List[str],
        coord_bounds: Dict,
        marker_files: List[Dict[str, str]],
        project_suffix: str = "geneExpr",
    ):
        """Create cellbrowser.conf configuration file"""
        config_file = output_dir / "cellbrowser.conf"
        project_name = f"{self.config['output_file_prefix']}_{project_suffix}"

        default_bounds = {"minX": 0, "minY": 0, "maxX": 2000, "maxY": 2000}
        spatial_bounds = coord_bounds.get("spatial", default_bounds)
        umap_bounds = coord_bounds.get("umap", None)

        coords_block = f"""[
    {{
        "file": "spatial_coords.tsv",
        "shortLabel": "spatial",
        "flipY": 1,
        "useRaw": 1,
        "minX": {spatial_bounds["minX"]:.0f},
        "minY": {spatial_bounds["minY"]:.0f},
        "maxX": {spatial_bounds["maxX"]:.0f},
        "maxY": {spatial_bounds["maxY"]:.0f}
    }}"""
        if umap_bounds is not None:
            coords_block += f""",
    {{
        "file": "umap_coords.tsv",
        "shortLabel": "umap",
        "flipY": 1,
        "useRaw": 1,
        "minX": {umap_bounds["minX"]:.0f},
        "minY": {umap_bounds["minY"]:.0f},
        "maxX": {umap_bounds["maxX"]:.0f},
        "maxY": {umap_bounds["maxY"]:.0f}
    }}"""
        coords_block += "]"

        # Format marker files list
        markers_block = "["
        for i, marker_file in enumerate(marker_files):
            if i > 0:
                markers_block += ", "
            markers_block += f'{{"file": "{marker_file["file"]}", "shortLabel": "{marker_file["shortLabel"]}"}}'
        markers_block += "]"

        config_content = f"""
# https://github.com/maximilianh/cellBrowser/blob/master/src/cbPyLib/cellbrowser/sampleConfig/cellbrowser.conf

name='{project_name}'
shortLabel='{project_name}'
exprMatrix='matrix.mtx.gz'
tags = ["10x", 'VisiumHD']
meta='meta.tsv'
geneIdType='auto'
defCatPal = 'tol-rainbow'
defColorField='{color_field}'
labelField='{label_field}'
enumFields={enum_fields}
coords={coords_block}
markers = {markers_block}
quickGenesFile="quickGenes.tsv"
alpha = 1
radius = 3
"""

        with open(config_file, "w") as f:
            f.write(config_content)


class Cell2LocationCellBrowserExporter(CellBrowserExporter):
    """
    Cell Browser exporter for cell2location cell type abundance data.
    Creates AnnData object from CSV abundance data and exports to Cell Browser format.
    """

    def __init__(self, config: Dict[str, Any]):
        super().__init__(config)

        # Get cell2location specific config
        self.cellbrowser_dir = Path(self.cellbrowser_config.get(
            "cell_abundances_dir", "output/cellbrowser_cellTypes"
        ))
        self.cellbrowser_dir.mkdir(exist_ok=True)
        self.abundances_csv_path = self.cellbrowser_config.get(
            "cell_abundances_csv",
            "output/cell_type_prediction/cell2location_results/cell_abundances_and_clusters.csv",
        )
        self.abundances_csv_path = Path(self.abundances_csv_path)

    def load_abundances_data(self) -> ad.AnnData:
        """Load cell type abundances from CSV and create AnnData object"""
        print(f"Loading cell type abundances from: {self.abundances_csv_path}")

        # Load CSV data
        df = pd.read_csv(self.abundances_csv_path, header=0, index_col=0)

        # Separate metadata columns and abundance data
        metadata_cols = ["region_cluster", "batch"]
        obs_data = {}

        # Extract metadata columns if they exist
        for col in metadata_cols:
            if col in df.columns:
                obs_data[col] = df[col]

        obs_df = pd.DataFrame(obs_data, index=df.index)

        # Get abundance data (exclude metadata columns)
        abundance_cols = [col for col in df.columns if col not in metadata_cols]
        X_df = df[abundance_cols].copy()

        # Make cell type names legal by replacing spaces with underscores
        X_df.columns = [col.replace(" ", "_") for col in X_df.columns]

        # Create AnnData object
        adata = ad.AnnData(
            X=X_df.values, obs=obs_df, var=pd.DataFrame(index=X_df.columns)
        )
        adata.obs_names = df.index
        adata.var_names = X_df.columns

        print(f"Created AnnData with {adata.n_obs} cells and {adata.n_vars} cell types")
        return adata

    def load_spatial_coordinates_from_h5ad(
        self, adata_abundances: ad.AnnData
    ) -> ad.AnnData:
        """Load spatial coordinates from the original h5ad file and match to abundance data"""
        try:
            # Find and load the original h5ad file to get spatial coordinates
            h5ad_path = self._find_h5ad_file()
            adata_spatial = sc.read_h5ad(h5ad_path)

            if "spatial" in adata_spatial.obsm:
                # Match spatial coordinates to abundance data cells by cell names
                common_cells = adata_abundances.obs_names.intersection(
                    adata_spatial.obs_names
                )

                if len(common_cells) > 0:
                    # Reorder spatial coordinates to match abundance data
                    spatial_coords_matched = np.zeros((len(adata_abundances), 2))

                    for i, cell_name in enumerate(adata_abundances.obs_names):
                        if cell_name in adata_spatial.obs_names:
                            idx = adata_spatial.obs_names.get_loc(cell_name)
                            spatial_coords_matched[i] = adata_spatial.obsm["spatial"][
                                idx
                            ]

                    adata_abundances.obsm["spatial"] = spatial_coords_matched
                    print(f"Matched spatial coordinates for {len(common_cells)} cells")

                    # Also copy UMAP coordinates if available
                    if "X_umap" in adata_spatial.obsm:
                        umap_coords_matched = np.zeros((len(adata_abundances), 2))
                        for i, cell_name in enumerate(adata_abundances.obs_names):
                            if cell_name in adata_spatial.obs_names:
                                idx = adata_spatial.obs_names.get_loc(cell_name)
                                umap_coords_matched[i] = adata_spatial.obsm["X_umap"][
                                    idx
                                ]
                        adata_abundances.obsm["X_umap"] = umap_coords_matched
                        print("Also matched UMAP coordinates")
                else:
                    print(
                        "Warning: No common cells found between abundance and spatial data"
                    )
            else:
                print("No spatial coordinates found in original h5ad file")

        except Exception as e:
            print(f"Could not load spatial coordinates: {e}")

        return adata_abundances

    def export_for_cellbrowser(self):
        """
        Export cell type abundance data for UCSC Cell Browser.
        """
        # Load abundance data instead of h5ad
        adata = self.load_abundances_data()

        # Load spatial coordinates from original h5ad file
        adata = self.load_spatial_coordinates_from_h5ad(adata)

        # Load and merge cluster data
        cluster_data = self.load_cluster_data()
        if cluster_data:
            adata = self.merge_cluster_data(adata, cluster_data)

        # Get configuration from cellbrowser section
        color_field = self.cellbrowser_config.get(
            "color_field_cell2location", "region_cluster"
        )
        label_field = self.cellbrowser_config.get(
            "label_field_cell2location", "region_cluster"
        )
        enum_fields = self.cellbrowser_config.get(
            "enum_fields_cell2location", ["region_cluster", "batch"]
        )

        # Add cluster fields to enum_fields
        cluster_fields = list(cluster_data.keys()) if cluster_data else []
        enum_fields.extend(cluster_fields)
        enum_fields = list(set(enum_fields))  # Remove duplicates

        # Process data based on batch presence
        if "batch" in adata.obs.columns:
            print("Found batch column - combining all batches with shifted coordinates")
            # Shift spatial coordinates for tiled layout
            adata = self._shift_spatial_coordinates(adata)
            # Export all data as single combined batch
            self._export_batch(adata, color_field, label_field, enum_fields)
        else:
            # Export all data as single batch
            print("No batch column found - exporting all data")
            self._export_batch(adata, color_field, label_field, enum_fields)

    def _create_desc_config(self, output_dir: Path, data_type: str = "gene expression"):
        super()._create_desc_config(
            output_dir, data_type="cell type population density"
        )

    def _create_cellbrowser_config(
        self,
        output_dir: Path,
        color_field: str,
        label_field: str,
        enum_fields: List[str],
        coord_bounds: Dict,
        marker_files: List[Dict[str, str]],
        project_suffix: str = "cellTypes",
    ):
        super()._create_cellbrowser_config(
            output_dir,
            color_field,
            label_field,
            enum_fields,
            coord_bounds,
            marker_files,
            project_suffix=project_suffix,
        )

    def _export_quick_genes(self, adata: ad.AnnData, output_dir: Path):
        """Export quick cell types list (override for cell types instead of genes)"""
        quick_genes_file = output_dir / "quickGenes.tsv"

        # Get all available cell type names
        all_cell_types = adata.var_names.tolist()

        with open(quick_genes_file, "w") as f:
            for i, cell_type in enumerate(all_cell_types):
                # Write as cellType|cellType to suppress Cell Browser warning
                f.write(f"{cell_type}|{cell_type}\t{i}\n")

        print(f"Exported {len(all_cell_types)} quick cell types")


def main():
    """Main function to run the Cell Browser export"""
    # Load configuration
    config_file = Path("config/batch_config.yaml")
    if not config_file.exists():
        print(f"Configuration file not found: {config_file}")
        return

    with open(config_file, "r") as f:
        config = yaml.safe_load(f)

    # Get cellbrowser_html_output_dir from config, fallback to default if missing
    cellbrowser_html_output_dir = Path(config.get("cellbrowser_html_output_dir", "output/html_reports/cellbrowser")).resolve()

    # ===========================================================
    # ===== Create an exporter for gene expression data =====
    print("Creating Cell Browser exporter for gene expression levels...")
    gene_exporter = CellBrowserExporter(config)

    try:
        gene_exporter.export_for_cellbrowser()
        print("Cell Browser export completed successfully!")
        subprocess.run(
            [
                "cbBuild",
                "-o", cellbrowser_html_output_dir
            ],
            check=True,
            cwd=str(gene_exporter.cellbrowser_dir)
        )
    except Exception as e:
        print(f"Export failed: {e}")
        raise

    # ===========================================================
    # ==== create another exporter for cell abundances ====
    print("Creating Cell Browser exporter for cell type abundances...")
    cell_type_exporter = Cell2LocationCellBrowserExporter(config)

    try:
        cell_type_exporter.export_for_cellbrowser()
        print("Cell2Location Cell Browser export completed successfully!")
        subprocess.run(
            [
                "cbBuild",
                "-o", cellbrowser_html_output_dir
            ],
            check=True,
            cwd=str(cell_type_exporter.cellbrowser_dir)
        )
    except Exception as e:
        print(f"Export failed: {e}")
        raise

    # ===========================================================
    # === Check for combine_cell_types and export combined cell types ===
    # Load config/cellType_config.yaml for combine_cell_types and combinations

    celltype_config_path = Path("config/cellType_config.yaml")
    if not celltype_config_path.exists():
        print("Could not find config/cellType_config.yaml, skipping combined cell type export.")
        return

    with open(celltype_config_path, "r") as f:
        celltype_config = yaml.safe_load(f)

    combine_cell_types = celltype_config.get("shared", {}).get("combine_cell_types", False)
    cell_type_combinations = celltype_config.get("shared", {}).get("cell_type_combinations", {})

    if combine_cell_types and cell_type_combinations:
        print("combine_cell_types is enabled in config/cellType_config.yaml, exporting combined cell type abundances for Cell Browser...")

        # Load the original cell abundance CSV
        abundances_csv_path = config.get("cellbrowser", {}).get(
            "cell_abundances_csv",
            "output/cell_type_prediction/cell2location_results/cell_abundances_and_clusters.csv",
        )
        abundances_csv_path = Path(abundances_csv_path)
        if not abundances_csv_path.exists():
            print(f"Cell abundance CSV not found: {abundances_csv_path}")
            return

        df = pd.read_csv(abundances_csv_path, header=0, index_col=0)
        # Identify metadata columns
        metadata_cols = [ meta_col for meta_col in ["region_cluster", "batch"] if meta_col in df.columns ]
        # Combine cell types as in pred_cell2location_2_map_spatial.py
        combined_cell_types = []
        combined_subtypes = set()
        for combined, subtypes in cell_type_combinations.items():
            valid_subtypes = [s for s in subtypes if s in df.columns]
            if not valid_subtypes:
                continue
            df[combined] = df[valid_subtypes].sum(axis=1)
            combined_cell_types.append(combined)
            combined_subtypes.update(valid_subtypes)

        # Also keep other cell subtypes (not those being combined)
        all_cell_subtypes = [col for col in df.columns if col not in metadata_cols]
        other_cell_subtypes = [col for col in all_cell_subtypes if col not in combined_subtypes and col not in combined_cell_types]

        # Only keep metadata columns + combined cell types + other cell subtypes
        keep_cols = metadata_cols + combined_cell_types + other_cell_subtypes
        df_combined = df[keep_cols].copy()

        # Write to a new CSV for combined cell types
        combined_csv_path = abundances_csv_path.parent / "cell_abundances_and_clusters_combined.csv"
        df_combined.to_csv(combined_csv_path)
        print(f"Combined cell type abundance CSV written to: {combined_csv_path}")

        # Patch config for the combined run
        config_combined = copy.deepcopy(config)
        # Set new CSV path and output dir
        config_combined.setdefault("cellbrowser", {})
        config_combined["cellbrowser"]["cell_abundances_csv"] = str(combined_csv_path)
        # Use a new output dir for combined cell types
        orig_dir = config_combined["cellbrowser"].get("cell_abundances_dir", "output/cellbrowser_cellTypes")
        combined_dir = orig_dir + "_combined"
        config_combined["cellbrowser"]["cell_abundances_dir"] = combined_dir

        # Only use combined cell types + other cell subtypes as quickGenes
        config_combined["cellbrowser"]["quick_genes"] = combined_cell_types + other_cell_subtypes

        # Run Cell2LocationCellBrowserExporter for combined cell types
        print("Creating Cell Browser exporter for combined cell type abundances...")

        # Patch Cell2LocationCellBrowserExporter to use project_suffix="cellTypes_combine"
        class CombinedCell2LocationCellBrowserExporter(Cell2LocationCellBrowserExporter):
            def _create_cellbrowser_config(
                self,
                output_dir: Path,
                color_field: str,
                label_field: str,
                enum_fields: List[str],
                coord_bounds: Dict,
                marker_files: List[Dict[str, str]],
            ):
                super()._create_cellbrowser_config(
                    output_dir,
                    color_field,
                    label_field,
                    enum_fields,
                    coord_bounds,
                    marker_files,
                    project_suffix="cellTypes_combine",
                )

        cell_type_exporter_combined = CombinedCell2LocationCellBrowserExporter(config_combined)
        try:
            cell_type_exporter_combined.export_for_cellbrowser()
            print("Cell2Location Cell Browser export for combined cell types completed successfully!")
            subprocess.run(
                [
                    "cbBuild",
                    "-o", cellbrowser_html_output_dir
                ],
                check=True,
                cwd=str(cell_type_exporter_combined.cellbrowser_dir)
            )
        except Exception as e:
            print(f"Export failed for combined cell types: {e}")
            raise

    # ========== MOVE ALL DELETE PROMPTS TO HERE ==========
    # Prompt to ask user whether to delete the .cellbrowser_dir for gene_exporter
    resp = input(f"Do you want to delete the directory '{gene_exporter.cellbrowser_dir}'? [y/N]: ")
    if resp.lower() == "y":
        shutil.rmtree(gene_exporter.cellbrowser_dir)
        print(f"Deleted directory: {gene_exporter.cellbrowser_dir}")

    # Prompt to ask user whether to delete the .cellbrowser_dir for cell_type_exporter
    resp = input(f"Do you want to delete the directory '{cell_type_exporter.cellbrowser_dir}'? [y/N]: ")
    if resp.lower() == "y":
        shutil.rmtree(cell_type_exporter.cellbrowser_dir)
        print(f"Deleted directory: {cell_type_exporter.cellbrowser_dir}")

    # If combined exporter was created, prompt to delete its directory
    if combine_cell_types and cell_type_combinations:
        resp = input(f"Do you want to delete the directory '{cell_type_exporter_combined.cellbrowser_dir}'? [y/N]: ")
        if resp.lower() == "y":
            shutil.rmtree(cell_type_exporter_combined.cellbrowser_dir)
            print(f"Deleted directory: {cell_type_exporter_combined.cellbrowser_dir}")


if __name__ == "__main__":
    main()
