# -*- coding: utf-8 -*-

import seaborn as sns
import os
from scipy import sparse
import tempfile
import zarr
import dask.array as da
from typing import Optional, List, Dict, Tuple
import concurrent.futures
from anndata.io import read_text

import logging
import gc
import psutil
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import plotly.express as px
from plotly.offline import plot
import spatialdata_io as sdio
from pathlib import Path
import yaml
from scipy.sparse import csr_matrix, issparse
import warnings
import scanpy.external as sce 
from contextlib import contextmanager

# Configure warnings
warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", category=FutureWarning, module="dask.dataframe")
warnings.filterwarnings("ignore", category=FutureWarning, module="anndata")

class MemoryManager:
    """Memory monitoring and management utilities."""
    
    @staticmethod
    def get_memory_usage() -> float:
        """Get current memory usage in GB."""
        return psutil.virtual_memory().used / (1024**3)
    
    @staticmethod
    def get_available_memory() -> float:
        """Get available memory in GB."""
        return psutil.virtual_memory().available / (1024**3)
    
    @contextmanager
    def memory_monitor(self, operation_name: str):
        """Context manager to monitor memory usage during operations."""
        start_memory = self.get_memory_usage()
        logging.info(f"Starting {operation_name} - Memory usage: {start_memory:.2f} GB")
        try:
            yield
        finally:
            end_memory = self.get_memory_usage()
            logging.info(f"Completed {operation_name} - Memory usage: {end_memory:.2f} GB (Δ: {end_memory-start_memory:+.2f} GB)")
            gc.collect()

class VisiumProcessor:
    """
    Optimized processor for large Visium HD datasets with on-disk operations.
    """
    
    def _load_config(self, config_path: str) -> dict:
        """Load configuration with validation."""
        try:
            with open(config_path, 'r') as f:
                config = yaml.load(f, Loader=yaml.FullLoader)
            
            # Validate required keys
            required_keys = ['project_dir', 'raw_data_dir', 'batch_file_names', 'batch_names']
            for key in required_keys:
                if key not in config:
                    raise ValueError(f"Missing required config key: {key}")
            
            return config
        except Exception as e:
            logging.error(f"Error loading config: {e}")
            raise

    def __init__(self, config_path: str):
        self.config = self._load_config(config_path)
        self.project_dir = Path(self.config['project_dir'])
        self.output_dirs = {k: Path(v) for k, v in self.config['output_dirs'].items()}
        self._setup_directories()
        self.memory_manager = MemoryManager()
        self.sdata_objects = []  # Store SpatialData objects for on-disk access
        self.batch_metadata = []
        self.merged_adata_path = None
        
        # Configure scanpy for on-disk operations
        # Use values from config['on_disk'] if present, otherwise fallback
        on_disk = self.config.get('on_disk', {})
        sc.settings.cache = True
        sc.settings.max_memory = on_disk.get('max_memory_gb', 45)
        
        self._setup_logging()

    def _setup_logging(self):
        """Setup enhanced logging."""
        log_file = self.project_dir / "visium_optimized_processing.log"
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger(__name__)

    def _setup_directories(self):
        """Create necessary directories."""
        self.zarr_cache_dir = self.project_dir / "zarr_cache"
        self.zarr_cache_dir.mkdir(parents=True, exist_ok=True)
        
        for dir_path in self.output_dirs.values():
            dir_path.mkdir(parents=True, exist_ok=True)

    def load_data_optimized(self):
        """Load data with enhanced memory management and on-disk operations."""
        self.logger.info("=== Loading Visium Data with Optimization ===")
        
        for i, (batch_file_name, batch_name) in enumerate(zip(
            self.config['batch_file_names'], 
            self.config['batch_names']
        )):
            with self.memory_manager.memory_monitor(f"Loading batch {batch_name}"):
                try:
                    self._load_single_batch(batch_file_name, batch_name, i)
                except Exception as e:
                    self.logger.error(f"Failed to load batch {batch_name}: {e}")
                    continue

        self._create_batch_summary()

    def _load_single_batch(self, batch_file_name: str, batch_name: str, batch_idx: int):
        """Load a single batch with optimization."""
        self.logger.info(f"Loading batch {batch_idx+1}: {batch_name}")
        
        full_path = Path(self.config['raw_data_dir']) / batch_file_name
        zarr_cache_path = self.zarr_cache_dir / f"{batch_name}_cache.zarr"
        
        # Check if cached version exists
        if zarr_cache_path.exists() and self.config['on_disk'].get('use_zarr_cache', True):
            self.logger.info(f"Loading cached data for {batch_name}")
            adata = ad.read_zarr(zarr_cache_path)
        else:
            # Load fresh data
            if self.config.get('VisiumHD', False):
                adata = self._load_visium_hd(full_path, batch_name)
            else:
                adata = self._load_standard_visium(full_path, batch_name, batch_idx)
            
            # Cache the processed data
            if self.config['on_disk'].get('use_zarr_cache', True):
                adata.write_zarr(zarr_cache_path)
                self.logger.info(f"Cached data for {batch_name}")

        # Basic preprocessing with memory management
        self._preprocess_batch(adata, batch_name)
        
        # Store metadata
        self.batch_metadata.append({
            'batch': batch_name,
            'n_cells': adata.n_obs,
            'n_genes': adata.n_vars,
            'total_counts': float(adata.X.sum()) if issparse(adata.X) else float(np.sum(adata.X)),
            'zarr_path': zarr_cache_path
        })

    def _load_visium_hd(self, full_path: Path, batch_name: str) -> ad.AnnData:
        """Load Visium HD data using spatialdata-io with on-disk capabilities."""
        full_path_outs = full_path / "outs"
        
        # Use spatialdata-io with on-disk mode
        sdata = sdio.visium_hd(
            path=full_path_outs,
            dataset_id=batch_name,
            load_all_images=False,  # Don't load images to save memory
        )
        
        # Extract AnnData object
        adata_key_candidates = [k for k in sdata.tables.keys() 
                               if any(term in k.lower() for term in ['counts', 'visiumhd', 'table'])]
        adata_key = adata_key_candidates[0] if adata_key_candidates else list(sdata.tables.keys())[0]
        
        adata = sdata.tables[adata_key]
        adata.obs['batch'] = batch_name
        adata.obs['orig.ident'] = batch_name
        
        # Store spatial data reference
        self.sdata_objects.append((batch_name, sdata))
        
        return adata

    def _load_standard_visium(self, full_path: Path, batch_name: str, batch_idx: int) -> ad.AnnData:
        """Load standard Visium data."""
        h5_file_names = self.config.get('batch_file_h5_names', [None] * len(self.config['batch_file_names']))
        h5_file_name = h5_file_names[batch_idx]
        
        if h5_file_name:
            adata = sq.read.visium(str(full_path.parent), count_file=h5_file_name)
        else:
            adata = sq.read.visium(str(full_path))
        
        adata.obs['batch'] = batch_name
        adata.obs['orig.ident'] = batch_name
        
        return adata

    def _preprocess_batch(self, adata: ad.AnnData, batch_name: str):
        """Preprocess individual batch with memory optimization."""
        # Store raw counts
        if 'counts' not in adata.layers:
            adata.layers['counts'] = adata.X.copy()
        
        # Normalize and log-transform
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        
        # Store normalized data
        adata.raw = adata
        
        # Find highly variable genes with chunking for large datasets
        if adata.n_obs > 50000:
            # Use chunked HVG calculation for very large datasets
            self._calculate_hvg_chunked(adata)
        else:
            sc.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=2000)

    def _calculate_hvg_chunked(self, adata: ad.AnnData):
        """Calculate highly variable genes using chunked operations."""
        chunk_size = self.config['on_disk']['chunk_size']
        n_chunks = (adata.n_obs - 1) // chunk_size + 1
        
        # Calculate mean and variance per gene across chunks
        means = np.zeros(adata.n_vars)
        variances = np.zeros(adata.n_vars)
        
        for i in range(n_chunks):
            start_idx = i * chunk_size
            end_idx = min((i + 1) * chunk_size, adata.n_obs)
            
            chunk_data = adata.X[start_idx:end_idx]
            if issparse(chunk_data):
                chunk_means = np.array(chunk_data.mean(axis=0)).flatten()
                chunk_vars = np.array(chunk_data.power(2).mean(axis=0)).flatten() - chunk_means**2
            else:
                chunk_means = np.mean(chunk_data, axis=0)
                chunk_vars = np.var(chunk_data, axis=0)
            
            weight = (end_idx - start_idx) / adata.n_obs
            means += weight * chunk_means
            variances += weight * chunk_vars
        
        # Calculate highly variable genes based on mean-variance relationship
        # Simplified version - in practice, use scanpy's method on subsampled data
        mean_var_ratio = variances / (means + 1e-8)
        hvg_threshold = np.percentile(mean_var_ratio, 90)
        adata.var['highly_variable'] = mean_var_ratio > hvg_threshold

    def merge_data_optimized(self):
        """Merge data using on-disk operations."""
        self.logger.info("=== Merging Data with Optimization ===")
        
        with self.memory_manager.memory_monitor("Data merging"):
            # Create merged zarr store
            merged_zarr_path = self.zarr_cache_dir / "merged_data.zarr"
            
            # Load all batches and concatenate
            adata_list = []
            for metadata in self.batch_metadata:
                adata = ad.read_zarr(metadata['zarr_path'])
                adata_list.append(adata)
            
            # Concatenate with memory management
            merged_adata = ad.concat(
                adata_list, 
                label='batch', 
                keys=[m['batch'] for m in self.batch_metadata],
                index_unique='_',
                merge='same'
            )
            
            # Save merged data
            merged_adata.write_zarr(merged_zarr_path)
            self.merged_adata_path = merged_zarr_path
            
            # Clean up
            del adata_list, merged_adata
            gc.collect()

    def integrate_and_cluster_optimized(self):
        """Optimized integration and clustering with on-disk operations."""
        if not self.merged_adata_path:
            raise ValueError("Merged data not available")
        
        self.logger.info("=== Integration and Clustering with Optimization ===")
        
        with self.memory_manager.memory_monitor("Integration and clustering"):
            # Load merged data
            adata = ad.read_zarr(self.merged_adata_path)
            
            # Find HVGs on merged data
            sc.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=2000)
            
            # Optimized scaling
            self._scale_data_optimized(adata)
            
            # PCA with memory management
            sc.tl.pca(adata, n_comps=self.config['on_disk']['integration_pca_dims'], 
                     use_highly_variable=True, svd_solver='arpack')
            
            # Integration
            self._run_integration(adata)
            
            # Clustering
            self._run_clustering(adata)
            
            # UMAP
            sc.tl.umap(adata, min_dist=0.3, n_components=3, spread=1.0)
            
            # Save final result
            final_path = self.zarr_cache_dir / f"{self.config['output_file_prefix']}_final.zarr"
            adata.write_zarr(final_path)
            self.final_adata_path = final_path
            
            self.logger.info(f"Saved final results to {final_path}")

    def _scale_data_optimized(self, adata: ad.AnnData):
        """Optimized scaling using chunked operations with dynamic chunk size."""
        self.logger.info("Scaling data with dynamic chunked operations (low-memory mode)")
        gc.collect()

        hvg_mask = adata.var['highly_variable']
        hvg_indices = np.where(hvg_mask)[0]

        # Calculate scaling parameters using HVGs only
        if issparse(adata.X):
            from sklearn.utils.sparsefuncs import mean_variance_axis
            means, variances = mean_variance_axis(adata.X[:, hvg_indices], axis=0)
            stds = np.sqrt(variances)
        else:
            means = np.mean(adata.X[:, hvg_indices], axis=0)
            stds = np.std(adata.X[:, hvg_indices], axis=0)

        stds = np.clip(stds, 1e-8, None)

        # Use AnnData backed mode if possible
        if hasattr(adata, 'filename') and adata.filename is not None:
            adata.file.close()
            adata = ad.read_zarr(self.merged_adata_path, backed='r+')
            self.logger.info("AnnData loaded in backed mode for scaling.")

        # Dynamically determine chunk_size based on available memory
        mem = psutil.virtual_memory()
        available_gb = mem.available / (1024 ** 3)
        # Use up to 30% of available memory for a chunk
        target_gb = max(available_gb * 0.3, 1.0)  # at least 1GB

        # Estimate bytes per cell for HVG columns
        n_hvg = len(hvg_indices)
        dtype_size = np.dtype(adata.X.dtype).itemsize if hasattr(adata.X, 'dtype') else 8
        bytes_per_cell = n_hvg * dtype_size
        chunk_size_dynamic = int((target_gb * (1024 ** 3)) // bytes_per_cell)
        chunk_size_config = self.config['on_disk'].get('scaling_chunk_size', 5000)
        chunk_size = max(chunk_size_config, min(chunk_size_dynamic, adata.n_obs))

        self.logger.info(f"Dynamic chunk size for scaling: {chunk_size} (config: {chunk_size_config}, available_gb: {available_gb:.2f})")

        n_chunks = (adata.n_obs - 1) // chunk_size + 1

        for i in range(n_chunks):
            start_idx = i * chunk_size
            end_idx = min((i + 1) * chunk_size, adata.n_obs)
            # Extract chunk matrix
            if issparse(adata.X):
                chunk_hvg = adata.X[start_idx:end_idx, hvg_indices].toarray()
            else:
                chunk_hvg = np.array(adata.X[start_idx:end_idx, hvg_indices])
            chunk_hvg = (chunk_hvg - means) / stds
            chunk_hvg = np.clip(chunk_hvg, -10, 10)
            # Assign scaled chunk back
            if issparse(adata.X):
                for j, idx in enumerate(hvg_indices):
                    adata.X[start_idx:end_idx, idx] = chunk_hvg[:, j][:, np.newaxis]
            else:
                adata.X[start_idx:end_idx, hvg_indices] = chunk_hvg
            self.logger.info(f"Scaled chunk {i+1}/{n_chunks} ({end_idx-start_idx} cells)")
            gc.collect()

        gc.collect()
        self.logger.info("Scaling complete.")

    def _run_integration(self, adata: ad.AnnData):
        """Run Harmony integration with error handling."""
        try:
            sce.pp.harmony_integrate(adata, key='batch', basis='X_pca')
            self.integration_key = 'X_pca_harmony'
            self.logger.info("Harmony integration completed")
        except Exception as e:
            self.logger.warning(f"Harmony failed: {e}. Using PCA only.")
            self.integration_key = 'X_pca'

    def _run_clustering(self, adata: ad.AnnData):
        """Run clustering with optimized parameters."""
        n_neighbors = self.config.get('Cluster_n_neighbors', 15)
        resolution = self.config.get('Cluster_resolution', 1.0)
        
        sc.pp.neighbors(
            adata, 
            n_neighbors=n_neighbors, 
            n_pcs=self.config['on_disk']['integration_pca_dims'],
            use_rep=self.integration_key,
            method='umap'  # Use UMAP for neighbor calculation
        )
        
        sc.tl.leiden(adata, resolution=resolution, key_added='leiden_cluster')
        n_clusters = adata.obs['leiden_cluster'].nunique()
        self.logger.info(f"Found {n_clusters} clusters")

    def _create_batch_summary(self):
        """Create comprehensive batch summary."""
        if not self.batch_metadata:
            return
        
        df = pd.DataFrame(self.batch_metadata)
        
        # Add quality metrics
        low_count_threshold = self.config.get('low_count_threshold', 50000)
        df['low_counts'] = df['total_counts'] < low_count_threshold
        
        # Save summary
        summary_path = self.output_dirs['clustering'] / "batch_summary_optimized.csv"
        df.to_csv(summary_path, index=False)
        
        self.logger.info(f"Batch summary saved to {summary_path}")
        print("\nOptimized Batch Summary:")
        print(df)

class OptimizedVisualizer:
    """Optimized visualizer with memory-efficient plotting."""
    
    def __init__(self, processor: VisiumProcessor):
        self.processor = processor
        self.config = processor.config
        self.output_dirs = processor.output_dirs
        self.logger = processor.logger
        self.memory_manager = processor.memory_manager

    def create_visualizations(self):
        """Create all visualizations with memory management."""
        if not hasattr(self.processor, 'final_adata_path'):
            self.logger.error("No final data available for visualization")
            return
        
        self.logger.info("=== Creating Optimized Visualizations ===")
        
        with self.memory_manager.memory_monitor("Visualization"):
            # Load final data
            adata = ad.read_zarr(self.processor.final_adata_path)
            
            # Create UMAP plots
            self._create_umap_plots(adata)
            
            # Create batch-specific visualizations
            self._create_batch_visualizations(adata)

    def _create_umap_plots(self, adata: ad.AnnData):
        """Create optimized UMAP plots."""
        # 2D UMAP
        fig, ax = plt.subplots(figsize=(10, 8))
        sc.pl.umap(adata, color='leiden_cluster', ax=ax, show=False, 
                   frameon=False, legend_loc='on data')
        plt.savefig(self.output_dirs['clustering_umap_spatial'] / 
                   f"{self.config['output_file_prefix']}_umap_2d_optimized.png", 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        # 3D interactive plot
        if adata.obsm['X_umap'].shape[1] >= 3:
            self._create_3d_umap(adata)

    def _create_3d_umap(self, adata: ad.AnnData):
        """Create 3D UMAP plot with downsampling for performance."""
        # Downsample for interactive plotting
        max_cells = 10000
        if adata.n_obs > max_cells:
            sc.pp.subsample(adata, n_obs=max_cells, random_state=42)
        
        umap_coords = adata.obsm['X_umap'][:, :3]
        df_plot = pd.DataFrame(umap_coords, columns=['UMAP1', 'UMAP2', 'UMAP3'])
        df_plot['Cluster'] = adata.obs['leiden_cluster'].astype(str)
        df_plot['Batch'] = adata.obs['batch']
        
        fig = px.scatter_3d(df_plot, x='UMAP1', y='UMAP2', z='UMAP3', 
                           color='Cluster', hover_data=['Batch'],
                           title='3D UMAP - Optimized')
        
        output_path = self.output_dirs['clustering_umap_spatial'] / \
                     f"{self.config['output_file_prefix']}_umap_3d_optimized.html"
        plot(fig, filename=str(output_path), auto_open=False)
        
        self.logger.info(f"Saved 3D UMAP to {output_path}")

    def _create_batch_visualizations(self, adata: ad.AnnData):
        """Create batch-specific visualizations with memory management."""
        for batch_name in adata.obs['batch'].unique():
            with self.memory_manager.memory_monitor(f"Batch viz: {batch_name}"):
                self._process_single_batch_viz(adata, batch_name)

    def _process_single_batch_viz(self, adata: ad.AnnData, batch_name: str):
        """Process visualizations for a single batch."""
        # Subset data
        batch_mask = adata.obs['batch'] == batch_name
        adata_batch = adata[batch_mask].copy()
        
        if adata_batch.n_obs == 0:
            return
        
        # Create output directory
        batch_dir = self.output_dirs['clustering_umap_spatial'] / batch_name
        batch_dir.mkdir(exist_ok=True)
        
        # Find markers with memory management
        self._find_markers_optimized(adata_batch, batch_dir, batch_name)
        
        # Create spatial plots if coordinates available
        if 'spatial' in adata_batch.obsm:
            self._create_spatial_plots(adata_batch, batch_dir, batch_name)

    def _find_markers_optimized(self, adata: ad.AnnData, output_dir: Path, batch_name: str):
        """Find markers with memory optimization."""
        try:
            # Downsample if too large
            if adata.n_obs > 5000:
                sc.pp.subsample(adata, n_obs=5000, random_state=42)
            
            sc.tl.rank_genes_groups(adata, 'leiden_cluster', 
                                   method='wilcoxon', use_raw=True)
            
            # Save results
            markers = sc.get.rank_genes_groups_df(adata, group=None)
            markers.to_csv(output_dir / f"markers_{batch_name}_optimized.csv", index=False)
            
            self.logger.info(f"Saved markers for {batch_name}")
            
        except Exception as e:
            self.logger.error(f"Error finding markers for {batch_name}: {e}")

    def _create_spatial_plots(self, adata: ad.AnnData, output_dir: Path, batch_name: str):
        """Create spatial plots with optimization."""
        try:
            # Interactive spatial plot
            spatial_coords = adata.obsm['spatial']
            df_spatial = pd.DataFrame(spatial_coords, columns=['X', 'Y'])
            df_spatial['Cluster'] = adata.obs['leiden_cluster'].astype(str)
            
            fig = px.scatter(df_spatial, x='X', y='Y', color='Cluster',
                           title=f"Spatial Clusters - {batch_name} (Optimized)")
            fig.update_yaxes(autorange="reversed")
            
            spatial_path = output_dir / f"spatial_{batch_name}_optimized.html"
            plot(fig, filename=str(spatial_path), auto_open=False)
            
            self.logger.info(f"Saved spatial plot for {batch_name}")
            
        except Exception as e:
            self.logger.error(f"Error creating spatial plot for {batch_name}: {e}")


def main():
    """Main execution with enhanced error handling."""
    config_path = "config/batch_config.yaml"
    
    try:
        # Initialize processor
        processor = VisiumProcessor(config_path)

        merged_zarr_path = processor.zarr_cache_dir / "merged_data.zarr"
        load_from_merge = processor.config.get('on_disk', {}).get('load_from_merge', False)

        if merged_zarr_path.exists() and load_from_merge:
            processor.logger.info("Merged data found and load_from_merge=True; skipping load_data_optimized.")
            processor.merged_adata_path = merged_zarr_path
        else:
            processor.load_data_optimized()
            processor.merge_data_optimized()

        processor.integrate_and_cluster_optimized()
        
        # Create visualizations
        visualizer = OptimizedVisualizer(processor)
        visualizer.create_visualizations()
        
        processor.logger.info("=== Optimized Processing Complete ===")
        print("\n" + "="*50)
        print("OPTIMIZED PROCESSING COMPLETE")
        print("="*50)
        
    except Exception as e:
        logging.error(f"Critical error in main execution: {e}")
        raise

if __name__ == "__main__":
    main()
