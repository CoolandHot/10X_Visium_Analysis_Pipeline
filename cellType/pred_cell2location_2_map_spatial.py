'''
Refactored Cell2Location Spatial Mapping Script.

This script provides a class-based, modular structure for running the
Cell2Location pipeline. It maintains all original functionality, including
dataset splitting for distributed training with sample batch preservation and
cell overlap.

=== OVERVIEW ===

Cell2Location is a Bayesian model for comprehensive mapping of tissue architecture 
via integrative analysis of single-cell and spatial transcriptomics data. This script
handles the spatial mapping phase, which uses pre-trained reference cell type 
signatures to deconvolve spatial transcriptomics data.

Key Features:
- Full dataset processing or distributed batch processing
- Dataset splitting with stratified sampling and overlap handling
- Automatic model training/loading with incremental saves
- Comprehensive visualization generation
- Memory-efficient batch merging
- Robust error handling and data validation

=== MODES OF OPERATION ===

1. FULL DATASET MODE (Default)
   - Processes the entire spatial dataset in one go
   - Suitable for smaller datasets or high-memory systems
   - Command: python pred_cell2location_2_map_spatial.py

2. DATASET SPLITTING MODE
   - Splits large datasets into smaller batches for distributed processing
   - Preserves sample batch proportions and adds cell overlap between splits
   - Command: python pred_cell2location_2_map_spatial.py --split-only

3. BATCH PROCESSING MODE
   - Processes individual batches on separate machines/processes
   - Each batch includes overlapping cells for continuity
   - Command: python pred_cell2location_2_map_spatial.py --batch N

4. BATCH MERGING MODE
   - Combines results from all processed batches into unified dataset
   - Handles overlapping cells and restores spatial metadata
   - Command: python pred_cell2location_2_map_spatial.py --merge

5. VISUALIZATION MODE
   - Generates plots from existing merged results without reprocessing
   - Useful for re-running visualizations after fixing errors
   - Command: python pred_cell2location_2_map_spatial.py --visualize

=== TYPICAL WORKFLOW ===

For Large Datasets (Distributed Processing):
1. python pred_cell2location_2_map_spatial.py --split-only
2. python pred_cell2location_2_map_spatial.py --batch 0  # Machine 1
   python pred_cell2location_2_map_spatial.py --batch 1  # Machine 2
   python pred_cell2location_2_map_spatial.py --batch N  # Machine N+1
3. python pred_cell2location_2_map_spatial.py --merge
4. python pred_cell2location_2_map_spatial.py --visualize  # If needed

For Small/Medium Datasets (Single Machine):
1. python pred_cell2location_2_map_spatial.py

=== CONFIGURATION ===

The script reads configuration from config/cellType_config.yaml, which should contain:

paths:
  spatial_data: "path/to/spatial_data.h5ad"
  reference_output: "path/to/reference_signatures"
  spatial_output: "path/to/output"
  cell_abundance_heatmap_path: "path/to/plots"

dataset:
  sample_batch_key: "orig.ident"  # Column name for sample batches
  available_slices: ["slice1", "slice2", ...]
  splitting:
    enabled: true
    batch_size: 10000
    overlap_cells: 1000
    split_data_dir: "path/to/split_data"

cell2location:
  N_cells_per_location: 30
  detection_alpha: 20
  max_epochs: 30000
  # ... other parameters

=== INPUT REQUIREMENTS ===

1. Spatial transcriptomics data (AnnData format)
   - Raw count data (not log-transformed)
   - Spatial coordinates in .obsm['spatial']
   - Sample batch information in .obs[sample_batch_key]

2. Pre-trained reference signatures
   - From pred_cell2location_1_reference.py
   - CSV file with cell type signatures

=== OUTPUT STRUCTURE ===

For full dataset processing:
output_dir/
├── sp.h5ad                           # Base processed data
├── sp_model_increments.pkl           # Model posterior additions
├── sp_final_increments.pkl           # Final processing additions
├── model.pt                          # Trained Cell2Location model
├── spatial_training_history.png      # Training convergence plot
├── cell_abundances_and_clusters.csv  # Results summary
└── cell_abundance_heatmap/           # Visualization outputs
    ├── slice1_cell_abundance_all.pdf
    ├── slice1_selected_cell_types.pdf
    ├── by_cell_type/
    │   ├── CellType1_all_batches.pdf
    │   └── ...
    └── umap_clustering.pdf

For batch processing:
output_dir/
├── batch_0/                          # First batch results
├── batch_1/                          # Second batch results
├── ...
└── merged_results/                   # Combined results
    ├── merged_adata.h5ad
    ├── merged_cell_abundances_and_clusters.csv
    └── [visualization files]

=== ERROR HANDLING ===

The script includes robust error handling for:
- Missing input files
- Data format inconsistencies
- Memory limitations during processing
- NaN values in plotting data
- Missing spatial coordinates
- Incomplete batch results

Common issues and solutions:
- "cannot convert float NaN to integer": Fixed with NaN handling in plotting
- "Could not find 'umap'": UMAP coordinates missing, will be skipped
- Memory errors: Use batch processing mode
- Missing reference signatures: Check reference model training

=== TECHNICAL DETAILS ===

Data Processing:
- Automatic detection and handling of log-transformed data
- Mitochondrial gene filtering
- Gene intersection between spatial and reference data
- Data validation with quality metrics

Model Training:
- Incremental saving to handle interruptions
- Automatic model loading if already trained
- GPU/CPU acceleration support
- Configurable training parameters

Visualization:
- Multi-page PDF outputs for large cell type sets
- Adaptive subplot layouts based on sample count
- Consistent color scaling across samples
- Cell2Location's advanced spatial plotting integration

Memory Management:
- Lean data structures for batch merging
- Incremental processing to minimize memory usage
- Garbage collection at strategic points
- Spatial metadata preservation during merging

Original reference documentation:
https://github.com/BayraktarLab/cell2location/blob/master/docs/notebooks/cell2location_short_demo_colab.ipynb
https://github.com/BayraktarLab/cell2location/blob/master/docs/notebooks/cell2location_short_demo_downstream.ipynb
https://github.com/BayraktarLab/cell2location/blob/master/docs/notebooks/cell2location_tutorial.ipynb

Author: Refactored for modular, distributed processing
Version: 2.0
Last Updated: 2024
'''
import scipy.sparse
import scanpy as sc
import numpy as np
import torch
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
import os
import yaml
import cell2location
from cell2location.plt import plot_spatial
from pred_cell2location_utils import load_configuration, extract_adata_additions, restore_adata_additions, load_processed_spatial_data, load_incremental_adata
import anndata as ad
import argparse
import json
from copy import deepcopy

# Set Matplotlib defaults
rcParams = plt.rcParams
rcParams['pdf.fonttype'] = 42

class ConfigManager:
    """Handles loading and accessing the configuration from config/cellType_config.yaml."""
    def __init__(self):
        """Loads configuration and sets up paths."""
        print("Loading configuration...")
        self.config = load_configuration()
        batch_config = load_configuration(config_path='config/batch_config.yaml')
        self.paths = self.config['paths']
        self.config['dataset']['available_slices'] = batch_config.get('batch_names', [])
        self.dataset_config = self.config['dataset']
        self.splitting_config = self.dataset_config['splitting']

    def get_path(self, key):
        """Returns a path from the 'paths' section."""
        return self.paths[key]

    def get_config(self, section):
        """Returns a whole configuration section."""
        return self.config[section]


class DataProcessor:
    """Handles all data loading, validation, cleaning, and preparation."""
    def __init__(self, config_manager):
        self.config = config_manager.config
        self.dataset_config = self.config['dataset']
        self.paths = self.config['paths']
        self.ref_run_name_path = self.paths['reference_output']
    
    def load_initial_spatial_data(self):
        """Loads the initial GBM spatial data."""
        print("Loading initial spatial data...")
        adata_vis_initial = sc.read_h5ad(self.paths['spatial_data'])
        print(f"Initial spatial data shape: {adata_vis_initial.shape}")
        return adata_vis_initial

    def _validate_data_quality(self, adata, name="data"):
        """Validates data quality by checking for inf/nan and zero-sum features."""
        print(f"\nValidating {name} quality...")
        inf_count = np.isinf(adata.X.data).sum() if hasattr(adata.X, 'data') else np.isinf(adata.X).sum()
        nan_count = np.isnan(adata.X.data).sum() if hasattr(adata.X, 'data') else np.isnan(adata.X).sum()
        min_val = adata.X.min()
        max_val = adata.X.max()
        zero_genes = (np.array(adata.X.sum(axis=0)).flatten() == 0).sum()
        zero_cells = (np.array(adata.X.sum(axis=1)).flatten() == 0).sum()
        
        print(f"  - Infinite values in X: {inf_count}")
        print(f"  - NaN values in X: {nan_count}")
        print(f"  - Data range: [{min_val:.3f}, {max_val:.3f}]")
        print(f"  - Zero-sum genes: {zero_genes}")
        print(f"  - Zero-sum cells: {zero_cells}")

    def _ensure_count_data(self, adata):
        """Ensures that the AnnData object contains raw count data."""
        print("\nEnsuring data is in raw count format for Cell2Location...")
        if adata.raw is not None:
            print("Using raw counts from adata.raw")
            adata = adata.raw.to_adata()
        elif np.issubdtype(adata.X.dtype, np.floating) and adata.X.max() < 20:
             print("WARNING: Data appears to be log-transformed. Reverting with np.expm1.")
             adata.X = np.expm1(adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X)

        if scipy.sparse.issparse(adata.X):
            # Only clip negative values in-place for sparse matrix
            adata.X.data = np.clip(adata.X.data, 0, None)
        else:
            adata.X = np.clip(adata.X, 0, None)

        if not np.issubdtype(adata.X.dtype, np.integer):
            # For sparse, convert to dense before rounding and casting
            if scipy.sparse.issparse(adata.X):
                adata.X = np.round(adata.X.toarray()).astype(np.int32)
            else:
                adata.X = np.round(adata.X).astype(np.int32)
        print(f"Final data type: {adata.X.dtype}, Range: [{adata.X.min()}, {adata.X.max()}]")
        return adata

    def load_reference_and_prepare_spatial(self, adata_vis_input):
        """Loads reference signatures and prepares spatial data by finding intersecting genes."""
        print("\nLoading pre-trained reference model signatures...")
        inf_aver_ref = pd.read_csv(f"{self.ref_run_name_path}/reference_signatures.csv", index_col=0)
        
        # Clean reference signatures
        inf_aver_ref = inf_aver_ref.loc[~inf_aver_ref.index.duplicated(keep='first')]
        inf_aver_ref = inf_aver_ref.fillna(0).replace([np.inf, -np.inf], 0)

        # Prepare spatial data by removing mitochondrial genes
        adata_vis_input.var['MT_gene'] = [gene.startswith('MT-') or gene.startswith('mt-') for gene in adata_vis_input.var_names]
        adata_vis_input = adata_vis_input[:, ~adata_vis_input.var['MT_gene'].values].copy()
        
        adata_vis_input = self._ensure_count_data(adata_vis_input)
        self._validate_data_quality(adata_vis_input, "spatial input before filtering")
        
        sc.pp.filter_genes(adata_vis_input, min_cells=5)
        sc.pp.filter_cells(adata_vis_input, min_genes=10)
        self._validate_data_quality(adata_vis_input, "spatial cleaned")
        
        # Find and subset to intersecting genes
        intersect_genes = sorted(list(set(adata_vis_input.var_names) & set(inf_aver_ref.index)))
        if len(intersect_genes) < 100:
            raise ValueError(f"Found only {len(intersect_genes)} shared genes. Check gene name consistency.")
        print(f"Found {len(intersect_genes)} shared genes between spatial data and reference signatures.")

        adata_vis_processed = adata_vis_input[:, intersect_genes].copy()
        inf_aver_processed = inf_aver_ref.loc[intersect_genes, :].copy()

        # Final check to ensure gene order is identical
        if not np.all(adata_vis_processed.var_names == inf_aver_processed.index):
            raise ValueError("Gene names do not match exactly after intersection. This should not happen.")

        print("Data preparation complete.")
        cell2location.models.Cell2location.setup_anndata(adata=adata_vis_processed, batch_key=self.dataset_config['sample_batch_key'])
        return adata_vis_processed, inf_aver_processed


class ModelManager:
    """Handles model training, loading, saving, and post-processing."""
    def __init__(self, config_manager, visualizer, run_path):
        self.config = config_manager.config
        self.run_path = run_path
        self.visualizer = visualizer
        self.model_config = self.config['cell2location']
        self.clustering_config = self.config['clustering']
        self.dataset_config = self.config['dataset']
        
    def _create_and_train_model(self, adata_vis, inf_aver):
        """Creates and trains a new Cell2Location model."""
        print("Creating and training a new cell2location model...")
        mod = cell2location.models.Cell2location(
            adata_vis, 
            cell_state_df=inf_aver,
            N_cells_per_location=self.model_config['N_cells_per_location'], 
            detection_alpha=self.model_config['detection_alpha'] 
        )
        mod.view_anndata_setup()

        # Train model
        torch.set_float32_matmul_precision('medium')
        mod.train(
            max_epochs=self.model_config['max_epochs'], 
            batch_size=self.model_config['batch_size'], 
            train_size=self.model_config['train_size'], 
            accelerator=self.model_config['accelerator'], 
            device=self.model_config['device'] 
        )
        
        return mod

    def train_or_load_model(self, adata_vis_processed, inf_aver_processed, model_path, base_adata_path, increments_path):
        """Trains a new model or loads an existing one."""
        if os.path.exists(model_path):
            print(f"Found existing trained model. Loading from {self.run_path}.")
            adata_vis_base = sc.read_h5ad(base_adata_path)
            adata_vis_loaded = restore_adata_additions(adata_vis_base, increments_path)
            mod = cell2location.models.Cell2location.load(self.run_path, adata_vis_loaded)
            cell2location.models.Cell2location.setup_anndata(adata=mod.adata, batch_key=self.dataset_config['sample_batch_key'])
            return mod, adata_vis_loaded
        
        # Save base adata before training
        adata_vis_processed.write(base_adata_path)

        mod = self._create_and_train_model(adata_vis_processed, inf_aver_processed)
        
        print("Exporting posterior...")
        adata_vis_with_posterior = mod.export_posterior(
            adata_vis_processed.copy(), 
            sample_kwargs={
                'num_samples': self.model_config['num_samples'],
                'batch_size': mod.adata.n_obs,
                'accelerator': self.model_config['accelerator'],
                'device': self.model_config['device']}
        )
        
        print("Saving model and increments...")
        mod.save(self.run_path, overwrite=True)
        extract_adata_additions(adata_vis_processed, adata_vis_with_posterior, increments_path)
        
        self.visualizer.plot_training_history(mod)
        
        return mod, adata_vis_with_posterior

    def process_and_save_results(self, adata_vis_mapped, final_increments_path):
        """Performs clustering, UMAP, and saves results."""
        print("Performing post-training processing (clustering, UMAP)...")
        adata_before_processing = adata_vis_mapped.copy()
        
        # Add cell abundance to .obs
        factor_names = adata_vis_mapped.uns['mod']['factor_names']
        adata_vis_mapped.obs[factor_names] = adata_vis_mapped.obsm['q05_cell_abundance_w_sf']

        # Clustering
        sc.pp.neighbors(adata_vis_mapped, use_rep='q05_cell_abundance_w_sf', n_neighbors=self.clustering_config['n_neighbors'])
        sc.tl.leiden(adata_vis_mapped, resolution=self.clustering_config['leiden_resolution'])
        adata_vis_mapped.obs["region_cluster"] = adata_vis_mapped.obs["leiden"].astype("category")

        # UMAP
        sc.tl.umap(adata_vis_mapped, min_dist=self.clustering_config['umap_min_dist'], spread=self.clustering_config['umap_spread'])

        # Save results to CSV
        columns_to_export = factor_names + ["region_cluster", self.dataset_config['sample_batch_key']]
        csv_path = f"{self.run_path}/cell_abundances_and_clusters.csv"
        adata_vis_mapped.obs[columns_to_export].to_csv(csv_path)
        print(f"Saved cell abundances and clusters to {csv_path}")

        # Save final processing increments
        extract_adata_additions(adata_before_processing, adata_vis_mapped, final_increments_path)
        print(f"Saved final processing increments to {final_increments_path}")
        return adata_vis_mapped

class Visualizer:
    """Handles all plotting and visualization tasks."""
    def __init__(self, config_manager, run_path):
        self.config = config_manager.config
        self.run_path = run_path
        self.vis_config = self.config['visualization']
        self.dataset_config = self.config['dataset']
        self.paths = self.config['paths']

    @staticmethod
    def generate_visualizations_for_adata(adata, config, output_path, combine_cell_types=False, cell_type_combinations=None):
        """
        Shared visualization logic for merged results.
        Optionally combines cell types as specified.
        """
        from types import SimpleNamespace
        temp_config = SimpleNamespace()
        temp_config_dict = deepcopy(config)
        combined_cell_types = []
        combined_heatmap_path = None

        if combine_cell_types and cell_type_combinations:
            print("Combining subtypes into new cell types as specified in config (for visualization)...")
            def combine_cell_type_abundances(adata, cell_type_combinations):
                combined_cell_types = []
                for combined, subtypes in cell_type_combinations.items():
                    valid_subtypes = [s for s in subtypes if s in adata.obs.columns]
                    if not valid_subtypes:
                        continue
                    adata.obs[combined] = adata.obs[valid_subtypes].sum(axis=1)
                    combined_cell_types.append(combined)
                return combined_cell_types
            combined_cell_types = combine_cell_type_abundances(adata, cell_type_combinations)
            combined_heatmap_path = config['paths']['cell_abundance_heatmap_path'] + '_combined'
            os.makedirs(combined_heatmap_path, exist_ok=True)
            temp_config_dict['paths']['cell_abundance_heatmap_path'] = combined_heatmap_path
            # Only plot combined cell types
            if 'mod' not in adata.uns:
                adata.uns['mod'] = {}
            adata.uns['mod']['factor_names'] = combined_cell_types
        else:
            print("Using original cell subtypes for visualization (no combination specified in config)")

        temp_config.config = temp_config_dict
        temp_config.get_path = lambda x: temp_config_dict['paths'][x]
        visualizer = Visualizer(temp_config, output_path)

        print("Generating visualizations for merged data...")
        try:
            visualizer._plot_spatial_abundances(adata)
            print("✓ Spatial abundance plots generated successfully")
        except Exception as e:
            print(f"Warning: Could not generate spatial abundance plots: {e}")
        try:
            visualizer._plot_umap(adata)
            print("✓ UMAP plots generated successfully")
        except Exception as e:
            print(f"Warning: Could not generate UMAP plots: {e}")
        try:
            visualizer._plot_batches_by_cell_type(adata)
            print("✓ Batch comparison plots generated successfully")
        except Exception as e:
            print(f"Warning: Could not generate batch comparison plots: {e}")

        print(f"Visualization complete. Results available in: {combined_heatmap_path if combine_cell_types and cell_type_combinations else output_path}")

    def _get_circle_diameter(self, adata):
        """Determines circle diameter based on dataset size."""
        n_cells = adata.n_obs
        if n_cells <= 10000:
            return 6
        elif n_cells <= 15000:
            return 1
        elif n_cells <= 30000:
            return 3
        elif n_cells <= 60000:
            return 1
        else:
            return 0.5

    def plot_training_history(self, model):
        """Plots and saves the model training history."""
        print("Plotting training history...")
        model.plot_history(self.config['quality_control']['spatial_history_skip_epochs'])
        plt.legend(labels=['full data training'])
        plt.savefig(f'{self.run_path}/spatial_training_history.png', dpi=300, bbox_inches='tight')
        plt.close()

    def generate_all_visualizations(self, model, adata_vis, is_batch_mode=False):
        """Generates all QC, spatial, and UMAP plots."""
        print("\nGenerating all visualizations...")
        # self._plot_qc(model, adata_vis)
        
        if not is_batch_mode:
            # self._plot_spatial_abundances(adata_vis)
            self._plot_batches_by_cell_type(adata_vis)
            # self._plot_umap(adata_vis)
        else:
            print("Batch mode: Skipping spatial abundance and UMAP plots.")

    def _plot_qc(self, model, adata_vis):
        """Generates and saves quality control plots for the model."""
        print("Plotting model QC...")
        
        # Check if model has samples (required for QC plots)
        if not hasattr(model, 'samples') or model.samples is None:
            print("Model samples not available (likely loaded from saved model). Skipping QC plots.")
            return
            
        model.plot_QC()
        plt.savefig(f'{self.run_path}/spatial_qc_plots.png', dpi=300, bbox_inches='tight')
        plt.close()

    def _plot_spatial_abundances(self, adata_vis):
        """Plots spatial abundance for all and selected cell types for each slice."""
        print("Plotting spatial cell type abundances...")
        all_cell_types = adata_vis.uns['mod']['factor_names']
        
        for slice_name in self.dataset_config['available_slices']:
            if slice_name not in adata_vis.obs[self.dataset_config['sample_batch_key']].unique():
                print(f"Slice {slice_name} not found in data, skipping.")
                continue
            
            if 'spatial' not in adata_vis.uns or slice_name not in adata_vis.uns['spatial']:
                 print(f"Skipping plotting for {slice_name}, spatial metadata missing in .uns['spatial']")
                 continue

            slice_data = adata_vis[adata_vis.obs[self.dataset_config['sample_batch_key']] == slice_name].copy()
            if slice_data.n_obs == 0: 
                continue

            # Get dynamic circle diameter based on slice size
            circle_diameter = self._get_circle_diameter(slice_data)
            print(f"Using circle_diameter={circle_diameter} for {slice_name} ({slice_data.n_obs} cells)")

            for ct in all_cell_types:
                if ct in slice_data.obs.columns:
                    slice_data.obs[ct] = slice_data.obs[ct].fillna(0)

            # Plot all cell types in multi-page PDF
            pdf_path = f"{self.paths['cell_abundance_heatmap_path']}/{slice_name}_cell_abundance_all.pdf"
            os.makedirs(os.path.dirname(pdf_path), exist_ok=True)
            with PdfPages(pdf_path) as pdf:
                for i in range(0, len(all_cell_types), 4):
                    batch_cell_types = all_cell_types[i:i+4]
                    with mpl.rc_context({'axes.facecolor': 'black', 'figure.figsize': [12, 10]}):
                        sc.pl.spatial(slice_data, color=batch_cell_types, ncols=2, size=1.3, vmin=0,
                                      vmax='p99.2', cmap='magma', show=False, img_key=None, library_id=slice_name, spot_size=30)
                        plt.suptitle(f'Cell type abundance - {slice_name} (Page {i//4 + 1})', fontsize=16)
                        pdf.savefig(bbox_inches='tight', dpi=300)
                        plt.close()
            print(f"Saved all cell type plots to {pdf_path}")

            # Plot selected cell types using cell2location's advanced plotting
            selected_types = [ct for ct in self.config['shared']['selected_cell_types'] if ct in slice_data.obs.columns]
            if selected_types:
                for ct in selected_types:
                    slice_data.obs[ct] = slice_data.obs[ct].fillna(0)
                pdf_path_selected = f"{self.paths['cell_abundance_heatmap_path']}/{slice_name}_selected_cell_types.pdf"
                with PdfPages(pdf_path_selected) as pdf:
                    fig = plot_spatial(adata=slice_data, color=selected_types, labels=selected_types,
                                       show_img=True, style='dark_background', img_alpha=0,
                                       max_color_quantile=0.992, circle_diameter=circle_diameter, colorbar_position='right')
                    pdf.savefig(fig, bbox_inches='tight', dpi=300)
                    plt.close(fig)
                print(f"Saved selected cell type plots to {pdf_path_selected}")
    
    def _plot_umap(self, adata_vis):
        """Generates and saves UMAP plots."""
        if 'umap' in adata_vis.obsm or 'X_umap' in adata_vis.obsm:
            print("Plotting UMAP...")
            with mpl.rc_context({'axes.facecolor': 'white', 'figure.figsize': [12, 6]}):
                plot_colors = ['region_cluster', self.dataset_config['sample_batch_key']]
                sc.pl.umap(adata_vis, color=plot_colors, size=30, ncols=len(plot_colors), show=False)
                plt.savefig(f"{self.paths['cell_abundance_heatmap_path']}/umap_clustering.pdf", bbox_inches='tight')
                plt.close()

    def _plot_batches_by_cell_type(self, adata_vis):
        """Plots all batches for each cell type in separate PDFs with adaptive layouts."""
        print("Plotting all batches for each cell type...")
        all_cell_types = adata_vis.uns['mod']['factor_names']
        available_slices = [slice_name for slice_name in self.dataset_config['available_slices'] 
                           if slice_name in adata_vis.obs[self.dataset_config['sample_batch_key']].unique()]
        
        if len(available_slices) < 2:
            print("Less than 2 slices available, skipping batch comparison plots.")
            return
        
        # Determine subplot layout based on number of batches
        n_batches = len(available_slices)
        if n_batches == 2:
            ncols, nrows = 2, 1
        elif n_batches <= 4:
            ncols, nrows = 2, 2
        elif n_batches <= 6:
            ncols, nrows = 3, 2
        elif n_batches <= 9:
            ncols, nrows = 3, 3
        else:
            ncols, nrows = 4, 3  # For more than 9 batches
        
        print(f"Using {nrows}x{ncols} layout for {n_batches} batches")
        
        # Create output directory
        output_dir = f"{self.paths['cell_abundance_heatmap_path']}/by_cell_type"
        os.makedirs(output_dir, exist_ok=True)
        
        for cell_type in all_cell_types:
            if cell_type not in adata_vis.obs.columns:
                print(f"Cell type {cell_type} not found in data, skipping.")
                continue
                
            pdf_path = f"{output_dir}/{cell_type}_all_batches.pdf"
            
            # Determine global vmin/vmax for consistent color scale
            all_values = []
            valid_slice_data = {}
            
            for slice_name in available_slices:
                slice_data = adata_vis[adata_vis.obs[self.dataset_config['sample_batch_key']] == slice_name]
                if slice_data.n_obs > 0 and cell_type in slice_data.obs.columns and 'spatial' in slice_data.obsm:
                    # Get values and remove NaN values before appending
                    values = slice_data.obs[cell_type].values
                    values_clean = values[~np.isnan(values)]  # Remove NaN values
                    if len(values_clean) > 0:
                        all_values.extend(values_clean)
                        valid_slice_data[slice_name] = slice_data
                    else:
                        print(f"Warning: No valid (non-NaN) values for {cell_type} in {slice_name}")
                else:
                    print(f"Skipping {slice_name} - missing data or spatial coordinates")
            
            if not all_values:
                print(f"No valid data found for {cell_type}, skipping.")
                continue
                
            vmin = 0
            vmax = np.percentile(all_values, 99.2)
            cmap = 'magma'
            print(f"Color scale for {cell_type}: vmin={vmin:.3f}, vmax={vmax:.3f}")

            with PdfPages(pdf_path) as pdf:
                with mpl.rc_context({'axes.facecolor': 'black', 'figure.figsize': [4*ncols, 4*nrows]}):
                    fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 4*nrows))
                    
                    # Handle single subplot case
                    if n_batches == 1:
                        axes = [axes]
                    elif nrows == 1 or ncols == 1:
                        axes = axes.flatten()
                    else:
                        axes = axes.flatten()
                    
                    im_handles = []
                    for i, slice_name in enumerate(available_slices):
                        if slice_name not in valid_slice_data:
                            print(f"Skipping {slice_name} - no valid data")
                            continue
                            
                        if 'spatial' not in adata_vis.uns or slice_name not in adata_vis.uns['spatial']:
                            print(f"Skipping {slice_name}, spatial metadata missing")
                            continue
                            
                        if i >= len(axes):
                            print(f"Skipping {slice_name} - no more axes available")
                            break
                            
                        slice_data = valid_slice_data[slice_name].copy()
                        ax = axes[i]
                        
                        coords = slice_data.obsm.get('spatial', None)
                        if coords is None:
                            print(f"No spatial coordinates for {slice_name}, skipping.")
                            continue
                            
                        # Fill NaN values with 0 before plotting
                        values = slice_data.obs[cell_type].fillna(0).values
                        
                        print(f"Plotting {slice_name}: {len(coords)} points, value range: {values.min():.3f} to {values.max():.3f}")
                        
                        img = ax.scatter(
                            coords[:, 0], coords[:, 1],
                            c=values, cmap=cmap, vmin=vmin, vmax=vmax, s=3
                        )
                        ax.set_title(f'{slice_name}', fontsize=14)
                        ax.set_xticks([])
                        ax.set_yticks([])
                        ax.invert_yaxis()  # Flip the y axis for Visium-style plots
                        im_handles.append(img)
                    
                    # Hide unused subplots
                    for i in range(len(valid_slice_data), len(axes)):
                        axes[i].set_visible(False)

                    plt.suptitle(f'{cell_type} - All Batches', fontsize=16)
                    plt.tight_layout(rect=[0, 0.05, 1, 0.95])

                    # Add a single horizontal colorbar at the bottom
                    if im_handles:
                        cbar = fig.colorbar(
                            im_handles[0], ax=axes, orientation='horizontal',
                            fraction=0.04, pad=0.08, aspect=40
                        )
                        cbar.set_label(f'{cell_type} abundance', fontsize=12)
                        cbar.ax.tick_params(labelsize=10)
                    
                    pdf.savefig(bbox_inches='tight', dpi=300)
                    plt.close(fig)
            
            print(f"Saved {cell_type} batch comparison to {pdf_path}")

class DatasetSplitter:
    """Manages splitting the dataset into batches for distributed processing."""
    def __init__(self, config_manager):
        self.config = config_manager.config
        self.split_config = self.config['dataset']['splitting']
        self.split_dir = self.split_config['split_data_dir']
        self.sample_batch_key = self.config['dataset']['sample_batch_key']

    def split_dataset(self, adata):
        """
        Splits dataset into multiple batches for distributed training,
        preserving sample batch distribution and adding overlapping cells.
        Each split overlaps with the next by sharing the last n cells.
        """
        print("Splitting dataset into batches with stratification and overlap...")
        os.makedirs(self.split_dir, exist_ok=True)

        # Get splitting parameters from the configuration
        batch_size = self.split_config['batch_size']
        overlap_cells = self.split_config['overlap_cells']

        if self.sample_batch_key not in adata.obs.columns:
            raise ValueError(f"Sample batch key '{self.sample_batch_key}' not found in adata.obs")

        # Get unique original batches and calculate the number of output splits
        unique_batches = adata.obs[self.sample_batch_key].unique()
        n_cells = adata.n_obs
        # Adjust calculation to account for overlaps
        effective_batch_size = batch_size - overlap_cells
        n_splits = max(1, int(np.ceil((n_cells - overlap_cells) / effective_batch_size)))

        print(f"Splitting {n_cells} cells into {n_splits} splits of max {batch_size} cells each.")
        print(f"Overlap cells between splits: {overlap_cells}")
        print(f"Effective batch size (excluding overlap): {effective_batch_size}")

        # Group cell indices by their original sample batch
        batch_indices = {}
        for batch_name in unique_batches:
            batch_indices[batch_name] = np.where(adata.obs[self.sample_batch_key] == batch_name)[0]

        # Create a global ordering of all cell indices while preserving batch proportions
        all_indices = []
        max_batch_size = max(len(indices) for indices in batch_indices.values())
        
        # Interleave indices from different batches to maintain proportions
        for i in range(max_batch_size):
            for batch_name in unique_batches:
                if i < len(batch_indices[batch_name]):
                    all_indices.append(batch_indices[batch_name][i])

        split_info = []
        current_start = 0
        
        for split_idx in range(n_splits):
            # Calculate the range of indices for this split
            if split_idx == 0:
                # First split: no overlap at the beginning
                split_start = current_start
                split_end = min(current_start + batch_size, len(all_indices))
            else:
                # Subsequent splits: start with overlap from previous split
                split_start = max(0, current_start - overlap_cells)
                split_end = min(current_start + effective_batch_size, len(all_indices))
            
            # Get the actual cell indices for this split
            split_cell_indices = all_indices[split_start:split_end]
            
            if not split_cell_indices:
                break
                
            # Create and save the AnnData object for the split
            split_adata = adata[split_cell_indices].copy()
            split_filename = f"{self.split_dir}/batch_{split_idx}.h5ad"
            split_adata.write(split_filename)
            
            # Record metadata about the split
            split_batch_counts = split_adata.obs[self.sample_batch_key].value_counts()
            
            # Calculate overlap information
            overlap_start = 0
            overlap_end = 0
            if split_idx > 0:
                overlap_start = overlap_cells
            if split_idx < n_splits - 1 and split_end < len(all_indices):
                overlap_end = min(overlap_cells, len(split_cell_indices))
            
            info = {
                'batch_idx': split_idx,
                'n_cells': len(split_cell_indices),
                'filename': split_filename,
                'cell_indices': [int(i) for i in split_cell_indices],
                'batch_distribution': {k: int(v) for k, v in split_batch_counts.items()},
                'overlap_with_previous': overlap_start if split_idx > 0 else 0,
                'overlap_with_next': overlap_end if split_idx < n_splits - 1 else 0,
                'global_start_idx': split_start,
                'global_end_idx': split_end
            }
            split_info.append(info)
            
            print(f"Saved split {split_idx} ({len(split_cell_indices)} cells) to {split_filename}")
            print(f"  Split {split_idx} batch distribution: {dict(split_batch_counts)}")
            if split_idx > 0:
                print(f"  Overlap with previous split: {overlap_start} cells")
            if split_idx < n_splits - 1 and split_end < len(all_indices):
                print(f"  Will overlap with next split: {overlap_end} cells")
            
            # Move to the next split position
            current_start += effective_batch_size
        
        # Save the metadata for all splits to a JSON file
        split_info_file = f"{self.split_dir}/batch_info.json"
        with open(split_info_file, 'w') as f:
            json.dump(split_info, f, indent=2)
        
        print(f"Split information saved to {split_info_file}")
        print("Dataset splitting complete.")
        
        # Print overlap summary
        print("\nOverlap Summary:")
        for i, info in enumerate(split_info):
            if i > 0:
                prev_info = split_info[i-1]
                actual_overlap = len(set(info['cell_indices']) & set(prev_info['cell_indices']))
                print(f"  Split {i-1} -> Split {i}: {actual_overlap} overlapping cells")
        
        return split_info

    def load_batch(self, batch_number):
        """Loads a specific pre-split batch file."""
        batch_info_file = f"{self.split_dir}/batch_info.json"
        if not os.path.exists(batch_info_file):
            raise FileNotFoundError(f"Batch info file not found: {batch_info_file}. Run with --split-only first.")
        
        with open(batch_info_file, 'r') as f:
            all_batches_info = json.load(f)
        
        if batch_number >= len(all_batches_info):
            raise ValueError(f"Batch {batch_number} not found. Available batches: 0-{len(all_batches_info)-1}")
        
        batch_file = all_batches_info[batch_number]['filename']
        print(f"Loading batch {batch_number} from {batch_file}")
        return sc.read_h5ad(batch_file)


class BatchMerger:
    """Handles merging results from multiple trained batches."""
    def __init__(self, config_manager):
        self.config = config_manager.config
        self.dataset_config = self.config['dataset']
        self.split_config = self.dataset_config['splitting']
        self.spatial_output_path = config_manager.get_path('spatial_output')
        
    def merge_all_batches(self):
        """Merges all trained batch results into a single unified dataset."""
        print("=== Merging All Trained Batches ===")
        
        # Find all batch directories
        batch_dirs = self._find_batch_directories()
        if not batch_dirs:
            raise ValueError("No trained batch directories found")
        
        print(f"Found {len(batch_dirs)} batch directories to merge")
        
        # Load batch metadata to understand overlaps
        batch_metadata = self._load_batch_metadata()
        
        # Load and combine all batch results with proper overlap handling
        merged_adata = self._load_and_combine_batches(batch_dirs, batch_metadata)
        
        # Save merged results
        merged_output_path = f"{self.spatial_output_path}/merged_results"
        os.makedirs(merged_output_path, exist_ok=True)
        
        self._save_merged_results(merged_adata, merged_output_path)
        
        # Generate visualizations for merged data
        self._generate_merged_visualizations(merged_adata, merged_output_path)
        
        print(f"Merged results saved to: {merged_output_path}")
        
    def _load_batch_metadata(self):
        """Loads the batch metadata to understand overlaps."""
        split_dir = self.split_config['split_data_dir']
        batch_info_file = f"{split_dir}/batch_info.json"
        
        if not os.path.exists(batch_info_file):
            print("Warning: No batch metadata found. Will use simple duplicate removal.")
            return None
            
        with open(batch_info_file, 'r') as f:
            batch_metadata = json.load(f)
        
        print(f"Loaded metadata for {len(batch_metadata)} batches")
        return batch_metadata
        
    def _find_batch_directories(self):
        """Finds all batch_* directories with completed results."""
        batch_dirs = []
        
        if not os.path.exists(self.spatial_output_path):
            return batch_dirs
            
        for item in os.listdir(self.spatial_output_path):
            if item.startswith('batch_'):
                batch_path = os.path.join(self.spatial_output_path, item)
                # Check if this batch has completed processing
                final_increments_path = f"{batch_path}/sp_final_increments.pkl"
                if os.path.exists(final_increments_path):
                    batch_dirs.append(batch_path)
                else:
                    print(f"Warning: Batch {item} appears incomplete (missing final increments)")
        
        return sorted(batch_dirs)
    
    def _load_and_combine_batches(self, batch_dirs, batch_metadata):
        """Loads and combines data from all batch directories with memory-efficient incremental merging."""
        print("Loading batch results with memory-efficient incremental merging...")
        
        if not batch_dirs:
            raise ValueError("No valid batch data found to merge")
        
        # Load the first batch as the base and extract spatial information
        first_batch_dir = batch_dirs[0]
        batch_name = os.path.basename(first_batch_dir)
        print(f"Loading base batch: {batch_name}...")
        
        # Use load_incremental_adata for robust loading
        adata_base = load_incremental_adata(first_batch_dir)
        if adata_base is None:
            # fallback to manual restoration if needed
            final_increments_path = f"{first_batch_dir}/sp_final_increments.pkl"
            base_adata_path = f"{first_batch_dir}/sp.h5ad"
            model_increments_path = f"{first_batch_dir}/sp_model_increments.pkl"
            if not all(os.path.exists(p) for p in [base_adata_path, model_increments_path, final_increments_path]):
                raise ValueError(f"Missing required files for base batch {batch_name}")
            adata_base = sc.read_h5ad(base_adata_path)
            adata_with_posterior = restore_adata_additions(adata_base, model_increments_path)
            adata_base = restore_adata_additions(adata_with_posterior, final_increments_path)
        
        print(f"Base batch loaded: {adata_base.shape}")
        
        # Create lean version of first batch (without spatial info)
        merged_adata = self._create_lean_adata(adata_base)
        
        # Store other metadata from the first batch
        factor_names = merged_adata.uns.get('mod', {}).get('factor_names', None)
        
        # Clean up first batch from memory
        del adata_base
        import gc
        gc.collect()
        
        # Incrementally merge each subsequent batch (lean versions)
        for i, batch_dir in enumerate(batch_dirs[1:], 1):
            batch_name = os.path.basename(batch_dir)
            print(f"Loading and merging batch {i+1}/{len(batch_dirs)}: {batch_name}...")
            
            # Use load_incremental_adata for robust loading
            current_batch_full = load_incremental_adata(batch_dir)
            if current_batch_full is None:
                final_increments_path = f"{batch_dir}/sp_final_increments.pkl"
                base_adata_path = f"{batch_dir}/sp.h5ad"
                model_increments_path = f"{batch_dir}/sp_model_increments.pkl"
                if not all(os.path.exists(p) for p in [base_adata_path, model_increments_path, final_increments_path]):
                    print(f"Warning: Skipping {batch_name} - missing required files")
                    continue
                try:
                    adata_base = sc.read_h5ad(base_adata_path)
                    adata_with_posterior = restore_adata_additions(adata_base, model_increments_path)
                    current_batch_full = restore_adata_additions(adata_with_posterior, final_increments_path)
                except Exception as e:
                    print(f"Error processing batch {batch_name}: {e}")
                    continue

            print(f"  Current batch shape: {current_batch_full.shape}")
            
            # Create lean version (without spatial info)
            current_batch_lean = self._create_lean_adata(current_batch_full)
            
            # Use batch_metadata to remove overlap cells from current_batch_lean
            if batch_metadata is not None:
                # batch_metadata[i] corresponds to batch_dirs[i]
                # batch_dirs[0] -> batch_metadata[0], batch_dirs[1] -> batch_metadata[1], etc.
                # For batch i (current), remove overlap with previous batch (i.e., overlapping cells at the start)
                curr_meta = batch_metadata[i]
                overlap_with_previous = curr_meta.get('overlap_with_previous', 0)
                if overlap_with_previous > 0:
                    # Remove the first 'overlap_with_previous' cells from current_batch_lean
                    # Use the order in 'cell_indices' from batch_metadata
                    cell_indices = curr_meta.get('cell_indices', [])
                    if len(cell_indices) == current_batch_lean.n_obs:
                        # Remove the first N cells
                        # keep_indices = cell_indices[overlap_with_previous:]
                        keep_obs_names = [current_batch_lean.obs.index[j] for j in range(overlap_with_previous, len(cell_indices))]
                        current_batch_lean = current_batch_lean[keep_obs_names].copy()
                        print(f"  Removed {overlap_with_previous} overlapping cells from start of batch {i} ({batch_name})")
                    else:
                        print(f"  Warning: cell_indices length does not match n_obs for batch {i}, skipping overlap removal")
                else:
                    print(f"  No overlap to remove for batch {i} ({batch_name})")
            else:
                print("  No batch_metadata provided, skipping overlap removal.")

            print(f"  Pre-merge shape: {merged_adata.shape}")
            
            # Concatenate lean versions
            merged_adata = ad.concat([merged_adata, current_batch_lean], join='outer', index_unique=None)
            
            print(f"  Post-merge shape: {merged_adata.shape}")
            
            # Clear the current batch from memory
            del current_batch_full, current_batch_lean
            gc.collect()
        
        # Remove duplicates (simple approach - no need for complex metadata handling)
        if merged_adata.obs.index.duplicated().any():
            n_duplicates = merged_adata.obs.index.duplicated().sum()
            print(f"Removing {n_duplicates} duplicate observations (keeping first occurrence)")
            merged_adata = merged_adata[~merged_adata.obs.index.duplicated(keep='first')].copy()
        else:
            print("No duplicate observations detected in merged dataset")
        
        # load original base adata to restore additional metadata
        adata_origin = sc.read_h5ad(self.config['paths']['spatial_data'])
        merge_origin_adata = self._merge_obs_data(adata_origin, merged_adata)
        merge_origin_adata = self._merge_obsm_data(merge_origin_adata, merged_adata)
        
        # Restore factor names in uns if they were available
        if factor_names is not None:
            if 'mod' not in merge_origin_adata.uns:
                merge_origin_adata.uns['mod'] = {}
            merge_origin_adata.uns['mod']['factor_names'] = factor_names
        
        print(f"Final merged dataset shape: {merge_origin_adata.shape}")
        
        # Final garbage collection
        gc.collect()
        
        return merge_origin_adata

    def _merge_obs_data(self, target_adata, source_adata, exclude_cols=['region_cluster', 'leiden']):
        print("--- Starting .obs metadata merge ---")
        print(f"Original target shape: {target_adata.n_obs} obs, {target_adata.n_vars} vars")
        print(f"Source shape: {source_adata.n_obs} obs, {source_adata.n_vars} vars")

        # Use a copy of source .obs to avoid modifying the original object
        source_obs_df = source_adata.obs.copy()
        
        # Drop the excluded columns from the source dataframe if they exist
        cols_to_drop = [col for col in exclude_cols if col in source_obs_df.columns]
        if cols_to_drop:
            source_obs_df.drop(columns=cols_to_drop, inplace=True)
            print(f"\nExcluding columns from source: {cols_to_drop}")

        # Identify columns in the source that are not in the target
        original_target_cols = target_adata.obs.columns
        source_cols_to_add = source_obs_df.columns.difference(original_target_cols)
        
        if len(source_cols_to_add) == 0:
            print("\nNo new columns to add. Aborting .obs merge.")
            return target_adata

        print(f"Found {len(source_cols_to_add)} new columns to merge.")
        print("Example columns to be added:", list(source_cols_to_add[:5]))

        # Use a left join to merge the dataframes on their index
        merged_obs = target_adata.obs.join(source_obs_df[source_cols_to_add], how='left')
        target_adata.obs = merged_obs

        print(f"\n.obs merge complete. New .obs shape: {target_adata.obs.shape}")
        nan_counts = target_adata.obs[source_cols_to_add].isnull().sum().sum()
        print(f"Total NaN values in new .obs columns: {nan_counts} (expected for non-matching cells).")
        print("--- Finished .obs merge ---\n")
        
        return target_adata

    def _merge_obsm_data(self, target_adata, source_adata):
        """
        Merges .obsm data from a source to a target AnnData object.

        It aligns the arrays in .obsm based on the cell indices. For cells in the
        target but not the source, the corresponding rows in the new .obsm arrays
        will be filled with NaN.

        Args:
            target_adata (sc.AnnData): The AnnData object to modify in place.
            source_adata (sc.AnnData): The AnnData object providing the .obsm data.

        Returns:
            sc.AnnData: The modified target_adata object.
        """
        print("--- Starting .obsm data merge ---")
        
        source_keys = list(source_adata.obsm.keys())
        if not source_keys:
            print("Source AnnData has no .obsm data to merge. Skipping.")
            return target_adata
            
        print(f"Found .obsm keys in source: {source_keys}")

        # Get the cell indices for alignment
        target_index = target_adata.obs.index
        source_index = source_adata.obs.index

        for key in source_keys:
            if key == 'X_umap':
                print("Skipping 'X_umap' key as it is not needed for merging.")
                continue
            if key in target_adata.obsm:
                print(f"Key '{key}' already exists in target.obsm. Skipping to avoid overwrite.")
                continue
                
            print(f"Processing key: '{key}'")
            # Get the data from the source .obsm
            source_array = source_adata.obsm[key]
            
            # If already a DataFrame, use as is; otherwise, convert to DataFrame with original column names if possible
            if isinstance(source_array, pd.DataFrame):
                source_df = source_array
            else:
                # Try to get column names from .columns or .dtype.names
                colnames = None
                if hasattr(source_array, 'columns'):
                    colnames = list(source_array.columns)
                elif hasattr(source_array, 'dtype') and getattr(source_array.dtype, 'names', None):
                    colnames = list(source_array.dtype.names)
                elif hasattr(source_array, 'shape') and len(source_array.shape) == 2:
                    colnames = [f"{key}_{i}" for i in range(source_array.shape[1])]
                source_df = pd.DataFrame(source_array, index=source_index, columns=colnames)
            
            # Re-index the DataFrame based on the target's cell index.
            # This aligns the data and fills missing rows with NaN.
            reindexed_df = source_df.reindex(target_index)
            
            # Assign to target AnnData
            target_adata.obsm[key] = reindexed_df
            
            print(f"Successfully merged '{key}'. New shape: {target_adata.obsm[key].shape}")

        print("--- Finished .obsm merge ---")
        return target_adata

    def _create_lean_adata(self, adata):
        """Creates a lean version of AnnData without spatial information."""
        # Create a copy
        lean_adata = adata.copy()
        
        # Remove spatial information from uns
        if 'spatial' in lean_adata.uns:
            del lean_adata.uns['spatial']
        
        # Remove spatial coordinates from obsm
        removal_obsm_keys = [key for key in lean_adata.obsm.keys() if 'spatial' in key.lower() or '_umap' in key.lower()]
        for key in removal_obsm_keys:
            del lean_adata.obsm[key]

        # Remove any spatial-related keys from obsp
        removal_obs_keys = [key for key in lean_adata.obs_keys() if 'leiden' in key.lower() or 'region_cluster' in key.lower() or 'orig.ident' in key.lower()]
        for key in removal_obs_keys:
            del lean_adata.obs[key]

        return lean_adata

    def _save_merged_results(self, merged_adata, output_path):
        """Saves the merged results."""
        print("Saving merged results...")
        
        # Save the full merged adata
        merged_adata.write_h5ad(f"{output_path}/merged_adata.h5ad")
        
        # Save cell abundances and clusters to CSV
        factor_names = merged_adata.uns['mod']['factor_names']
        columns_to_export = factor_names + ["region_cluster", self.dataset_config['sample_batch_key']]
        
        # Only export columns that exist
        existing_columns = [col for col in columns_to_export if col in merged_adata.obs.columns]

        csv_path = f"{self.config['paths']['cell_abundance_results']}/cell_abundances_and_clusters.csv"
        merged_adata.obs[existing_columns].to_csv(csv_path)
        print(f"Saved merged cell abundances to {csv_path}")
        
    def _generate_merged_visualizations(self, merged_adata, output_path):
        """Generates visualizations for the merged dataset using shared logic."""
        config = self.config
        combine_cell_types = config.get('shared', {}).get('combine_cell_types', False)
        cell_type_combinations = config.get('shared', {}).get('cell_type_combinations', {})
        Visualizer.generate_visualizations_for_adata(
            merged_adata, config, output_path, combine_cell_types, cell_type_combinations
        )


class Cell2LocationPipeline:
    """Orchestrates the entire Cell2Location workflow."""
    def __init__(self):
        self.args = self._parse_args()
        self.config_manager = ConfigManager()
        self.data_processor = DataProcessor(self.config_manager)

    def _parse_args(self):
        parser = argparse.ArgumentParser(description='Cell2Location Spatial Mapping Pipeline')
        parser.add_argument('--batch', type=int, default=None, help='Batch number to process')
        parser.add_argument('--split-only', action='store_true', help='Only split the dataset')
        parser.add_argument('--merge', action='store_true', help='Merge all trained batches')
        parser.add_argument('--visualize', action='store_true', help='Generate visualizations for existing results')
        return parser.parse_args()

    def run(self):
        """Main entry point for the pipeline."""
        if self.args.visualize:
            self._visualize()
        elif self.args.merge:
            self._merge_batches()
        elif self.args.split_only:
            self._run_split_only()
        elif self.args.batch is not None:
            self._run_batch(self.args.batch)
        else:
            self._run_full()
    
    def _visualize(self):
        """Generates visualizations for existing results using shared logic."""
        print("=== Mode: Visualize Existing Merged Results ===")
        merged_output_path = f"{self.config_manager.get_path('spatial_output')}/merged_results"
        merged_adata_path = f"{merged_output_path}/merged_adata.h5ad"
        fallback_adata_path = f"{self.config_manager.get_path('spatial_output')}/sp.h5ad"

        # Try merged results first, then fallback to single-run output
        if os.path.exists(merged_adata_path):
            print(f"Loading merged dataset from {merged_adata_path}")
            merged_adata = load_processed_spatial_data(merged_output_path)
            output_path = merged_output_path
        else:
            if not os.path.exists(fallback_adata_path):
                raise FileNotFoundError(
                    f"No merged results found at {merged_adata_path} and no single-run results found at {fallback_adata_path}."
                    " Run with --merge or process the full dataset first."
                )
            print(f"No merged results found. Falling back to single-run output at {fallback_adata_path}")
            merged_adata = load_incremental_adata(self.config_manager.get_path('spatial_output'))
            if merged_adata is None:
                merged_adata = load_processed_spatial_data(self.config_manager.get_path('spatial_output'))
            output_path = self.config_manager.get_path('spatial_output')

        config = self.config_manager.config
        combine_cell_types = config.get('shared', {}).get('combine_cell_types', False)
        cell_type_combinations = config.get('shared', {}).get('cell_type_combinations', {})

        Visualizer.generate_visualizations_for_adata(
            merged_adata, config, output_path, combine_cell_types, cell_type_combinations
        )
        return

    def _merge_batches(self):
        """Merges all trained batch results."""
        merger = BatchMerger(self.config_manager)
        merger.merge_all_batches()

    def _run_split_only(self):
        """Handles the dataset splitting workflow."""
        print("=== Mode: Dataset Splitting Only ===")
        splitter = DatasetSplitter(self.config_manager)
        initial_adata = self.data_processor.load_initial_spatial_data()
        splitter.split_dataset(initial_adata)

    def _run_batch(self, batch_number):
        """Runs the complete pipeline for a single batch."""
        print(f"=== Mode: Batch Processing (Batch {batch_number}) ===")
        run_path = f"{self.config_manager.get_path('spatial_output')}/batch_{batch_number}"
        os.makedirs(run_path, exist_ok=True)
        
        splitter = DatasetSplitter(self.config_manager)
        initial_adata = splitter.load_batch(batch_number)
        
        self._execute_workflow(initial_adata, run_path)

    def _run_full(self):
        """Runs the complete pipeline on the entire dataset."""
        print("=== Mode: Full Dataset Processing ===")
        run_path = self.config_manager.get_path('spatial_output')
        os.makedirs(run_path, exist_ok=True)
        
        initial_adata = self.data_processor.load_initial_spatial_data()
        self._execute_workflow(initial_adata, run_path)
    
    def _execute_workflow(self, adata_initial, run_path):
        """Executes the core data processing, model training, and visualization steps."""
        visualizer = Visualizer(self.config_manager, run_path)
        model_manager = ModelManager(self.config_manager, visualizer, run_path)
        is_batch_mode = self.args.batch is not None
        
        # Define paths for intermediate and final files
        final_increments_path = f"{run_path}/sp_final_increments.pkl"
        base_adata_path = f"{run_path}/sp.h5ad"
        model_increments_path = f"{run_path}/sp_model_increments.pkl"
        model_file_path = f"{run_path}/model.pt"

        # Check if final results already exist to skip to visualization
        adata_final = load_incremental_adata(run_path)
        if adata_final is not None:
            print("Found final processed data. Loading from increments and visualizing.")
            model = cell2location.models.Cell2location.load(run_path, adata_final)
            visualizer.generate_all_visualizations(model, adata_final, is_batch_mode)
            return

        # Step 1: Prepare data
        adata_prepared, inf_aver = self.data_processor.load_reference_and_prepare_spatial(adata_initial)
        
        # Step 2: Train or load model
        model, adata_with_posterior = model_manager.train_or_load_model(
            adata_prepared, inf_aver, model_file_path, base_adata_path, model_increments_path
        )
        
        # Step 3: Post-processing (clustering, etc.)
        adata_final = model_manager.process_and_save_results(adata_with_posterior, final_increments_path)
        
        # Step 4: Visualization
        visualizer.generate_all_visualizations(model, adata_final, is_batch_mode)
        
        print(f"\nResults saved in: {run_path}")
        print("Next steps:")
        print("- Run pred_cell2location_3_DGE_after_deconv.py for DGE analysis")
        print("- Run pred_cell2location_4_NMF_analysis.py for NMF colocalization analysis")


if __name__ == '__main__':
    pipeline = Cell2LocationPipeline()
    pipeline.run()