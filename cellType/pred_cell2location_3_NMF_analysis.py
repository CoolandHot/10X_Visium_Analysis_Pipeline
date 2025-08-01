import scanpy as sc
import pandas as pd
import os
import matplotlib.pyplot as plt
import pickle
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from cell2location import run_colocation
import squidpy as sq


# Import utility functions from DGE script
from cellType.pred_cell2location_utils import (
    load_configuration, 
    load_processed_spatial_data, 
    extract_adata_additions, 
    restore_adata_additions,
)

def _validate_and_clean_nmf_input(adata, config):
    """
    Validates and cleans the input data for NMF analysis.
    
    Parameters:
    -----------
    adata : AnnData
        Spatial AnnData object with cell type abundances
    config : dict
        Configuration dictionary
        
    Returns:
    --------
    AnnData
        Cleaned AnnData object ready for NMF
    """
    print("\n--- Validating and Cleaning NMF Input Data ---")
    
    # Check for the expected abundance matrix
    abundance_key = 'q05_cell_abundance_w_sf'
    if abundance_key not in adata.obsm_keys():
        print(f"Warning: '{abundance_key}' not found in adata.obsm.")
        print(f"Available obsm keys: {list(adata.obsm.keys())}")
        # Try to find alternative abundance matrices
        abundance_keys = [k for k in adata.obsm.keys() if 'abundance' in k.lower()]
        if abundance_keys:
            abundance_key = abundance_keys[0]
            print(f"Using alternative abundance key: {abundance_key}")
        else:
            raise ValueError("No cell abundance matrix found in adata.obsm")
    
    abundance_matrix = adata.obsm[abundance_key]
    original_is_dataframe = isinstance(abundance_matrix, pd.DataFrame)
    
    print(f"Abundance matrix shape: {abundance_matrix.shape}")
    print(f"Abundance matrix type: {type(abundance_matrix)}")
    print(f"Original is DataFrame: {original_is_dataframe}")
    
    # Store column names if it's a DataFrame
    if original_is_dataframe:
        column_names = abundance_matrix.columns.tolist()
        abundance_matrix_values = abundance_matrix.values
    else:
        column_names = None
        # Convert to dense array if sparse
        if hasattr(abundance_matrix, 'toarray'):
            abundance_matrix_values = abundance_matrix.toarray()
        else:
            abundance_matrix_values = np.asarray(abundance_matrix)
    
    # Ensure we have a numpy array for processing
    abundance_matrix_values = np.asarray(abundance_matrix_values)
    
    # Check for NaN values
    nan_count = np.isnan(abundance_matrix_values).sum()
    inf_count = np.isinf(abundance_matrix_values).sum()
    negative_count = (abundance_matrix_values < 0).sum()
    
    print("Data quality check:")
    print(f"  NaN values: {nan_count}")
    print(f"  Inf values: {inf_count}")
    print(f"  Negative values: {negative_count}")
    
    # Handle potential issues with min/max on arrays with NaN
    if nan_count > 0 or inf_count > 0:
        print(f"  Min value: {np.nanmin(abundance_matrix_values)}")
        print(f"  Max value: {np.nanmax(abundance_matrix_values)}")
        print(f"  Mean value: {np.nanmean(abundance_matrix_values)}")
    else:
        print(f"  Min value: {abundance_matrix_values.min()}")
        print(f"  Max value: {abundance_matrix_values.max()}")
        print(f"  Mean value: {abundance_matrix_values.mean()}")
    
    # Create a copy of adata to work with
    adata_clean = adata.copy()
    
    # Handle NaN and Inf values by removing affected cells
    if nan_count > 0 or inf_count > 0:
        print("Removing cells with NaN/Inf values...")
        
        # Find rows (cells) with NaN or Inf values
        has_nan = np.isnan(abundance_matrix_values).any(axis=1)
        has_inf = np.isinf(abundance_matrix_values).any(axis=1)
        bad_cells = has_nan | has_inf
        
        n_bad_cells = np.sum(bad_cells)
        print(f"  Cells with NaN: {np.sum(has_nan)}")
        print(f"  Cells with Inf: {np.sum(has_inf)}")
        print(f"  Total cells to remove: {n_bad_cells}")
        
        if n_bad_cells > 0:
            # Keep only good cells
            good_cells = ~bad_cells
            adata_clean = adata_clean[good_cells].copy()
            abundance_matrix_values = abundance_matrix_values[good_cells]
            
            print(f"  Removed {n_bad_cells} cells")
            print(f"  Remaining cells: {adata_clean.n_obs}")
            
            # Verify no NaN/Inf remain
            print("After removal:")
            print(f"  NaN values: {np.isnan(abundance_matrix_values).sum()}")
            print(f"  Inf values: {np.isinf(abundance_matrix_values).sum()}")
    
    # Handle negative values (set to 0 for abundance data)
    negative_count_remaining = (abundance_matrix_values < 0).sum()
    if negative_count_remaining > 0:
        print(f"Setting {negative_count_remaining} negative values to 0...")
        abundance_matrix_values = np.maximum(abundance_matrix_values, 0)
    
    # Check if all values are zero (would cause NMF to fail)
    row_sums = abundance_matrix_values.sum(axis=1)
    col_sums = abundance_matrix_values.sum(axis=0)
    zero_rows = np.sum(row_sums == 0)
    zero_cols = np.sum(col_sums == 0)
    
    if zero_rows > 0:
        print(f"Warning: {zero_rows} observations have all-zero abundances")
        # Remove cells with all-zero abundances
        non_zero_cells = row_sums > 0
        if np.sum(non_zero_cells) < adata_clean.n_obs:
            print(f"Removing {np.sum(~non_zero_cells)} cells with all-zero abundances")
            adata_clean = adata_clean[non_zero_cells].copy()
            abundance_matrix_values = abundance_matrix_values[non_zero_cells]
            print(f"Remaining cells after zero removal: {adata_clean.n_obs}")
    
    if zero_cols > 0:
        print(f"Warning: {zero_cols} cell types have all-zero abundances")
    
    # Add small epsilon only if the entire matrix is zero
    if abundance_matrix_values.max() == 0:
        print("All values are zero! Adding small epsilon...")
        abundance_matrix_values += 1e-10
    
    # Convert back to DataFrame if original was DataFrame, otherwise keep as numpy array
    if original_is_dataframe and column_names is not None:
        # Create DataFrame with proper column names and cell indices
        abundance_matrix_clean = pd.DataFrame(
            abundance_matrix_values,
            index=adata_clean.obs_names,
            columns=column_names
        )
        print("Converted cleaned matrix back to DataFrame for cell2location compatibility")
    else:
        abundance_matrix_clean = abundance_matrix_values
        print("Keeping as numpy array")
    
    # Update the abundance matrix in adata_clean
    adata_clean.obsm[abundance_key] = abundance_matrix_clean
    
    print(f"Final data shape: {abundance_matrix_clean.shape}")
    print(f"Final data type: {type(abundance_matrix_clean)}")
    if hasattr(abundance_matrix_clean, 'values'):
        print(f"Final data range: [{abundance_matrix_clean.values.min():.6f}, {abundance_matrix_clean.values.max():.6f}]")
    else:
        print(f"Final data range: [{abundance_matrix_clean.min():.6f}, {abundance_matrix_clean.max():.6f}]")
    print("--- Data validation and cleaning completed ---\n")
    
    return adata_clean

def run_nmf_colocalization(adata, config, n_fact_range_list):
    """
    Runs NMF colocalization analysis.
    
    Parameters:
    -----------
    adata : AnnData
        Spatial AnnData object with cell type abundances
    config : dict
        Configuration dictionary
    n_fact_range_list : list
        List of factor numbers to try for NMF
        
    Returns:
    --------
    tuple
        (results_dict, updated_adata) or (None, adata) if failed
    """
    print("\n--- Running NMF Colocalization Analysis ---")
    
    nmf_config = config['nmf_analysis']
    dataset_config = config['dataset']
    
    # Ensure the sample column exists
    sample_col = dataset_config['sample_batch_key']
    if sample_col not in adata.obs.columns:
        raise ValueError(f"Sample column '{sample_col}' not found in adata.obs. Please check config/cellType_config.yaml.")
    
    # Validate and clean input data
    try:
        adata_clean = _validate_and_clean_nmf_input(adata, config)
        print(f"Data cleaned successfully. Shape: {adata_clean.shape}")
        print(f"Sample column '{sample_col}' still present: {sample_col in adata_clean.obs.columns}")
    except Exception as e:
        print(f"Error during data validation/cleaning: {e}")
        import traceback
        traceback.print_exc()
        return None, adata
    
    nmf_export_path = config['paths']['nmf_analysis_path']
    os.makedirs(nmf_export_path, exist_ok=True)
    nmf_export_args = {'path': nmf_export_path}

    print(f"NMF parameters: n_fact={n_fact_range_list}, restarts={nmf_config['n_restarts']}")
    print(f"Using abundance matrix: {list(adata_clean.obsm.keys())}")

    try:
        res_dict, adata_updated = run_colocation(
            adata_clean,
            model_name='CoLocatedGroupsSklearnNMF',
            train_args={
                'n_fact': n_fact_range_list,
                'sample_name_col': sample_col,
                'n_restarts': nmf_config['n_restarts']
            },
            model_kwargs={
                'alpha': nmf_config['alpha'],
                'init': nmf_config['init'],
                "nmf_kwd_args": {"tol": nmf_config['tol']}
            },
            export_args=nmf_export_args,
        )
        print("NMF analysis completed successfully.")
        return res_dict, adata_updated
    
    except Exception as e:
        print(f"Error during NMF colocalization: {e}")
        import traceback
        traceback.print_exc()
        print("Skipping NMF analysis.")
        return None, adata

def analyze_nmf_results(adata, nmf_results_dict, config):
    """
    Analyze and visualize NMF results.
    
    Parameters:
    -----------
    adata : AnnData
        AnnData object with NMF results
    nmf_results_dict : dict
        Results dictionary from NMF analysis
    config : dict
        Configuration dictionary
    run_name_path : str
        Path to save analysis results
    """
    print("\n--- Analyzing NMF Results ---")
    
    # Create output directory for NMF analysis
    nmf_analysis_path = config['paths']['nmf_analysis_path']
    os.makedirs(nmf_analysis_path, exist_ok=True)
    
    # Look for NMF group columns in adata.obs
    nmf_cols = [col for col in adata.obs.columns if 'mean_nUMI_factorsfact' in col]
    
    if not nmf_cols:
        print("No NMF group columns found in adata.obs. Analysis incomplete.")
        return
    
    print(f"Found NMF columns: {nmf_cols}")
    
    # Generate summary plots
    _generate_nmf_summary_plots(adata, nmf_cols, nmf_analysis_path, config)
    
    # Save NMF group assignments to CSV
    _save_nmf_results_to_csv(adata, nmf_cols, nmf_analysis_path, config)

def _generate_nmf_summary_plots(adata, nmf_cols, output_path, config):
    """
    Generate summary plots for NMF results.
    """
    try:
        sample_batch_key = config['dataset']['sample_batch_key']
        available_slices = config['dataset']['available_slices']
        
        # Plot NMF groups spatially for each slice
        has_spatial = 'spatial' in adata.uns and adata.uns['spatial']
        
        if has_spatial:
            for slice_name in available_slices:
                # Use proper pandas syntax to check values
                if slice_name in adata.obs[sample_batch_key].unique():
                    slice_data = adata[adata.obs[sample_batch_key] == slice_name].copy()
                    slice_data.uns['spatial'] = {
                        k: v for k, v in slice_data.uns['spatial'].items() if k == slice_name
                    }
                    
                    if slice_data.n_obs > 0:
                        # Plot NMF groups for this slice
                        pdf_filename = f'{output_path}/nmf_spatial_{slice_name}.pdf'
                        with PdfPages(pdf_filename) as pdf:
                            for nmf_col in nmf_cols:
                                if nmf_col in slice_data.obs.columns:
                                    try:
                                        fig, ax = plt.subplots(1, 1, figsize=(10, 8))
                                        sq.pl.spatial_scatter(
                                            slice_data,
                                            color=nmf_col,
                                            title=f'{slice_name}: {nmf_col}',
                                            ax=ax,
                                            legend_loc='right margin'
                                        )
                                        pdf.savefig(fig, bbox_inches='tight')
                                        plt.close(fig)
                                    except Exception as e:
                                        print(f"Could not plot {nmf_col} for {slice_name}: {e}")
                        
                        print(f"Saved NMF spatial plots for {slice_name} to {pdf_filename}")
        
        # Generate UMAP plots colored by NMF groups
        if 'X_umap' in adata.obsm:
            pdf_filename = f'{output_path}/nmf_umap_analysis.pdf'
            with PdfPages(pdf_filename) as pdf:
                for nmf_col in nmf_cols:
                    try:
                        sc.pl.umap(adata, color=nmf_col, title=f'UMAP: {nmf_col}', show=False)
                        pdf.savefig(plt.gcf(), bbox_inches='tight')
                        plt.close()
                    except Exception as e:
                        print(f"Could not generate UMAP plot for {nmf_col}: {e}")
            
            print(f"Saved NMF UMAP plots to {pdf_filename}")
    
    except Exception as e:
        print(f"Error generating NMF summary plots: {e}")

def _save_nmf_results_to_csv(adata, nmf_cols, output_path, config):
    """
    Save NMF group assignments to CSV files.
    """
    try:
        sample_batch_key = config['dataset']['sample_batch_key']
        
        # Prepare columns to export
        columns_to_export = nmf_cols.copy()
        if sample_batch_key in adata.obs.columns:
            columns_to_export.append(sample_batch_key)
        
        # Add spatial coordinates if available
        if 'spatial' in adata.obsm:
            spatial_df = pd.DataFrame(
                adata.obsm['spatial'], 
                index=adata.obs_names, 
                columns=['spatial_x', 'spatial_y']
            )
            adata.obs = adata.obs.join(spatial_df)
            columns_to_export.extend(['spatial_x', 'spatial_y'])
        
        # Export to CSV
        valid_columns = [col for col in columns_to_export if col in adata.obs.columns]
        if valid_columns:
            csv_path = f"{output_path}/nmf_group_assignments.csv"
            adata.obs[valid_columns].to_csv(csv_path)
            print(f"Saved NMF group assignments to {csv_path}")
        
    except Exception as e:
        print(f"Error saving NMF results to CSV: {e}")

def run_nmf_analysis_pipeline(config, adata_vis_mapped):
    """
    Main pipeline for NMF analysis.
    
    Parameters:
    -----------
    config : dict
        Configuration dictionary
    run_name_path : str
        Path to spatial mapping results
    """
    print("\n=== Starting NMF Analysis Pipeline ===")
    

    # Setup NMF parameters
    nmf_dict_path = config['paths']['nmf_models']
    os.makedirs(nmf_dict_path, exist_ok=True)

    n_fact_range = config['nmf_analysis']['n_fact_range']
    # Handle range [start, end], [start, end, step], and single value cases
    if len(n_fact_range) == 3:
        # [start, end, step] format
        n_fact_list = list(range(n_fact_range[0], n_fact_range[1] + 1, n_fact_range[2]))
    elif len(n_fact_range) == 2 and n_fact_range[0] != n_fact_range[1]:
        # [start, end] format with step=1
        n_fact_list = list(range(n_fact_range[0], n_fact_range[1] + 1))
    else:
        # Single value case
        n_fact_list = [n_fact_range[0]]
    
    # Run NMF for each factor number
    for n_fact in n_fact_list:
        print(f"\n--- Running NMF with n_fact={n_fact} ---")

        nmf_result_dict_path = os.path.join(nmf_dict_path, f"nmf_result_dict_n_fact{n_fact}.pkl")
        nmf_adata_additions_path = os.path.join(nmf_dict_path, f"nmf_adata_additions_n_fact{n_fact}.pkl")

        # Try to load existing results
        nmf_results_dict = None
        if os.path.exists(nmf_result_dict_path) and os.path.exists(nmf_adata_additions_path):
            print("Loading existing NMF results...")
            try:
                adata_vis_mapped = restore_adata_additions(adata_vis_mapped, nmf_adata_additions_path)
                with open(nmf_result_dict_path, 'rb') as f:
                    nmf_results_dict = pickle.load(f)
                print("Successfully loaded existing NMF results.")
            except Exception as e:
                print(f"Error loading NMF files: {e}. Re-running NMF analysis.")
                nmf_results_dict = None

        # Run NMF if not loaded
        if nmf_results_dict is None:
            print("Running new NMF colocalization analysis...")
            nmf_results_dict, adata_vis_mapped_nmf = run_nmf_colocalization(
                adata_vis_mapped, config, [n_fact]
            )

            if nmf_results_dict is not None:
                # Save results
                print("Saving NMF results...")
                _ = extract_adata_additions(adata_vis_mapped, adata_vis_mapped_nmf, nmf_adata_additions_path)
                
                with open(nmf_result_dict_path, 'wb') as f:
                    pickle.dump(nmf_results_dict, f)
                
                # Update adata for analysis
                adata_vis_mapped = adata_vis_mapped_nmf
            else:
                print("NMF analysis failed. Skipping this n_fact value.")
                continue

        # no need to analyse as the `run_colocation` function already does this
        # Analyze results
        # if nmf_results_dict is not None:
        #     analyze_nmf_results(adata_vis_mapped, nmf_results_dict, config)

    print("\n=== NMF Analysis Pipeline Completed ===")

if __name__ == '__main__':
    # Load configuration
    config = load_configuration('config/cellType_config.yaml')
    run_name = config['paths']['spatial_output']
    # Load processed spatial data
    adata_vis_mapped = load_processed_spatial_data(run_name)
    # Run NMF analysis pipeline
    run_nmf_analysis_pipeline(config, adata_vis_mapped)

    import glob
    import shutil
    # rename folder as `run_colocation` mistakenly create a folder name concatenated with `CoLocatedGroupsSklearnNMF`
    nmf_analysis_path = config['paths']['nmf_analysis_path']
    # Move files from folders named CoLocatedGroupsSklearnNMF* to nmf_analysis_path, preserving structure
    for folder in glob.glob(f"{nmf_analysis_path}CoLocatedGroupsSklearnNMF*"):
        for root, dirs, files in os.walk(folder):
            rel_path = os.path.relpath(root, folder)
            dest_dir = os.path.join(nmf_analysis_path, rel_path) if rel_path != '.' else nmf_analysis_path
            os.makedirs(dest_dir, exist_ok=True)
            for file in files:
                src_file = os.path.join(root, file)
                dest_file = os.path.join(dest_dir, file)
                shutil.move(src_file, dest_file)
        shutil.rmtree(folder)
        print(f"Moved and deleted folder: {folder}")


    # Delete {nmf_analysis_path}/anndata folder if it exists
    # this is just a backup of the adata object
    anndata_folder = os.path.join(nmf_analysis_path, "anndata")
    if os.path.isdir(anndata_folder):
        shutil.rmtree(anndata_folder)
        print(f"Deleted folder and contents: {anndata_folder}")
    
    print("\nNMF analysis completed.")
