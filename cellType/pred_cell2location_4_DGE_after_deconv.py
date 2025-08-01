import scanpy as sc
import pandas as pd
import re
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import squidpy as sq
import gc
from cellType.pred_cell2location_utils import load_configuration, load_processed_spatial_data

def _cleanup_memory():
    """Force garbage collection to free memory."""
    gc.collect()
    
def _close_all_figures():
    """Close all matplotlib figures to free memory."""
    plt.close('all')
    gc.collect()

def _get_abundance_column_name(adata, cell_type_name, abundance_col_prefix="q05_"):
    """
    Helper function to determine the correct abundance column name.
    
    Parameters:
    -----------
    adata : AnnData
        Spatial AnnData object
    cell_type_name : str
        Name of the cell type
    abundance_col_prefix : str
        Prefix for abundance columns
        
    Returns:
    --------
    str or None
        Name of the abundance column if found, None otherwise
    """
    abundance_col = cell_type_name
    if abundance_col not in adata.obs.columns:
        abundance_col_alt = f"{abundance_col_prefix}{cell_type_name}"
        if abundance_col_alt in adata.obs.columns:
            abundance_col = abundance_col_alt
        else:
            return None
    return abundance_col

def _create_abundance_groups(adata, abundance_col, cell_type_name, low_quantile=0.25, high_quantile=0.75):
    """
    Helper function to create abundance groups (High/Mid/Low) for a cell type.
    
    Parameters:
    -----------
    adata : AnnData
        Spatial AnnData object
    abundance_col : str
        Name of the abundance column
    cell_type_name : str
        Name of the cell type
    low_quantile : float
        Quantile threshold for low abundance
    high_quantile : float
        Quantile threshold for high abundance
        
    Returns:
    --------
    tuple
        (group_column_name, group_counts)
    """
    low_quantile_val = adata.obs[abundance_col].quantile(low_quantile)
    high_quantile_val = adata.obs[abundance_col].quantile(high_quantile)

    group_col_name = f'{cell_type_name}_abundance_group'
    adata.obs[group_col_name] = 'Mid'
    adata.obs.loc[adata.obs[abundance_col] <= low_quantile_val, group_col_name] = 'Low'
    adata.obs.loc[adata.obs[abundance_col] >= high_quantile_val, group_col_name] = 'High'

    group_counts = adata.obs[group_col_name].value_counts()
    return group_col_name, group_counts

def _perform_dge_analysis(adata, group_column, use_raw=True):
    """
    Helper function to perform DGE analysis using scanpy.
    
    Parameters:
    -----------
    adata : AnnData
        Spatial AnnData object with group assignments
    group_column : str
        Column name containing group assignments
    cell_type_name : str
        Name of the cell type being analyzed
    config : dict
        Configuration dictionary
    use_raw : bool
        Whether to use raw counts for DGE
        
    Returns:
    --------
    bool
        True if DGE was successful, False otherwise
    """    
    # Perform DGE
    sc.tl.rank_genes_groups(adata,
                            groupby=group_column,
                            groups=['High'],
                            reference='Low',
                            method='wilcoxon',
                            use_raw=use_raw)
    
    # Check if results exist
    if ('names' in adata.uns['rank_genes_groups'] and 
        'High' in pd.DataFrame(adata.uns['rank_genes_groups']['names']).columns):
        return True
    return False

def _save_dge_results(adata, group_key, comparison_name, output_path):
    """
    Helper function to save DGE results to CSV - generalized for any comparison type.
    
    Parameters:
    -----------
    adata : AnnData
        AnnData object with DGE results
    group_key : str
        Key for the group in rank_genes_groups results
    comparison_name : str
        Descriptive name for the comparison (used in filename)
    config : dict
        Configuration dictionary
    output_path : str
        Directory to save results
    file_prefix : str
        Prefix for the output filename
    """
    
    dge_results_path = f"{output_path}/{comparison_name}.csv"
    
    dge_data = {
        'names': adata.uns['rank_genes_groups']['names'][group_key],
        'scores': adata.uns['rank_genes_groups']['scores'][group_key],
        'pvals': adata.uns['rank_genes_groups']['pvals'][group_key],
        'pvals_adj': adata.uns['rank_genes_groups']['pvals_adj'][group_key],
        'logfoldchanges': adata.uns['rank_genes_groups']['logfoldchanges'][group_key]
    }
    
    dge_results_df = pd.DataFrame(dge_data)
    dge_results_df.to_csv(dge_results_path, index=False)
    print(f"DGE results saved to: {dge_results_path}")

def _plot_spatial_comparison(adata, color_col, title, ax, vmin=None, vmax=None, cmap='magma'):
    """
    Helper function to create a single spatial scatter plot.
    
    Returns the image object for colorbar creation.
    """
    try:
        if 'spatial' in adata.obsm and len(adata) > 0:
            im = sq.pl.spatial_scatter(
                adata,
                color=color_col,
                title=title,
                cmap=cmap,
                ax=ax,
                vmin=vmin,
                vmax=vmax,
                colorbar=False,
                img_alpha =0.2,
            )
            ax.set_title(title, fontsize=10)
            return im
        else:
            ax.text(0.5, 0.5, f'No spatial data\n({len(adata)} spots)', 
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_title(title, fontsize=10)
            return None
    except Exception as e:
        ax.text(0.5, 0.5, f'Spatial plot error:\n{str(e)}', 
               ha='center', va='center', transform=ax.transAxes)
        ax.set_title(title, fontsize=10)
        return None

def _generate_dge_plots(adata, comparison_name, plot_config, output_path, file_prefix="dge"):
    """
    Generalized function to generate DGE visualization plots.
    
    Parameters:
    -----------
    adata : AnnData
        AnnData object with DGE results
    comparison_name : str
        Name for the comparison (used in filename)
    plot_config : dict
        Configuration for plotting with keys:
        - 'spatial_plots': list of dicts with keys 'data', 'color_col', 'title'
        - 'color_range': dict with 'vmin', 'vmax' (optional)
        - 'colorbar_label': str
        - 'group_key': str (key for violin plot)
    output_path : str
        Directory to save plots
    file_prefix : str
        Prefix for the output filename
    """
    pdf_plot_path = f"{output_path}/{file_prefix}_{comparison_name}_report.pdf"
    
    with PdfPages(pdf_plot_path) as pdf:
        # Page 1: Spatial plots
        spatial_plots = plot_config.get('spatial_plots', [])
        if spatial_plots:
            n_plots = len(spatial_plots)
            fig_spatial, axes = plt.subplots(1, n_plots, figsize=(7*n_plots, 6))
            if n_plots == 1:
                axes = [axes]
            
            # Get color range
            color_range = plot_config.get('color_range', {})
            vmin = color_range.get('vmin')
            vmax = color_range.get('vmax')
            
            im = None
            for i, plot_info in enumerate(spatial_plots):
                plot_data = plot_info['data']
                color_col = plot_info['color_col']
                title = plot_info['title']
                
                current_im = _plot_spatial_comparison(
                    plot_data, color_col, title, axes[i], vmin, vmax
                )
                if current_im is not None:
                    im = current_im
            
            # Add colorbar if we have a valid image
            if im is not None:
                plt.tight_layout()
                cbar_ax = fig_spatial.add_axes([0.92, 0.15, 0.02, 0.7])
                cbar = fig_spatial.colorbar(im, cax=cbar_ax)
                cbar.set_label(plot_config.get('colorbar_label', 'Abundance'), 
                              rotation=270, labelpad=15)
            
            pdf.savefig(fig_spatial, bbox_inches='tight')
            plt.close(fig_spatial)
        else:
            # Create placeholder if no spatial plots
            fig_placeholder = plt.figure(figsize=(8, 6))
            plt.text(0.5, 0.5, f'No spatial plots configured\nfor {comparison_name}', 
                    ha='center', va='center', fontsize=12)
            plt.title(f'Spatial plots - {comparison_name}', fontsize=10)
            pdf.savefig(fig_placeholder)
            plt.close(fig_placeholder)

        # Page 2: Violin plot
        group_key = plot_config.get('group_key')
        if group_key:
            try:
                sc.pl.rank_genes_groups_violin(adata, groups=[group_key], n_genes=10, show=False)
                pdf.savefig(plt.gcf(), bbox_inches='tight')
                plt.close()
            except Exception as e:
                print(f"Could not generate violin plot for {comparison_name}: {e}")
                fig_placeholder_violin = plt.figure(figsize=(10, 6))
                plt.text(0.5, 0.5, f"Violin plot could not be generated\nError: {str(e)}", 
                        ha='center', va='center', fontsize=12)
                plt.title(f'DGE Violin Plot - {comparison_name}', fontsize=10)
                pdf.savefig(fig_placeholder_violin)
                plt.close(fig_placeholder_violin)
        else:
            # Create placeholder if no group key
            fig_placeholder_violin = plt.figure(figsize=(10, 6))
            plt.text(0.5, 0.5, "No group key specified\nfor violin plot", 
                    ha='center', va='center', fontsize=12)
            plt.title(f'DGE Violin Plot - {comparison_name}', fontsize=10)
            pdf.savefig(fig_placeholder_violin)
            plt.close(fig_placeholder_violin)
    
    print(f"DGE report saved to: {pdf_plot_path}")

def dge_within_samples(adata_full, cell_type_name, config, use_raw=True):
    """
    Performs DGE between spots with high vs. low abundance of a specific cell type,
    processing each sample/slice individually.
    
    This function was previously named 'dge_by_cell_type_abundance'.
    """
    sample_batch_key = config['dataset']['sample_batch_key']
    combine_cell_types = config['shared'].get('combine_cell_types', False)
    if combine_cell_types:
        deg_output_path = config['paths']['dge_results_within_sample_combined']
    else:
        deg_output_path = config['paths']['dge_results_within_sample']
    os.makedirs(deg_output_path, exist_ok=True)

    # Determine samples to process
    if sample_batch_key not in adata_full.obs.columns:
        print(f"Sample batch key '{sample_batch_key}' not found in adata.obs. Performing DGE on the whole dataset.")
        unique_samples = ['all_samples']
    else:
        unique_samples = adata_full.obs[sample_batch_key].unique()

    # Get abundance column name
    abundance_col = _get_abundance_column_name(adata_full, cell_type_name)
    if abundance_col is None:
        print(f"Abundance column for '{cell_type_name}' not found in adata.obs. Available: {adata_full.obs.columns.tolist()}")
        return

    print(f"\n--- Performing Within-Sample DGE for cell type: {cell_type_name} ---")
    print(f"Using abundance column: {abundance_col}")

    # Process each sample
    for sample_id in unique_samples:
        print(f"\n--- Processing sample: {sample_id} ---")
        
        # Create sample subset
        if sample_id == 'all_samples':
            adata = adata_full.copy()
        else:
            adata = adata_full[adata_full.obs[sample_batch_key] == sample_id].copy()
            adata.uns['spatial'] = {
                k: v for k, v in adata.uns['spatial'].items() if k == sample_id
            }

        # Check minimum spots
        if adata.shape[0] < 10:
            print(f"Skipping sample {sample_id} due to insufficient spots ({adata.shape[0]}).")
            continue

        # Handle raw data subsetting
        effective_use_raw = _setup_raw_data_for_sample(adata, adata_full, sample_id, use_raw)

        # Create abundance groups
        group_col_name, group_counts = _create_abundance_groups(adata, abundance_col, cell_type_name)
        print(f"Group counts for {cell_type_name} in {sample_id}: \n{group_counts}")

        # ---- compute true population ratio for this cell type ----
        cell_types = adata_full.uns.get('mod', {}).get('factor_names', [])
        if cell_types is not None and len(cell_types) > 0:
            batch_total_population    = adata.obs[cell_types].sum().sum()
        else:
            batch_total_population    = adata.obs[abundance_col].sum()

        # Check if we have enough spots for DGE
        if ('High' not in group_counts or 'Low' not in group_counts or 
            group_counts['High'] < 3 or group_counts['Low'] < 3):
            print(f"Not enough spots in 'High' or 'Low' abundance groups for {sample_id}. Skipping DGE.")
            continue

        # Perform DGE analysis
        if _perform_dge_analysis(adata, group_col_name, effective_use_raw):
            # Display top results
            result_df = pd.DataFrame(adata.uns['rank_genes_groups']['names'])['High']
            print(f"\nTop DEGs for spots with High vs. Low {cell_type_name} abundance in {sample_id}:")
            print(result_df.head(10))

            # Save results
            comparison_name = f"{sample_id}_{cell_type_name.replace(' ', '_').replace('/', '_')}_high_vs_low"
            _save_dge_results(adata, 'High', comparison_name, deg_output_path)
            
            # Generate plots
            vmin = adata.obs[abundance_col].min()
            vmax = adata.obs[abundance_col].max()
            
            # Create temporary columns for plotting
            high_mask = adata.obs[f'{cell_type_name}_abundance_group'] == 'High'
            batch_cell_type_high_population = adata.obs.loc[high_mask, abundance_col].sum()
            pop_ratio_high = batch_cell_type_high_population / batch_total_population if batch_total_population > 0 else 0
            low_mask = adata.obs[f'{cell_type_name}_abundance_group'] == 'Low'
            batch_cell_type_low_population = adata.obs.loc[low_mask, abundance_col].sum()
            pop_ratio_low = batch_cell_type_low_population / batch_total_population if batch_total_population > 0 else 0

            adata.obs[f'{abundance_col}_high_only'] = adata.obs[abundance_col].copy()
            adata.obs.loc[~high_mask, f'{abundance_col}_high_only'] = float('nan')
            
            adata.obs[f'{abundance_col}_low_only'] = adata.obs[abundance_col].copy()
            adata.obs.loc[~low_mask, f'{abundance_col}_low_only'] = float('nan')
            
            plot_config = {
                'spatial_plots': [
                    {
                        'data': adata,
                        'color_col': f'{abundance_col}_high_only',
                        'title': (
                            f'{sample_id}: {cell_type_name} - High abundance '
                            f'(pop_ratio={pop_ratio_high * 100:.2f}%)'
                        )
                    },
                    {
                        'data': adata,
                        'color_col': f'{abundance_col}_low_only',
                        'title': (
                            f'{sample_id}: {cell_type_name} - Low abundance '
                            f'(pop_ratio={pop_ratio_low * 100:.2f}%)'
                        )
                    }
                ],
                'color_range': {'vmin': vmin, 'vmax': vmax},
                'colorbar_label': f'{cell_type_name} abundance',
                'group_key': 'High'
            }
            
            _generate_dge_plots(adata, comparison_name, plot_config, deg_output_path)
        else:
            print(f"DGE 'High' group results not found for {cell_type_name} in {sample_id}. Skipping saving and plotting.")
        
        # Clean up memory after each sample
        _close_all_figures()
        del adata
        _cleanup_memory()
        print(f"Memory cleaned up after processing {sample_id}")

def _setup_raw_data_for_sample(adata, adata_full, sample_id, use_raw):
    """
    Helper function to setup raw data for a sample subset.
    
    Returns:
    --------
    bool
        Whether raw data is effectively being used
    """
    if use_raw and adata_full.raw is not None:
        try:
            adata.raw = adata_full.raw[adata.obs_names, :].to_adata()
            print("Using adata.raw.X for DGE (raw data subset for current sample).")
            return True
        except Exception as e:
            print(f"Warning: Could not subset raw data for sample {sample_id}: {e}")
            print("Using adata.X for DGE instead.")
            return False
    else:
        print("Using adata.X for DGE (adata.raw is None or use_raw=False).")
        return False
def dge_across_samples(adata_full, cell_type_name, config, use_raw=True):
    """
    Performs DGE between samples for high abundance groups and low abundance groups separately.
    Compares sample1_high vs sample2_high and sample1_low vs sample2_low.
    """
    abundance_col_prefix = "q05_"
    sample_batch_key = config['dataset']['sample_batch_key']
    cross_sample_comparisons = config['shared']['cross_sample_comparisons']
    quantile_threshold = config['shared']['quantiles']['abundance_threshold']
    deg_across_sample_output_path = config['paths']['dge_results_across_sample']
    combine_cell_types = config['shared'].get('combine_cell_types', False)
    if combine_cell_types:
        deg_across_sample_output_path = config['paths']['dge_results_across_sample_combined']
    os.makedirs(deg_across_sample_output_path, exist_ok=True)

    if sample_batch_key not in adata_full.obs.columns:
        print(f"Sample batch key '{sample_batch_key}' not found in adata.obs. Cannot perform cross-sample DGE.")
        return

    print(f"\n--- Performing Cross-Sample DGE for cell type: {cell_type_name} ---")
    
    # Determine abundance column name
    abundance_col = cell_type_name
    if abundance_col not in adata_full.obs.columns:
        abundance_col_alt = f"{abundance_col_prefix}{cell_type_name}"
        if abundance_col_alt in adata_full.obs.columns:
            abundance_col = abundance_col_alt
        else:
            print(f"Abundance column for '{cell_type_name}' not found in adata.obs. Available: {adata_full.obs.columns.tolist()}")
            return
    
    print(f"Using abundance column: {abundance_col}")

    # Process each comparison pair
    for comparison in cross_sample_comparisons:
        if len(comparison) != 2:
            print(f"Skipping invalid comparison: {comparison}. Must be a pair of samples.")
            continue
            
        sample1, sample2 = comparison
        print(f"\n--- Comparing {sample1} vs {sample2} ---")
        
        # Filter data for each sample
        sample1_mask = adata_full.obs[sample_batch_key] == sample1
        sample2_mask = adata_full.obs[sample_batch_key] == sample2
        
        if not sample1_mask.any():
            print(f"No spots found for sample '{sample1}'. Skipping comparison.")
            continue
        if not sample2_mask.any():
            print(f"No spots found for sample '{sample2}'. Skipping comparison.")
            continue

        cell_types = adata_full.uns.get('mod', {}).get('factor_names', [])
        if cell_types is not None and len(cell_types) > 0:
            batch1_total_population    = adata_full.obs.loc[sample1_mask, cell_types].sum().sum()
            batch2_total_population    = adata_full.obs.loc[sample2_mask, cell_types].sum().sum()
        else:
            batch1_total_population    = adata_full.obs.loc[sample1_mask, abundance_col].sum()
            batch2_total_population    = adata_full.obs.loc[sample2_mask, abundance_col].sum()

        # Create abundance groups for each sample
        sample1_vals = adata_full.obs.loc[sample1_mask, abundance_col]
        sample2_vals = adata_full.obs.loc[sample2_mask, abundance_col]
        sample1_high_thresh = sample1_vals.quantile(quantile_threshold)
        sample1_low_thresh  = sample1_vals.quantile(1 - quantile_threshold)
        sample2_high_thresh = sample2_vals.quantile(quantile_threshold)
        sample2_low_thresh  = sample2_vals.quantile(1 - quantile_threshold)

        sample1_high_mask = sample1_mask & (adata_full.obs[abundance_col] >= sample1_high_thresh)
        sample1_low_mask  = sample1_mask & (adata_full.obs[abundance_col] <= sample1_low_thresh)
        sample2_high_mask = sample2_mask & (adata_full.obs[abundance_col] >= sample2_high_thresh)
        sample2_low_mask  = sample2_mask & (adata_full.obs[abundance_col] <= sample2_low_thresh)
        
        # Check if we have enough spots in each group
        groups_info = {
            f'{sample1}_high': sample1_high_mask.sum(),
            f'{sample1}_low': sample1_low_mask.sum(),
            f'{sample2}_high': sample2_high_mask.sum(),
            f'{sample2}_low': sample2_low_mask.sum()
        }

        population_info = {
            f'{sample1}_high': sample1_vals[sample1_vals >= sample1_high_thresh].sum(),
            f'{sample1}_low': sample1_vals[sample1_vals < sample1_low_thresh].sum(),
            f'{sample2}_high': sample2_vals[sample2_vals >= sample2_high_thresh].sum(),
            f'{sample2}_low': sample2_vals[sample2_vals < sample2_low_thresh].sum(),
            'batch1_total_population': batch1_total_population,
            'batch2_total_population': batch2_total_population
        }
        
        print(f"Group sizes: {groups_info}")
        
        # Perform DGE for high abundance groups
        if groups_info[f'{sample1}_high'] >= 3 and groups_info[f'{sample2}_high'] >= 3:
            print(f"Performing DGE: {sample1}_high vs {sample2}_high")
            _perform_cross_sample_dge_comparison(
                adata_full, sample1_high_mask, sample2_high_mask, 
                f"{sample1}_high", f"{sample2}_high", 
                cell_type_name, abundance_col, use_raw,
                deg_across_sample_output_path, "high_abundance", population_info
            )
            # Clean up memory after each comparison
            _close_all_figures()
            _cleanup_memory()
        else:
            print(f"Insufficient spots for high abundance comparison: {sample1}_high ({groups_info[f'{sample1}_high']}) vs {sample2}_high ({groups_info[f'{sample2}_high']})")
        
        # Perform DGE for low abundance groups
        if groups_info[f'{sample1}_low'] >= 3 and groups_info[f'{sample2}_low'] >= 3:
            print(f"Performing DGE: {sample1}_low vs {sample2}_low")
            _perform_cross_sample_dge_comparison(
                adata_full, sample1_low_mask, sample2_low_mask,
                f"{sample1}_low", f"{sample2}_low",
                cell_type_name, abundance_col, use_raw,
                deg_across_sample_output_path, "low_abundance", population_info
            )
            # Clean up memory after each comparison
            _close_all_figures()
            _cleanup_memory()
        else:
            print(f"Insufficient spots for low abundance comparison: {sample1}_low ({groups_info[f'{sample1}_low']}) vs {sample2}_low ({groups_info[f'{sample2}_low']})")

        # Clean up memory after each comparison pair
        print(f"Memory cleaned up after comparing {sample1} vs {sample2}")

def _perform_cross_sample_dge_comparison(adata_full, group1_mask, group2_mask, 
                                       group1_name, group2_name, cell_type_name, 
                                       abundance_col, use_raw, output_path, 
                                       abundance_level, population_info):
    """
    Helper function to perform a single cross-sample DGE comparison.
    """    
    # Create subset with only the two groups - use view instead of copy to save memory
    combined_mask = group1_mask | group2_mask
    adata_subset = adata_full[combined_mask].copy()
    
    # Handle raw data subsetting
    effective_use_raw = _setup_raw_data_for_sample(adata_subset, adata_full, "cross_sample", use_raw)
    
    # Create group labels
    adata_subset.obs['comparison_group'] = 'other'
    adata_subset.obs.loc[group1_mask[combined_mask], 'comparison_group'] = group1_name
    adata_subset.obs.loc[group2_mask[combined_mask], 'comparison_group'] = group2_name
    
    ratio1 = population_info[group1_name] / population_info['batch1_total_population'] if population_info['batch1_total_population'] > 0 else 0
    ratio2 = population_info[group2_name] / population_info['batch2_total_population'] if population_info['batch2_total_population'] > 0 else 0

    # Perform DGE
    try:
        sc.tl.rank_genes_groups(adata_subset,
                               groupby='comparison_group',
                               groups=[group1_name],
                               reference=group2_name,
                               method='wilcoxon',
                               use_raw=effective_use_raw)
        
        print(f"Top DEGs for {group1_name} vs {group2_name}:")
        if 'names' in adata_subset.uns['rank_genes_groups'] and \
           group1_name in pd.DataFrame(adata_subset.uns['rank_genes_groups']['names']).columns:
            result_df = pd.DataFrame(adata_subset.uns['rank_genes_groups']['names'])[group1_name]
            print(result_df.head(10))
            
            # Save DGE results
            safe_cell_type = cell_type_name.replace(' ', '_').replace('/', '_')
            comparison_name = f"cross_sample_{group1_name}_vs_{group2_name}_{safe_cell_type}_{abundance_level}"
            _save_dge_results(adata_subset, group1_name, comparison_name, output_path)
            
            # Generate plots with memory optimization
            vmin = adata_subset.obs[abundance_col].min()
            vmax = adata_subset.obs[abundance_col].max()
            
            # Create subsets for plotting - use views when possible
            group1_indices = adata_subset.obs['comparison_group'] == group1_name
            group2_indices = adata_subset.obs['comparison_group'] == group2_name
            
            group1_subset = adata_subset[group1_indices].copy()
            group2_subset = adata_subset[group2_indices].copy()
            
            # Update spatial info for subsets
            sample1_id = re.match(r'^(.+)_(high|low)', group1_name).group(1)
            sample2_id = re.match(r'^(.+)_(high|low)', group2_name).group(1)
            
            if 'spatial' in adata_subset.uns:
                group1_subset.uns['spatial'] = {
                    k: v for k, v in adata_subset.uns['spatial'].items() if k == sample1_id
                }
                group2_subset.uns['spatial'] = {
                    k: v for k, v in adata_subset.uns['spatial'].items() if k == sample2_id
                }
            
            plot_config = {
                'spatial_plots': [
                    {
                        'data': group1_subset,
                        'color_col': abundance_col,
                        'title': f'{group1_name}: {cell_type_name} (pop_ratio={ratio1 * 100:.2f}%)'
                    },
                    {
                        'data': group2_subset,
                        'color_col': abundance_col,
                        'title': f'{group2_name}: {cell_type_name} (pop_ratio={ratio2 * 100:.2f}%)'
                    }
                ],
                'color_range': {'vmin': vmin, 'vmax': vmax},
                'colorbar_label': f'{cell_type_name} abundance',
                'group_key': group1_name
            }
            _generate_dge_plots(adata_subset, comparison_name, plot_config, output_path)
            
            # Clean up temporary subsets
            del group1_subset, group2_subset
        
        # Clean up main subset
        del adata_subset
        _cleanup_memory()
        
    except Exception as e:
        print(f"Error performing DGE for {group1_name} vs {group2_name}: {e}")
        _cleanup_memory()

def combine_cell_type_abundances(adata, cell_type_combinations):
    """
    Combine subtypes into new cell types by summing their abundance columns.
    Adds new columns to adata.obs for each combined cell type.
    Returns a list of new combined cell type names.
    """
    combined_cell_types = []
    for combined, subtypes in cell_type_combinations.items():
        valid_subtypes = [s for s in subtypes if s in adata.obs.columns]
        if not valid_subtypes:
            continue
        adata.obs[combined] = adata.obs[valid_subtypes].sum(axis=1)
        combined_cell_types.append(combined)
    return combined_cell_types



# --- Main Execution ---
if __name__ == '__main__':
    # --- 1. Load Configuration and Data ---
    config = load_configuration('config/cellType_config.yaml')
    run_name = config['paths']['spatial_output']
    os.makedirs(run_name, exist_ok=True)
    adata_vis_mapped = load_processed_spatial_data(run_name)

    # --- Combine cell types if enabled ---
    combine_cell_types = config['shared'].get('combine_cell_types', False)
    cell_type_combinations = config['shared'].get('cell_type_combinations', {})
    all_cell_types = config['shared'].get('all_cell_types', [])
    sample_batch_key = config['dataset'].get('sample_batch_key', 'batch')
    output_path = config['paths'].get('results_folder', '.')

    if combine_cell_types and cell_type_combinations:
        print("Combining subtypes into new cell types as specified in config...")
        combined_cell_types = combine_cell_type_abundances(adata_vis_mapped, cell_type_combinations)
        
        # Get all individual cell types that are part of combinations
        cell_types_in_combinations = set()
        for combined, subtypes in cell_type_combinations.items():
            cell_types_in_combinations.update(subtypes)
        
        # Get target cell types that are NOT part of any combination
        original_target_cell_types = config['shared']['target_cell_types']
        if isinstance(original_target_cell_types, str):
            original_target_cell_types = [original_target_cell_types]
        
        remaining_individual_cell_types = [
            ct for ct in original_target_cell_types 
            if ct not in cell_types_in_combinations
        ]
        
        # Use combined cell types plus remaining individual cell types for DGE
        target_cell_types_for_dge = combined_cell_types + remaining_individual_cell_types
        print(f"Processing combined cell types: {combined_cell_types}")
        print(f"Processing remaining individual cell types: {remaining_individual_cell_types}")
    else:
        # Use original target cell types
        target_cell_types_for_dge = config['shared']['target_cell_types']
        # Handle both single string (backward compatibility) and list formats
        if isinstance(target_cell_types_for_dge, str):
            target_cell_types_for_dge = [target_cell_types_for_dge]


    use_raw_cell_type_dge = config['dge_analysis'].get('use_raw_counts_cell_type_dge', False)

    # Check if the target cell types are available
    if 'mod' in adata_vis_mapped.uns and 'factor_names' in adata_vis_mapped.uns['mod']:
        print(f"Available cell types (factor_names): {adata_vis_mapped.uns['mod']['factor_names']}")
        for target_cell_type in target_cell_types_for_dge:
            if target_cell_type not in adata_vis_mapped.obs.columns:
                print(f"Warning: '{target_cell_type}' not found in adata.obs columns.")

    # Process each target cell type
    for i, target_cell_type_for_dge in enumerate(target_cell_types_for_dge):
        print(f"\n=== Processing DGE for cell type {i+1}/{len(target_cell_types_for_dge)}: {target_cell_type_for_dge} ===")
        
        dge_within_samples(adata_vis_mapped,
                          cell_type_name=target_cell_type_for_dge,
                          config=config,
                          use_raw=use_raw_cell_type_dge)

        # Clean up memory after within-sample analysis
        _close_all_figures()
        _cleanup_memory()
        print(f"Memory cleaned up after within-sample DGE for {target_cell_type_for_dge}")

        # --- 3. Cross-Sample DGE Analysis ---
        if config['shared']['cross_sample_comparisons']:
            print(f"\n--- Starting Cross-Sample DGE Analysis for {target_cell_type_for_dge} ---")
            dge_across_samples(adata_vis_mapped,
                              cell_type_name=target_cell_type_for_dge,
                              config=config,
                              use_raw=use_raw_cell_type_dge)
            
            # Clean up memory after cross-sample analysis
            _close_all_figures()
            _cleanup_memory()
            print(f"Memory cleaned up after cross-sample DGE for {target_cell_type_for_dge}")
        else:
            print(f"No cross-sample comparisons configured. Skipping cross-sample DGE analysis for {target_cell_type_for_dge}.")

    print("\n--- DGE Analysis Script Finished ---")
    # Final cleanup
    _close_all_figures()
    _cleanup_memory()
