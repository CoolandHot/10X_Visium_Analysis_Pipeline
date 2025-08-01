"""
This script generates spatial heatmaps of cell type abundance for a specified cell type
using processed spatial transcriptomics data. It supports multiple batches/samples and
outputs the plots to a single PDF file.

Usage:
    python pred_cell2location_plot_spatial_cell_abundance.py --config CONFIG_PATH --cell_type CELL_TYPE [--output OUTPUT_PDF]

Arguments:
    --config     Path to the YAML configuration file (default: config/cellType_config.yaml)
    --cell_type  Name of the cell type to plot (required; must be a column in adata.obs)
    --combined   If set, treat cell_type as a combined cell type based on config combinations
    --output     Output PDF file for the plots (default: spatial_cell_abundance.pdf)

"""

#%%
import argparse
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import squidpy as sq
from pred_cell2location_utils import load_configuration, load_processed_spatial_data
#%%
def create_combined_cell_types(adata, config):
    """Create combined cell types based on config combinations and remove original cell types."""
    if 'shared' in config and 'cell_type_combinations' in config['shared']:
        combinations = config['shared']['cell_type_combinations']
        for combined_name, component_types in combinations.items():
            # Check if all component types exist in adata.obs
            existing_components = [ct for ct in component_types if ct in adata.obs.columns]
            if existing_components:
                # Sum the abundances of existing component types
                adata.obs[combined_name] = adata.obs[existing_components].sum(axis=1)
                print(f"Created combined cell type '{combined_name}' from {len(existing_components)} components")
        # Remove all original cell type columns used in combinations
        all_components = set(ct for cts in combinations.values() for ct in cts)
        cols_to_remove = [ct for ct in all_components if ct in adata.obs.columns]
        adata.obs.drop(columns=cols_to_remove, inplace=True)

def plot_spatial_cell_abundance(adata, cell_type, sample_batch_key, output_pdf):
    batches = adata.obs[sample_batch_key].unique()
    n_batch = len(batches)
    if n_batch <= 4:
        ncols, nrows = 2, 2
        figsize = (7*ncols, 6*nrows)
    elif n_batch <= 6:
        ncols, nrows = 3, 2
        figsize = (7*ncols, 6*nrows)
    else:
        ncols = min(n_batch, 4)
        nrows = (n_batch + ncols - 1) // ncols
        figsize = (7*ncols, 6*nrows)
        
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    axes = axes.flatten() if n_batch > 1 else [axes]
    vmin = adata.obs[cell_type].min()
    vmax = adata.obs[cell_type].max()
    for i, batch in enumerate(batches):
        adata_batch = adata[adata.obs[sample_batch_key] == batch].copy()
        adata_batch.obs = adata_batch.obs.drop(columns=[sample_batch_key, 'library_id'], errors='ignore')
        if 'spatial' in adata_batch.uns and batch in adata_batch.uns['spatial']:
            adata_batch.uns['spatial'] = {k: v for k, v in adata_batch.uns['spatial'].items() if k == batch}
        if adata_batch.shape[0] == 0:
            axes[i].text(0.5, 0.5, f'No data for batch {batch}', 
                         ha='center', va='center', transform=axes[i].transAxes)
            axes[i].set_title(f"{batch}: {cell_type}")
            continue
        
        if 'spatial' in adata_batch.obsm and len(adata_batch) > 0:
            sq.pl.spatial_scatter(
                adata_batch,
                color=cell_type,
                title=f"{batch}: {cell_type}",
                cmap='Purples',
                ax=axes[i],
                vmin=vmin,
                vmax=vmax,
                img=False,
                colorbar=False,
            )
        else:
            axes[i].text(0.5, 0.5, f'No spatial data\n({len(adata_batch)} spots)', 
                         ha='center', va='center', transform=axes[i].transAxes)
            axes[i].set_title(f"{batch}: {cell_type}")
    # Hide unused axes
    for j in range(i+1, len(axes)):
        axes[j].axis('off')
    # Add colorbar
    im = axes[0].collections[0] if hasattr(axes[0], 'collections') and axes[0].collections else None
    if im is not None:
        cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
        fig.colorbar(im, cax=cbar_ax, label=f"{cell_type} abundance")
    plt.tight_layout(rect=[0, 0, 0.9, 1])
    with PdfPages(output_pdf) as pdf:
        # Save with white page background, only axes/figure area is dark
        pdf.savefig(fig, facecolor='white')  # remove bbox_inches
    plt.close(fig)
    print(f"Spatial cell abundance heatmap saved to: {output_pdf}")

#%%
parser = argparse.ArgumentParser(description="Plot spatial cell abundance heatmap for a given cell type.")
parser.add_argument('--config', type=str, default='config/cellType_config.yaml', help='Path to config/cellType_config.yaml')
parser.add_argument('--cell_type', type=str, required=True, help='Cell type to plot')
parser.add_argument('--output', type=str, default='spatial_cell_abundance.pdf', help='Output PDF file')
parser.add_argument('--combined', action='store_true', help='If set, treat cell_type as a combined cell type')
args = parser.parse_args()

config = load_configuration(args.config)
run_name = config['paths']['spatial_output']
adata = load_processed_spatial_data(run_name)
sample_batch_key = config['dataset']['sample_batch_key']

if args.combined:
    # Only create combined cell types if requested
    create_combined_cell_types(adata, config)
    combined_types = []
    if 'shared' in config and 'cell_type_combinations' in config['shared']:
        combined_types = list(config['shared']['cell_type_combinations'].keys())
    if args.cell_type not in combined_types or args.cell_type not in adata.obs.columns:
        print(f"Combined cell type '{args.cell_type}' not found.")
        print("Available combined cell types:")
        print(f"  {combined_types}")
        exit(1)
else:
    # Only allow original cell types
    individual_types = [col for col in adata.obs.columns if 'shared' in config and 'all_cell_types' in config['shared'] and col in config['shared']['all_cell_types']]
    if args.cell_type not in individual_types:
        print(f"Cell type '{args.cell_type}' not found in adata.obs columns.")
        print("Available individual cell types:")
        print(f"  {individual_types}")
        exit(1)

plot_spatial_cell_abundance(adata, args.cell_type, sample_batch_key, args.output)
