import scanpy as sc
import matplotlib.pyplot as plt
import os

import cell2location
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel
from cellType.pred_cell2location_utils import load_configuration

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42  # enables correct plotting of text for PDFs

# Load configuration
config = load_configuration('config/cellType_config.yaml')

# Set up results folder
results_folder = config['paths']['results_folder']
ref_run_name = config['paths']['reference_output']

# Create directories if they don't exist
os.makedirs(results_folder, exist_ok=True)
os.makedirs(ref_run_name, exist_ok=True)

print("Loading GBM reference data...")
# Load your GBM reference data
adata_ref = sc.read_h5ad(config['paths']['reference_data'])
print(f"GBM reference data shape: {adata_ref.shape}")
print(f"Available samples: {adata_ref.obs[config['dataset']['reference_batch_key']].unique()}")
if len(config['dataset']['reference_covariate_keys']) > 0:
    print(f"Available covariates: {adata_ref.obs[config['dataset']['reference_covariate_keys'][0]].unique()}")

print(f"Available cell types: {adata_ref.obs[config['dataset']['reference_labels_key']].unique()}")

# Gene filtering for reference data
print("\nFiltering genes in reference data...")
selected = filter_genes(
    adata_ref, 
    cell_count_cutoff=config['gene_filtering']['cell_count_cutoff'],
    cell_percentage_cutoff2=config['gene_filtering']['cell_percentage_cutoff2'],
    nonz_mean_cutoff=config['gene_filtering']['nonz_mean_cutoff']
)
adata_ref = adata_ref[:, selected].copy()
print(f"Reference data after filtering: {adata_ref.shape}")

# Setup anndata for reference regression model
print("\nSetting up reference regression model...")
cell2location.models.RegressionModel.setup_anndata(
    adata=adata_ref,
    batch_key=config['dataset']['reference_batch_key'],
    labels_key=config['dataset']['reference_labels_key'],
    categorical_covariate_keys=config['dataset']['reference_covariate_keys']
)

# Create and train the regression model
print("Creating regression model...")
mod_ref = RegressionModel(adata_ref)
mod_ref.view_anndata_setup()

print("Training regression model...")
mod_ref.train(
    max_epochs=config['reference_model']['max_epochs'],
    accelerator=config['reference_model']['accelerator'],
    device=config['reference_model']['device']
)

# Export posterior and save reference model
print("Exporting reference model results...")
adata_ref = mod_ref.export_posterior(
    adata_ref, 
    sample_kwargs={
        'num_samples': config['reference_model']['num_samples'],
        'batch_size': config['reference_model']['batch_size'],
        'accelerator': config['reference_model']['accelerator'],
        'device': config['reference_model']['device']
    }
)

# Save reference model
mod_ref.save(f"{ref_run_name}", overwrite=True)
adata_ref.write(f"{ref_run_name}/sc_ref.h5ad")

# Plot training history
print("Plotting training history...")
mod_ref.plot_history(config['quality_control']['plot_history_skip_epochs'])
plt.savefig(f'{ref_run_name}/training_history.png', bbox_inches='tight')

# Plot QC
print("Plotting reference model QC...")
mod_ref.plot_QC()
plt.savefig(f'{ref_run_name}/qc_plots.png', bbox_inches='tight')


# Extract reference signatures
print("Extracting reference signatures...")
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']
print(f"Reference signatures shape: {inf_aver.shape}")
print("Reference signatures preview:")
print(inf_aver.iloc[0:5, 0:5])

# Save reference signatures separately for easy loading
inf_aver.to_csv(f"{ref_run_name}/reference_signatures.csv")

print("Reference model training completed successfully!")
print(f"Reference model saved in: {ref_run_name}")
print("Key outputs:")
print(f"- Trained model: {ref_run_name}/")
print(f"- Reference data: {ref_run_name}/sc_ref.h5ad")
print(f"- Reference signatures: {ref_run_name}/reference_signatures.csv")
