#!/usr/bin/env python3
"""
DIFFERENTIAL GENE EXPRESSION ANALYSIS AND VISUALIZATION BY CLUSTERS
Optimized Python version for large datasets using .h5ad format
"""

import scanpy as sc
import pandas as pd
import numpy as np
import yaml
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.sparse import csr_matrix
import squidpy as sq
import warnings
warnings.filterwarnings('ignore')

# Set scanpy settings
sc.settings.verbosity = 1  # verbosity level
sc.settings.set_figure_params(dpi=80, facecolor='white')

class DGEAnalyzer:
    def __init__(self, config_path='config/batch_config.yaml'):
        """Initialize with configuration"""
        self.config = self._load_config(config_path)
        self.project_dir = self.config['project_dir']
        self.output_file_prefix = self.config['output_file_prefix']
        self.cluster_method = self.config['cluster_method']
        self.visium_hd = self.config['VisiumHD']
        
        # Extract cluster assignments
        dge_config = self.config['Differential_Gene_Analysis']
        self.out_tumour_clusters = dge_config['outTumour_cluster_nums']
        self.in_tumour_clusters = dge_config['inTumour_cluster_nums']
        self.edge_tumour_clusters = dge_config['edgeTumour_cluster_nums']
        self.across_batch_comparisons = dge_config['across_batch_comparisons']
        self.batch_names = self.config['batch_names']
        
        # Output directories
        self.output_dirs = self.config['output_dirs']
        self._create_directories()
        
    def _load_config(self, config_path):
        """Load YAML configuration"""
        with open(config_path, 'r') as f:
            return yaml.safe_load(f)
    
    def _create_directories(self):
        """Create output directories"""
        for dir_path in [
            f"{self.output_dirs['deg_clusters']}dge_results_across_sample/",
            f"{self.output_dirs['deg_clusters']}dge_results_within_sample/"
        ]:
            Path(dir_path).mkdir(parents=True, exist_ok=True)
    
    def load_data(self):
        """Load and prepare the AnnData object"""
        print("Loading data...")
        h5ad_path = f"{self.project_dir}rds_data/{self.output_file_prefix}_clustered_12k.h5ad"
        self.adata = sc.read_h5ad(h5ad_path)
        
        # Validate cluster method
        if self.cluster_method not in self.adata.obs.columns:
            # Try loading from CSV
            cluster_csv_path = f"{self.output_dirs['clustering']}{self.cluster_method}.csv"
            if Path(cluster_csv_path).exists():
                cluster_data = pd.read_csv(cluster_csv_path, index_col=0)
                self.adata.obs[self.cluster_method] = cluster_data.iloc[:, 0].astype(str)
            else:
                raise ValueError(f"Cluster method {self.cluster_method} not found")
        
        # Set cluster identities
        self.adata.obs['cluster'] = self.adata.obs[self.cluster_method].astype(str)
        print(f"Loaded data with {self.adata.n_obs} cells and {self.adata.n_vars} genes")
    
    def assign_regions(self):
        """Assign region labels based on cluster numbers"""
        print("Assigning regions...")
        
        def assign_region(cluster):
            # Handle NaN or missing values
            if pd.isna(cluster) or cluster == 'nan' or str(cluster).lower() == 'nan':
                return "Unknown"
            
            try:
                cluster_num = int(float(str(cluster)))  # Convert to float first, then int to handle string numbers
                if cluster_num in self.out_tumour_clusters:
                    return "outTumour"
                elif cluster_num in self.in_tumour_clusters:
                    return "inTumour"
                elif cluster_num in self.edge_tumour_clusters:
                    return "edgeTumour"
                else:
                    return str(cluster_num)  # Return as string for consistency
            except (ValueError, TypeError):
                # If conversion fails, return the original value as string
                return str(cluster)
        
        self.adata.obs['cell_region_cluster'] = self.adata.obs['cluster'].apply(assign_region)
        self.adata.obs['batch_cluster'] = (
            self.adata.obs['batch'].astype(str) + "_" + 
            self.adata.obs['cell_region_cluster'].astype(str)
        )
        
        print("Region assignment completed")
        print(self.adata.obs['cell_region_cluster'].value_counts())
    
    def run_differential_expression(self, group1, group2, group_by_col, min_cells=3):
        """
        Run differential expression analysis between two groups
        Key difference: Uses scanpy's rank_genes_groups with group1 as test group vs group2 as reference
        This means positive avg_log2FC indicates higher expression in group1
        """
        # Filter for cells in the two groups
        mask = self.adata.obs[group_by_col].isin([group1, group2])
        if mask.sum() == 0:
            return None
            
        adata_subset = self.adata[mask].copy()
        
        # Check minimum cell count
        group_counts = adata_subset.obs[group_by_col].value_counts()
        if any(group_counts < min_cells):
            print(f"Skipping comparison {group1} vs {group2}: insufficient cells")
            return None
        
        # Ensure both groups exist in the subset
        unique_groups = adata_subset.obs[group_by_col].unique()
        if group1 not in unique_groups or group2 not in unique_groups:
            print(f"Skipping comparison {group1} vs {group2}: groups not found in subset")
            return None
        
        # Convert to categorical with only the present categories to avoid issues
        adata_subset.obs[group_by_col] = adata_subset.obs[group_by_col].astype('category')
        adata_subset.obs[group_by_col] = adata_subset.obs[group_by_col].cat.remove_unused_categories()
        
        # Run differential expression
        sc.tl.rank_genes_groups(
            adata_subset,
            groupby=group_by_col,
            groups=[group1],  # Test group
            reference=group2,  # Reference group
            method='wilcoxon',
            pts=True,
            use_raw=False,
            n_genes=adata_subset.n_vars
        )
        
        result = sc.get.rank_genes_groups_df(adata_subset, group=group1)
        
        # Get cell masks for both groups
        group1_cells = (adata_subset.obs[group_by_col] == group1).values
        group2_cells = (adata_subset.obs[group_by_col] == group2).values
        
        # Calculate log fold changes and percentages
        result_enhanced = []
        for _, row in result.iterrows():
            gene = row['names']

            if gene not in adata_subset.var_names:
                continue
                
            result_enhanced.append({
                'gene': gene,
                'avg_log2FC': row['logfoldchanges'],
                'pct.1': row['pct_nz_group'] if 'pct_nz_group' in row else None,
                'p_val_adj': row['pvals_adj'],
                'p_val': row['pvals'],
                'zscore': row['scores'] if 'scores' in row else None
            })
        
        result_df = pd.DataFrame(result_enhanced)
        result_df = result_df.sort_values('avg_log2FC', ascending=False)
        
        return result_df
    
    def export_scaled_data(self, cells, region, batch_id):
        """Export scaled data for specific cells"""
        output_path = f"{self.output_dirs['deg_clusters']}gene_counts_{region}_{batch_id}_scaled.csv"
        
        if Path(output_path).exists():
            return
        
        # Get scaled data for specific cells
        cell_indices = [i for i, cell in enumerate(self.adata.obs_names) if cell in cells]
        if not cell_indices:
            return
            
        # Use scaled data if available, otherwise scale on the fly
        if 'X_scaled' in self.adata.layers:
            scaled_data = self.adata.layers['X_scaled'][cell_indices, :].toarray() if hasattr(self.adata.layers['X_scaled'], 'toarray') else self.adata.layers['X_scaled'][cell_indices, :]
        else:
            # Scale data on the fly
            scaled_data = stats.zscore(self.adata.X[cell_indices, :].toarray() if hasattr(self.adata.X, 'toarray') else self.adata.X[cell_indices, :], axis=0)
        
        # Create DataFrame
        scaled_df = pd.DataFrame(
            scaled_data,
            index=[self.adata.obs_names[i] for i in cell_indices],
            columns=self.adata.var_names
        )
        scaled_df.index.name = 'barcode'
        
        # Save to CSV
        scaled_df.to_csv(output_path)
        print(f"Exported scaled data: {output_path}")
    
    def run_across_sample_analysis(self):
        """Run across-sample differential expression analysis"""
        print("Running across-sample analysis...")
        
        for compare_pair in self.across_batch_comparisons:
            batch_id1, batch_id2 = compare_pair
            print(f"Comparing {batch_id1} vs {batch_id2}")
            
            for region in ["inTumour", "outTumour", "edgeTumour"]:
                group1 = f"{batch_id1}_{region}"
                group2 = f"{batch_id2}_{region}"
                
                # Export scaled data
                group1_cells = self.adata.obs[self.adata.obs['batch_cluster'] == group1].index
                group2_cells = self.adata.obs[self.adata.obs['batch_cluster'] == group2].index
                
                if len(group1_cells) > 0:
                    self.export_scaled_data(group1_cells, region, batch_id1)
                if len(group2_cells) > 0:
                    self.export_scaled_data(group2_cells, region, batch_id2)
                
                # Run differential expression
                result = self.run_differential_expression(
                    group1, group2, 'batch_cluster'
                )
                
                if result is not None:
                    output_path = f"{self.output_dirs['deg_clusters']}dge_results_across_sample/{region}_{batch_id1}_vs_{batch_id2}.csv"
                    result.to_csv(output_path, index=False)
                    print(f"Saved: {output_path}")
    
    def run_within_sample_analysis(self):
        """Run within-sample differential expression analysis"""
        print("Running within-sample analysis...")
        
        for batch_group in self.batch_names:
            print(f"Processing batch: {batch_group}")
            
            for region1 in ["inTumour", "edgeTumour"]:
                for region2 in ["outTumour"]:
                    group1 = f"{batch_group}_{region1}"
                    group2 = f"{batch_group}_{region2}"
                    
                    # Export scaled data for group1
                    group1_cells = self.adata.obs[self.adata.obs['batch_cluster'] == group1].index
                    if len(group1_cells) > 0:
                        self.export_scaled_data(group1_cells, region1, batch_group)
                    
                    # Run differential expression
                    result = self.run_differential_expression(
                        group1, group2, 'batch_cluster'
                    )
                    
                    if result is not None:
                        output_path = f"{self.output_dirs['deg_clusters']}dge_results_within_sample/{region1}_against_{region2}_{batch_group}.csv"
                        result.to_csv(output_path, index=False)
                        print(f"Saved: {output_path}")
    
    def merge_results(self):
        """Merge all differential expression results"""
        print("Merging results...")
        
        # Find all result files
        result_files = []
        for pattern in ["dge_results_across_sample/*.csv", "dge_results_within_sample/*.csv"]:
            result_files.extend(Path(self.output_dirs['deg_clusters']).glob(pattern))
        
        merged_data = []
        
        for file_path in result_files:
            if file_path.stat().st_size == 0:  # Skip empty files
                continue
                
            try:
                data = pd.read_csv(file_path)
                if data.empty:
                    continue
                    
                data['File_Name'] = file_path.name
                
                # Parse file name for metadata
                if "_against_" in file_path.name:
                    data['group'] = 'within_sample'
                    # Extract region from filename like "inTumour_against_outTumour_batch.csv"
                    parts = file_path.stem.split('_')
                    if len(parts) >= 4:
                        data['region'] = f"{parts[0]}_against_{parts[2]}"
                        data['comparison'] = '_'.join(parts[3:])
                    else:
                        data['region'] = file_path.stem
                        data['comparison'] = file_path.stem
                else:
                    data['group'] = 'across_sample'
                    # Extract region from filename like "region_batch1_vs_batch2.csv"
                    parts = file_path.stem.split('_')
                    if len(parts) >= 1:
                        data['region'] = parts[0]
                        data['comparison'] = '_'.join(parts[1:])
                    else:
                        data['region'] = file_path.stem
                        data['comparison'] = file_path.stem
                
                merged_data.append(data)
                
            except Exception as e:
                print(f"Error processing {file_path}: {e}")
                continue
        
        if merged_data:
            # Combine all data
            final_data = pd.concat(merged_data, ignore_index=True)
            
            # Reorder columns
            column_order = ['gene', 'avg_log2FC', 'group', 'region', 'comparison']
            remaining_cols = [col for col in final_data.columns if col not in column_order and col != 'File_Name']
            final_columns = column_order + remaining_cols
            
            final_data = final_data[final_columns]
            final_data = final_data.sort_values(['gene', 'avg_log2FC'], ascending=[True, False])
            
            # Save merged results
            output_path = f"{self.output_dirs['deg_clusters']}merged_dge_on_clusters.csv"
            final_data.to_csv(output_path, index=False)
            print(f"Merged results saved: {output_path}")
        else:
            print("No valid result files found for merging")
    
    def create_spatial_visualizations(self):
        """Create spatial visualizations for regions using squidpy"""
        print("Creating spatial visualizations...")
        
        # Define colors for the main regions only
        region_colors = {
            'outTumour': '#1F78B4',    # Blue
            'inTumour': '#E31A1C',     # Red
            'edgeTumour': '#33A02C'    # Green
        }
        
        try:
            if 'spatial' in self.adata.obsm.keys():
                # Check if we have library_id information
                has_library_id = 'library_id' in self.adata.uns.get('spatial', {})
                
                for compare_pair in self.across_batch_comparisons:
                    batch_id1, batch_id2 = compare_pair
                    output_dir = f"{self.output_dirs['deg_clusters']}"
                    
                    for batch_id in [batch_id1, batch_id2]:
                        batch_mask = self.adata.obs['batch'] == batch_id
                        if not batch_mask.any():
                            continue
                            
                        adata_batch = self.adata[batch_mask].copy()
                        
                        # Filter to only include main tumor regions for plotting
                        main_regions = ['outTumour', 'inTumour', 'edgeTumour']
                        main_region_mask = adata_batch.obs['cell_region_cluster'].isin(main_regions)
                        adata_main_regions = adata_batch[main_region_mask].copy()
                        
                        # Get available regions in this batch
                        batch_regions = [r for r in main_regions if r in adata_main_regions.obs['cell_region_cluster'].unique()]
                        
                        if not batch_regions:
                            print(f"No main tumor regions found in batch {batch_id}, skipping...")
                            continue
                        
                        # Calculate number of subplots needed (All Regions + individual regions)
                        n_plots = len(batch_regions) + 1  # +1 for combined plot
                        n_cols = min(4, n_plots)
                        n_rows = (n_plots + n_cols - 1) // n_cols
                        
                        fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 5*n_rows))
                        if n_plots == 1:
                            axes = [axes]
                        elif n_rows == 1:
                            axes = axes.flatten()
                        else:
                            axes = axes.flatten()
                        
                        # Plot all main regions together
                        try:
                            sq.pl.spatial_scatter(
                                adata_main_regions,
                                color='cell_region_cluster',
                                palette=region_colors,
                                size=3,
                                ax=axes[0],
                                title=f'{batch_id} - All Regions',
                                library_key=batch_id if has_library_id else None
                            )
                        except ValueError as e:
                            if "library_key" in str(e):
                                # Fall back to matplotlib plotting if squidpy fails
                                coords = adata_main_regions.obsm['spatial']
                                colors = [region_colors.get(r, 'lightgray') for r in adata_main_regions.obs['cell_region_cluster']]
                                axes[0].scatter(coords[:, 0], coords[:, 1], c=colors, s=3)
                                axes[0].set_aspect('equal')
                                axes[0].axis('off')
                                axes[0].set_title(f'{batch_id} - All Regions')
                            else:
                                raise e
                        
                        # Plot individual main regions
                        for idx, region in enumerate(batch_regions):
                            if idx + 1 >= len(axes):
                                break
                                
                            region_cells = adata_batch.obs['cell_region_cluster'] == region
                            if not region_cells.any():
                                continue
                            
                            # Create binary mask for highlighting on the full batch data
                            adata_batch.obs[f'{region}_highlight'] = region_cells.astype(str)
                            
                            try:
                                sq.pl.spatial_scatter(
                                    adata_batch,
                                    color=f'{region}_highlight',
                                    palette={'True': region_colors[region], 'False': 'lightgray'},
                                    size=3,
                                    ax=axes[idx + 1],
                                    title=f'{batch_id} - {region}',
                                    legend_loc='none',
                                    library_key=batch_id if has_library_id else None
                                )
                            except ValueError as e:
                                if "library_key" in str(e):
                                    # Fall back to matplotlib plotting
                                    coords = adata_batch.obsm['spatial']
                                    colors = ['lightgray'] * len(coords)
                                    region_indices = region_cells.values
                                    for i, is_region in enumerate(region_indices):
                                        if is_region:
                                            colors[i] = region_colors[region]
                                    
                                    axes[idx + 1].scatter(coords[:, 0], coords[:, 1], c=colors, s=3)
                                    axes[idx + 1].set_aspect('equal')
                                    axes[idx + 1].axis('off')
                                    axes[idx + 1].set_title(f'{batch_id} - {region}')
                                else:
                                    raise e
                        
                        # Hide unused subplots
                        for idx in range(n_plots, len(axes)):
                            axes[idx].set_visible(False)
                        
                        plt.tight_layout()
                        output_path = f"{output_dir}{batch_id}_spatial_regions.pdf"
                        plt.savefig(output_path, bbox_inches='tight', dpi=300)
                        plt.close()
                        
                        print(f"Saved spatial plot: {output_path}")
            else:
                print("No spatial coordinates found in data - skipping spatial plots")
                
        except Exception as e:
            print(f"Error creating spatial visualizations: {e}")
            import traceback
            traceback.print_exc()

    def run_analysis(self):
        """Run the complete analysis pipeline"""
        print("Starting DGE analysis pipeline...")
        
        # Load data
        self.load_data()
        
        # Assign regions
        self.assign_regions()
        
        # Run differential expression analyses
        self.run_across_sample_analysis()
        self.run_within_sample_analysis()
        
        # Merge results
        self.merge_results()
        
        # Create visualizations
        self.create_spatial_visualizations()
        
        print("Analysis completed successfully!")

def main():
    """Main execution function"""
    analyzer = DGEAnalyzer()
    analyzer.run_analysis()

if __name__ == "__main__":
    main()