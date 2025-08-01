"""
===============================================================================
CLUSTER DGE INSPECTION PLOTS SCRIPT (Python Version)
===============================================================================
Purpose:
  Python version of the R script that generates inspection plots for merged 
  differential gene expression (DGE) results on clusters, cell subtypes, and 
  optionally combined cell types.

What are "Inspection Plots"?
  Diagnostic and summary visualizations designed to help users quickly assess,
  compare, and interpret the results of DGE analyses:
    - Scatterplots: Compare log2 fold changes for shared genes between comparisons
    - Stacked Bar Charts: Show number of up- and down-regulated genes per comparison
    - Volcano Plots: Display relationship between fold change and statistical significance

Usage:
  python 2_DGE_inspection_plots_from_csv.py

Dependencies:
  - pandas, numpy, matplotlib, seaborn, itertools, yaml
  - Expects merged DGE CSV files to be present in output directories

Author: Generated from R script
Date: 2024
===============================================================================
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
import yaml
import os
import warnings

class DGEInspectionPlots:
    def __init__(self):
        self.load_config()
        # Set ggplot style and better color palette
        plt.style.use('ggplot')
        self.colors = {
            'Up': '#E31A1C',      # Red
            'Down': '#1F78B4',    # Blue  
            'NS': '#CCCCCC',      # Light gray
            'scatter': '#D62728', # Dark red for scatter points
            'line': '#7F7F7F'     # Gray for reference lines
        }
        
    def load_config(self):
        """Load configuration from YAML file"""
        with open('config/batch_config.yaml', 'r') as f:
            self.config = yaml.safe_load(f)
        self.output_dirs = self.config['output_dirs']
        
    def standardize_dge_columns(self, dge):
        """Standardize column names for different DGE file formats"""
        # Check if this is the cellTypes format
        if 'names' in dge.columns and 'logfoldchanges' in dge.columns:
            # Map cellTypes format to standard format
            column_mapping = {
                'names': 'gene',
                'logfoldchanges': 'avg_log2FC',
                'pvals_adj': 'p_val_adj',
                'pvals': 'p_val'
            }
            dge = dge.rename(columns=column_mapping)
            
            # Add missing columns with default values if they don't exist
            if 'pct.1' not in dge.columns:
                dge['pct.1'] = np.nan
            if 'pct.2' not in dge.columns:
                dge['pct.2'] = np.nan
                
            # Rename cell_type to region for consistency (if region doesn't exist)
            if 'cell_type' in dge.columns and 'region' not in dge.columns:
                dge = dge.rename(columns={'cell_type': 'region'})
        
        # Ensure avg_log2FC is numeric
        dge['avg_log2FC'] = pd.to_numeric(dge['avg_log2FC'], errors='coerce')
        return dge
    
    def load_dge_data(self, dge_csv_path):
        """Load merged DGE data"""
        dge = pd.read_csv(dge_csv_path)
        # Standardize column names based on format
        dge = self.standardize_dge_columns(dge)
        return dge
    
    def plot_comparison_scatterplots(self, dge_subset, plot_data_dir, plot_pdf_dir, 
                                   lfc_threshold=0.25, pval_threshold=0.05):
        """Generate scatterplots for all pairs of comparisons (shared genes)"""
        comparisons = dge_subset['comparison'].unique()
        
        for comp1, comp2 in combinations(comparisons, 2):
            d1 = dge_subset[(dge_subset['comparison'] == comp1) & 
                           (dge_subset['p_val_adj'] < pval_threshold)]
            d2 = dge_subset[(dge_subset['comparison'] == comp2) & 
                           (dge_subset['p_val_adj'] < pval_threshold)]
            
            shared_genes = set(d1['gene']).intersection(set(d2['gene']))
            if len(shared_genes) == 0:
                continue
                
            d1_shared = d1[d1['gene'].isin(shared_genes)][['gene', 'avg_log2FC']].drop_duplicates('gene')
            d2_shared = d2[d2['gene'].isin(shared_genes)][['gene', 'avg_log2FC']].drop_duplicates('gene')
            
            plot_df = pd.merge(d1_shared, d2_shared, on='gene', suffixes=('_1', '_2'))
            plot_df['size'] = abs(plot_df['avg_log2FC_1']) + abs(plot_df['avg_log2FC_2'])
            
            # Save data
            plot_df.to_csv(os.path.join(plot_data_dir, f'scatter_{comp1}_against_{comp2}.csv'), 
                          index=False)
            
            # Create plot with ggplot styling
            fig, ax = plt.subplots(figsize=(8, 8))
            
            # Add reference lines
            ax.axhline(y=0, linestyle='--', color=self.colors['line'], alpha=0.6, linewidth=1)
            ax.axvline(x=0, linestyle='--', color=self.colors['line'], alpha=0.6, linewidth=1)
            
            # Create scatter plot with size mapping
            scatter = ax.scatter(plot_df['avg_log2FC_1'], plot_df['avg_log2FC_2'], 
                               s=plot_df['size']*15, alpha=0.7, color=self.colors['scatter'],
                               edgecolors='white', linewidth=0.5)
            
            # Add labels for most extreme points
            extreme_points = plot_df.nlargest(10, 'size')
            for _, row in extreme_points.iterrows():
                ax.annotate(row['gene'], (row['avg_log2FC_1'], row['avg_log2FC_2']), 
                           xytext=(8, 8), textcoords='offset points', fontsize=9,
                           bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8),
                           arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
            
            ax.set_xlabel(f'avg_log2FC ({comp1})', fontsize=12, fontweight='bold')
            ax.set_ylabel(f'avg_log2FC ({comp2})', fontsize=12, fontweight='bold')
            ax.set_title(f'Shared DE Genes: {comp1} vs {comp2}', fontsize=14, fontweight='bold', pad=20)
            
            # Improve grid and spines
            ax.grid(True, alpha=0.3)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            
            plt.tight_layout()
            plt.savefig(os.path.join(plot_pdf_dir, f'scatter_{comp1}_against_{comp2}.pdf'), 
                       dpi=300, bbox_inches='tight')
            plt.close()
    
    def plot_stacked_barchart(self, dge_subset, plot_data_dir, plot_pdf_dir,
                             lfc_threshold=0.25, pval_threshold=0.05):
        """Generate stacked bar chart: number of up/down-regulated genes per comparison"""
        def classify_regulation(row):
            if row['avg_log2FC'] >= lfc_threshold and row['p_val_adj'] < pval_threshold:
                return 'Up'
            elif row['avg_log2FC'] <= -lfc_threshold and row['p_val_adj'] < pval_threshold:
                return 'Down'
            else:
                return 'NS'
        
        # Make a copy to avoid SettingWithCopyWarning
        dge_subset = dge_subset.copy()
        dge_subset.loc[:, 'regulation'] = dge_subset.apply(classify_regulation, axis=1)
        
        bar_df = (dge_subset[dge_subset['regulation'] != 'NS']
                 .groupby(['comparison', 'regulation'])
                 .size()
                 .reset_index(name='n_genes'))
        
        # Save data
        bar_df.to_csv(os.path.join(plot_data_dir, 'stacked_barchart_data.csv'), index=False)
        
        # Create plot with ggplot styling
        pivot_df = bar_df.pivot(index='comparison', columns='regulation', values='n_genes').fillna(0)
        
        fig, ax = plt.subplots(figsize=(12, 7))
        
        # Create stacked bar chart with better colors
        colors_list = [self.colors['Down'], self.colors['Up']]
        pivot_df.plot(kind='bar', stacked=True, color=colors_list, ax=ax, 
                     edgecolor='white', linewidth=0.5)
        
        ax.set_title('Number of Up/Down-regulated Genes per Comparison', 
                    fontsize=14, fontweight='bold', pad=20)
        ax.set_xlabel('Comparison', fontsize=12, fontweight='bold')
        ax.set_ylabel('Number of Genes', fontsize=12, fontweight='bold')
        
        # Improve legend
        ax.legend(title='Regulation', title_fontsize=11, fontsize=10, 
                 loc='upper right', frameon=True, fancybox=True, shadow=True)
        
        # Rotate x-axis labels
        plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
        
        # Improve grid and spines
        ax.grid(True, alpha=0.3, axis='y')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        plt.tight_layout()
        plt.savefig(os.path.join(plot_pdf_dir, 'stacked_barchart.pdf'), 
                   dpi=300, bbox_inches='tight')
        plt.close()
    
    def plot_volcano_plots(self, dge_subset, plot_data_dir, plot_pdf_dir,
                          lfc_threshold=0.25, pval_threshold=0.05):
        """Generate volcano plots for each comparison"""
        def classify_regulation(row):
            if row['avg_log2FC'] >= lfc_threshold and row['p_val_adj'] < pval_threshold:
                return 'Up'
            elif row['avg_log2FC'] <= -lfc_threshold and row['p_val_adj'] < pval_threshold:
                return 'Down'
            else:
                return 'NS'
        
        comparisons = dge_subset['comparison'].unique()
        
        for comp in comparisons:
            comp_df = dge_subset[dge_subset['comparison'] == comp].copy()
            comp_df['regulation'] = comp_df.apply(classify_regulation, axis=1)
            
            # Fix: Replace zero or negative p_val_adj with small positive value to avoid log10 errors
            comp_df['p_val_adj'] = comp_df['p_val_adj'].clip(lower=1e-300)
            
            # Identify top 10 up- and down-regulated genes
            top_up = (comp_df[comp_df['regulation'] == 'Up']
                     .nlargest(10, 'avg_log2FC'))
            top_down = (comp_df[comp_df['regulation'] == 'Down']
                       .nsmallest(10, 'avg_log2FC'))
            
            # Create label column
            comp_df['label'] = ''
            comp_df.loc[comp_df['gene'].isin(top_up['gene']), 'label'] = comp_df['gene']
            comp_df.loc[comp_df['gene'].isin(top_down['gene']), 'label'] = comp_df['gene']
            
            # Save data
            comp_df.to_csv(os.path.join(plot_data_dir, f'volcano_{comp}.csv'), index=False)
            
            # Create plot with ggplot styling
            fig, ax = plt.subplots(figsize=(9, 7))
            
            # Plot points by regulation type with better styling
            for reg_type in ['NS', 'Down', 'Up']:  # Plot NS first so colored points are on top
                subset = comp_df[comp_df['regulation'] == reg_type]
                if len(subset) > 0:
                    alpha = 0.4 if reg_type == 'NS' else 0.7
                    size = 20 if reg_type != 'NS' else 15
                    ax.scatter(subset['avg_log2FC'], -np.log10(subset['p_val_adj']), 
                             color=self.colors[reg_type], alpha=alpha, s=size,
                             label=f'{reg_type} ({len(subset)})', 
                             edgecolors='white' if reg_type != 'NS' else 'none', 
                             linewidth=0.3)
            
            # Add significance and fold change threshold lines
            if pval_threshold > 0:
                ax.axhline(y=-np.log10(pval_threshold), linestyle='--', 
                          color=self.colors['line'], alpha=0.6, linewidth=1)
            ax.axvline(x=lfc_threshold, linestyle='--', 
                      color=self.colors['line'], alpha=0.6, linewidth=1)
            ax.axvline(x=-lfc_threshold, linestyle='--', 
                      color=self.colors['line'], alpha=0.6, linewidth=1)
            
            # Add labels for top genes
            labeled_genes = comp_df[comp_df['label'] != '']
            for _, row in labeled_genes.iterrows():
                ax.annotate(row['label'], (row['avg_log2FC'], -np.log10(row['p_val_adj'])), 
                           xytext=(8, 8), textcoords='offset points', fontsize=9,
                           bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8),
                           arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.1'))
            
            ax.set_xlabel('avg_log2FC', fontsize=12, fontweight='bold')
            ax.set_ylabel('-log10(p_val_adj)', fontsize=12, fontweight='bold')
            ax.set_title(f'Volcano Plot: {comp}', fontsize=14, fontweight='bold', pad=20)
            
            # Improve legend
            ax.legend(title='Regulation', title_fontsize=11, fontsize=10, 
                     loc='upper right', frameon=True, fancybox=True, shadow=True)
            
            # Improve grid and spines
            ax.grid(True, alpha=0.3)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            
            plt.tight_layout()
            plt.savefig(os.path.join(plot_pdf_dir, f'volcano_{comp}.pdf'), 
                       dpi=300, bbox_inches='tight')
            plt.close()
    
    def run_all_plots(self, dge_csv_path, plot_data_dir, plot_pdf_dir):
        """Main function to run all plots"""
        dge = self.load_dge_data(dge_csv_path)
        
        # Check for group column
        if 'group' not in dge.columns:
            raise ValueError("Column 'group' not found in DGE data. Please ensure the data contains a 'group' column.")
        
        groups = dge['group'].unique()
        print(f"Found groups: {', '.join(groups)}")
        
        for group_name in groups:
            print(f"Processing group: {group_name}")
            dge_subset = dge[dge['group'] == group_name]
            
            # Create group-specific directories
            group_plot_data_dir = os.path.join(plot_data_dir, group_name)
            group_plot_pdf_dir = os.path.join(plot_pdf_dir, group_name)
            os.makedirs(group_plot_data_dir, exist_ok=True)
            os.makedirs(group_plot_pdf_dir, exist_ok=True)
            
            # Run plots for this group subset
            self.plot_comparison_scatterplots(dge_subset, group_plot_data_dir, group_plot_pdf_dir)
            self.plot_stacked_barchart(dge_subset, group_plot_data_dir, group_plot_pdf_dir)
            self.plot_volcano_plots(dge_subset, group_plot_data_dir, group_plot_pdf_dir)

def main():
    """Main execution function"""
    plotter = DGEInspectionPlots()
    
    # Define DGE types to process
    dge_types = [
        {
            'csv_path': os.path.join(plotter.output_dirs['deg_clusters'], 'merged_dge_on_clusters.csv'),
            'plot_data_dir': os.path.join(plotter.output_dirs['deg_clusters'], 'inspection_plot_data/'),
            'plot_pdf_dir': os.path.join(plotter.output_dirs['deg_clusters'], 'inspection_plots_pdf/')
        },
        {
            'csv_path': os.path.join(plotter.output_dirs['deg_celltypes'], 'merged_dge_on_cellTypes.csv'),
            'plot_data_dir': os.path.join(plotter.output_dirs['deg_celltypes'], 'inspection_plot_data/'),
            'plot_pdf_dir': os.path.join(plotter.output_dirs['deg_celltypes'], 'inspection_plots_pdf/')
        }
    ]
    
    # Optionally add combine_celltypes if present
    if 'deg_combine_celltypes' in plotter.output_dirs:
        dge_types.append({
            'csv_path': os.path.join(plotter.output_dirs['deg_combine_celltypes'], 'merged_dge_on_cellTypes.csv'),
            'plot_data_dir': os.path.join(plotter.output_dirs['deg_combine_celltypes'], 'inspection_plot_data/'),
            'plot_pdf_dir': os.path.join(plotter.output_dirs['deg_combine_celltypes'], 'inspection_plots_pdf/')
        })
    
    # Process each DGE type
    for dge_info in dge_types:
        dge_csv_path = dge_info['csv_path']
        if not os.path.exists(dge_csv_path):
            warnings.warn(f"Merged DGE CSV file not found. Please generate it before running this script: {dge_csv_path}")
            continue
            
        print(f"Processing DGE file: {dge_csv_path}")
        
        # Create directories
        os.makedirs(dge_info['plot_data_dir'], exist_ok=True)
        os.makedirs(dge_info['plot_pdf_dir'], exist_ok=True)
        
        # Run all plots
        plotter.run_all_plots(dge_csv_path, dge_info['plot_data_dir'], dge_info['plot_pdf_dir'])

if __name__ == "__main__":
    main()
