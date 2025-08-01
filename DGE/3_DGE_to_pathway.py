#!/usr/bin/env python3
"""
Refactored pathway analysis script with class-based architecture
- Reads DGE CSV and runs KEGG enrichment analysis
- Creates visualizations and analyzes shared pathways
- Outputs organized data to CSV and logs to files
"""

import itertools
from matplotlib.backends.backend_pdf import PdfPages
import os
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import seaborn as sns
import logging
from typing import Dict, List, Set, Any
import yaml
import glob

class OutputManager:
    """Manages all file outputs including CSV, logs, and plots"""
    
    def __init__(self, deg_clusters_dir, pathway_results_dir, pathway_plots_dir):
        self.dirs = {
            'deg_clusters': deg_clusters_dir,
            'pathway_results': pathway_results_dir,
            'pathway_plots': pathway_plots_dir,
            'pathway_logs': pathway_results_dir
        }
        # Create directories
        for dir_path in self.dirs.values():
            os.makedirs(dir_path, exist_ok=True)
        self._setup_logging()
    
    def _setup_logging(self):
        """Setup logging configuration"""
        log_file = os.path.join(self.dirs['pathway_logs'], 'pathway_analysis.log')
        logging.basicConfig(
            filename=log_file,
            level=logging.INFO,
            format='%(message)s',
            filemode='w'
        )
        self.logger = logging.getLogger(__name__)
    
    def save_csv(self, data: pd.DataFrame, filename: str, subdir: str = 'pathway_results') -> str:
        """Save structured data to CSV"""
        filepath = os.path.join(self.dirs[subdir], f"{filename}.csv")
        data.to_csv(filepath, index=False)
        self.logger.info(f"CSV saved: {filepath}")
        return filepath
    
    def save_plot(self, figure, filename: str, subdir: str = 'pathway_plots') -> str:
        """Save plot to PDF"""
        filepath = os.path.join(self.dirs[subdir], f"{filename}.pdf")
        figure.savefig(filepath, format='pdf', bbox_inches='tight', dpi=300)
        plt.close(figure)
        self.logger.info(f"Plot saved: {filepath}")
        return filepath
    
    def log_info(self, message: str):
        """Log informational message"""
        self.logger.info(message)
    
    def log_analysis_summary(self, title: str, data: Dict[str, Any]):
        """Log analysis summary data"""
        self.logger.info(f"\n{'='*60}")
        self.logger.info(f"{title}")
        self.logger.info(f"{'='*60}")
        for key, value in data.items():
            self.logger.info(f"{key}: {value}")


class EnrichrAnalyzer:
    """Handles KEGG enrichment analysis using Enrichr"""
    
    def __init__(self, output_manager: OutputManager):
        self.output_manager = output_manager
        self.results = []
    
    def run_enrichment_analysis(self, df: pd.DataFrame) -> List[pd.DataFrame]:
        """Run Enrichr analysis for each group-region-comparison combination"""
        self.output_manager.log_info("Starting Enrichr analysis...")
        
        for (group, region, comparison), subdf in df.groupby(['group', 'region', 'comparison']):
            self.output_manager.log_info(f"Processing: {group} - {region} - {comparison}")
            
            # Filter by adjusted p-value
            subdf = subdf[subdf['p_val_adj'] < 0.05]
            
            if subdf.empty:
                self.output_manager.log_info(f"No significant genes for {group}-{region}-{comparison}")
                continue
            
            # Prepare gene list
            ranked = subdf.sort_values('avg_log2FC', ascending=False)
            gene_list = ranked['gene'].dropna().tolist()
            
            if not gene_list:
                self.output_manager.log_info(f"No valid genes for {group}-{region}-{comparison}")
                continue
            
            try:
                # Run Enrichr
                enr = gp.enrichr(
                    gene_list=gene_list,
                    gene_sets='KEGG_2019_Mouse',
                    organism='Mouse',
                    outdir=None,
                    cutoff=0.05
                )
                
                enrichr_df = enr.results
                enrichr_df = enrichr_df[enrichr_df['Adjusted P-value'] < 0.05]
                
                if not enrichr_df.empty:
                    # Create bar plot
                    self._create_enrichment_plot(enrichr_df, group, region, comparison)
                    
                    # Add metadata
                    enrichr_df = enrichr_df.copy()
                    enrichr_df['group'] = group
                    enrichr_df['region'] = region
                    enrichr_df['comparison'] = comparison
                    self.results.append(enrichr_df)
                    
                    self.output_manager.log_info(f"Found {len(enrichr_df)} significant pathways")
                else:
                    self.output_manager.log_info(f"No significant pathways for {group}-{region}-{comparison}")
                    
            except Exception as e:
                self.output_manager.log_info(f"Error in enrichment for {group}-{region}-{comparison}: {str(e)}")
        
        return self.results
    
    def _create_enrichment_plot(self, enrichr_df: pd.DataFrame, group: str, region: str, comparison: str):
        """Create and save enrichment bar plot"""
        fig, ax = plt.subplots(figsize=(10, 7))
        
        plot_df = enrichr_df.sort_values('Combined Score')
        sns.barplot(y='Term', x='Combined Score', data=plot_df, ax=ax)
        
        ax.set_title(f"KEGG Enrichment: {comparison} in {region} ({group})")
        ax.set_xlabel('Combined Score')
        ax.set_ylabel('Pathway')
        
        filename = f"enrichment_barplot_{group}_{region}_{comparison}"
        self.output_manager.save_plot(fig, filename)
    
    def get_merged_results(self) -> pd.DataFrame:
        """Get all results merged into single DataFrame"""
        if not self.results:
            raise ValueError("No enrichment results found")
        
        merged = pd.concat(self.results, ignore_index=True)
        merged = merged.sort_values(['Term', 'Combined Score'], ascending=[True, False])
        return merged


class VennAnalyzer:
    """Handles Venn diagram analysis for pathway overlaps"""
    
    def __init__(self, output_manager: OutputManager):
        self.output_manager = output_manager
    
    def create_sample_comparison_diagrams(self, df: pd.DataFrame):
        """Create Venn diagrams showing shared terms across regions for each group.
        For >3 regions, generate all 2- and 3-region combinations as separate Venn diagrams,
        arranged in a multi-page PDF (3x2 or 4x2 per page).
        """
        import math
        groups = df['group'].unique()
        group_summaries = []
        all_venn_plots = []
        for group in groups:
            group_data = df[df['group'] == group]
            regions = group_data['region'].unique()
            region_terms = {region: set(group_data[group_data['region'] == region]['Term']) for region in regions}
            # Collect summary data
            for region in regions:
                group_summaries.append({
                    'group': group,
                    'region': region,
                    'unique_terms': len(region_terms[region])
                })
            venn_jobs = []
            if len(regions) == 2:
                region_list = list(regions)
                def plot_venn2(ax, r1=region_list[0], r2=region_list[1], terms=region_terms, g=group):
                    venn2([terms[r1], terms[r2]], set_labels=[r1, r2], ax=ax)
                    ax.set_title(f'{g}: {r1} vs {r2}', fontsize=14)
                venn_jobs.append((f'{group}: {region_list[0]} vs {region_list[1]}', plot_venn2))
            elif len(regions) == 3:
                region_list = list(regions)
                def plot_venn3(ax, r1=region_list[0], r2=region_list[1], r3=region_list[2], terms=region_terms, g=group):
                    venn3([terms[r1], terms[r2], terms[r3]], set_labels=[r1, r2, r3], ax=ax)
                    ax.set_title(f'{g}: {r1}, {r2}, {r3}', fontsize=14)
                venn_jobs.append((f'{group}: {region_list[0]}, {region_list[1]}, {region_list[2]}', plot_venn3))
            elif len(regions) > 3:
                # All 2-region combinations
                region_pairs = list(itertools.combinations(regions, 2))
                for r1, r2 in region_pairs:
                    def plot_venn2(ax, r1=r1, r2=r2, terms=region_terms, g=group):
                        venn2([terms[r1], terms[r2]], set_labels=[r1, r2], ax=ax)
                        ax.set_title(f'{g}: {r1} vs {r2}', fontsize=12)
                    venn_jobs.append((f'{group}: {r1} vs {r2}', plot_venn2))
                # All 3-region combinations
                region_triplets = list(itertools.combinations(regions, 3))
                for r1, r2, r3 in region_triplets:
                    def plot_venn3(ax, r1=r1, r2=r2, r3=r3, terms=region_terms, g=group):
                        venn3([terms[r1], terms[r2], terms[r3]], set_labels=[r1, r2, r3], ax=ax)
                        ax.set_title(f'{g}: {r1}, {r2}, {r3}', fontsize=12)
                    venn_jobs.append((f'{group}: {r1}, {r2}, {r3}', plot_venn3))
            all_venn_plots.extend(venn_jobs)
        # Layout: 3x2 if <=6 per page, else 4x2 (8 per page)
        per_page = 6 if len(all_venn_plots) <= 24 else 8
        n_pages = math.ceil(len(all_venn_plots) / per_page)
        from matplotlib.backends.backend_pdf import PdfPages
        pdf_path = os.path.join(self.output_manager.dirs['pathway_plots'], "sample_comparison_diagrams.pdf")
        with PdfPages(pdf_path) as pdf:
            for page in range(n_pages):
                n_this_page = min(per_page, len(all_venn_plots) - page * per_page)
                nrows = 3 if per_page == 6 else 4
                ncols = 2
                fig, axes = plt.subplots(nrows, ncols, figsize=(16, nrows * 5))
                axes = axes.flatten()
                for i in range(per_page):
                    idx = page * per_page + i
                    if i >= n_this_page or idx >= len(all_venn_plots):
                        axes[i].set_visible(False)
                        continue
                    title, plot_func = all_venn_plots[idx]
                    plot_func(axes[i])
                plt.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)
        # Save summary data
        summary_df = pd.DataFrame(group_summaries)
        self.output_manager.save_csv(summary_df, "region_term_counts")
        return summary_df
    
    def create_regional_comparison_diagrams(self, df: pd.DataFrame) -> Dict[str, Any]:
        """Create Venn diagrams for comparisons within each group-region combination, arranging all in 3x2 grids per page."""
        group_region_combos = df[['group', 'region']].drop_duplicates()
        detailed_results = {}
        summary_data = []
        pdf_path = os.path.join(self.output_manager.dirs['pathway_plots'], "regional_comparison_diagrams.pdf")
        all_venn_plots = []  # List of (title, plotting function)
        for _, row in group_region_combos.iterrows():
            group, region = row['group'], row['region']
            subset = df[(df['group'] == group) & (df['region'] == region)]
            comparisons = subset['comparison'].unique()
            combo_key = f"{group}_{region}"
            detailed_results[combo_key] = {
                'group': group,
                'region': region,
                'comparisons': list(comparisons),
                'comparison_terms': {},
                'intersections': {},
                'unique_terms': {}
            }
            comparison_terms = {}
            for comparison in comparisons:
                terms = set(subset[subset['comparison'] == comparison]['Term'])
                comparison_terms[comparison] = terms
                detailed_results[combo_key]['comparison_terms'][comparison] = list(terms)
                summary_data.append({
                    'group': group,
                    'region': region,
                    'comparison': comparison,
                    'term_count': len(terms)
                })
            # Collect all venn plotting functions for this group-region
            venn_jobs = []
            if len(comparisons) < 2:
                def plot_empty(ax, title=f'{group} - {region}'):
                    ax.text(0.5, 0.5, f'Only {len(comparisons)} comparison\nin {group}-{region}', 
                            ha='center', va='center', transform=ax.transAxes, fontsize=12)
                    ax.set_title(title, fontsize=14)
                    ax.axis('off')
                venn_jobs.append((f'{group} - {region}', plot_empty))
            elif len(comparisons) <= 3:
                def plot_venn(ax, cterms=comparison_terms, dres=detailed_results[combo_key], comps=list(comparisons), g=group, r=region):
                    self._create_comparison_venn(ax, cterms, dres)
                    ax.set_title(f'{g} - {r}', fontsize=14, pad=20)
                venn_jobs.append((f'{group} - {region}', plot_venn))
            else:
                comp_triplets = list(itertools.combinations(comparisons, 3))
                for triplet in comp_triplets:
                    def plot_triplet(ax, triplet=triplet, cterms=comparison_terms, g=group, r=region):
                        venn3([cterms[triplet[0]], cterms[triplet[1]], cterms[triplet[2]]], set_labels=triplet, ax=ax)
                        ax.set_title(f"{g}-{r}: {triplet[0]}, {triplet[1]}, {triplet[2]}", fontsize=12)
                    venn_jobs.append((f"{group}-{region}: {triplet[0]}, {triplet[1]}, {triplet[2]}", plot_triplet))
            all_venn_plots.extend(venn_jobs)
        # Now arrange all venn plots into 3x2 grids per page
        per_page = 6
        n_pages = (len(all_venn_plots) + per_page - 1) // per_page
        with PdfPages(pdf_path) as pdf:
            for page in range(n_pages):
                fig, axes = plt.subplots(3, 2, figsize=(16, 18))
                axes = axes.flatten()
                for i in range(per_page):
                    idx = page * per_page + i
                    if idx >= len(all_venn_plots):
                        axes[i].set_visible(False)
                        continue
                    title, plot_func = all_venn_plots[idx]
                    plot_func(axes[i])
                plt.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)
        # Save summary data
        summary_df = pd.DataFrame(summary_data)
        self.output_manager.save_csv(summary_df, "comparison_term_counts")
        return detailed_results
    
    def _create_comparison_venn(self, ax, comparison_terms: Dict[str, Set], results: Dict[str, Any]):
        """Create individual Venn diagram for comparisons or a stacked bar chart for >3 comparisons."""
        comparisons = list(comparison_terms.keys())
        if len(comparisons) == 2:
            comp_a, comp_b = comparisons[0], comparisons[1]
            intersection = comparison_terms[comp_a] & comparison_terms[comp_b]
            unique_a = comparison_terms[comp_a] - comparison_terms[comp_b]
            unique_b = comparison_terms[comp_b] - comparison_terms[comp_a]
            
            # Store results
            results['intersections'][f"{comp_a} ∩ {comp_b}"] = list(intersection)
            results['unique_terms'][f"{comp_a} only"] = list(unique_a)
            results['unique_terms'][f"{comp_b} only"] = list(unique_b)
            
            # Use full labels, removing previous truncation
            venn = venn2([comparison_terms[comp_a], comparison_terms[comp_b]], 
                        set_labels=comparisons, ax=ax)
            
            # Add counts
            if venn.get_label_by_id('10'):
                venn.get_label_by_id('10').set_text(f'{len(unique_a)}')
            if venn.get_label_by_id('01'):
                venn.get_label_by_id('01').set_text(f'{len(unique_b)}')
            if venn.get_label_by_id('11'):
                venn.get_label_by_id('11').set_text(f'{len(intersection)}')
        
        elif len(comparisons) == 3:
            comp_a, comp_b, comp_c = comparisons[0], comparisons[1], comparisons[2]
            
            # Calculate intersections
            all_three = comparison_terms[comp_a] & comparison_terms[comp_b] & comparison_terms[comp_c]
            ab_only = (comparison_terms[comp_a] & comparison_terms[comp_b]) - comparison_terms[comp_c]
            ac_only = (comparison_terms[comp_a] & comparison_terms[comp_c]) - comparison_terms[comp_b]
            bc_only = (comparison_terms[comp_b] & comparison_terms[comp_c]) - comparison_terms[comp_a]
            a_only = comparison_terms[comp_a] - comparison_terms[comp_b] - comparison_terms[comp_c]
            b_only = comparison_terms[comp_b] - comparison_terms[comp_a] - comparison_terms[comp_c]
            c_only = comparison_terms[comp_c] - comparison_terms[comp_a] - comparison_terms[comp_b]
            
            # Store results
            results['intersections'][f"{comp_a} ∩ {comp_b} ∩ {comp_c}"] = list(all_three)
            results['intersections'][f"{comp_a} ∩ {comp_b} only"] = list(ab_only)
            results['intersections'][f"{comp_a} ∩ {comp_c} only"] = list(ac_only)
            results['intersections'][f"{comp_b} ∩ {comp_c} only"] = list(bc_only)
            results['unique_terms'][f"{comp_a} only"] = list(a_only)
            results['unique_terms'][f"{comp_b} only"] = list(b_only)
            results['unique_terms'][f"{comp_c} only"] = list(c_only)
            
            # Use full labels, removing previous truncation
            venn = venn3([comparison_terms[comp_a], comparison_terms[comp_b], comparison_terms[comp_c]], 
                        set_labels=comparisons, ax=ax)
            
            # Add counts to venn diagram
            venn_labels = ['100', '010', '001', '110', '101', '011', '111']
            counts = [len(a_only), len(b_only), len(c_only), len(ab_only), len(ac_only), len(bc_only), len(all_three)]
            
            for label, count in zip(venn_labels, counts):
                if venn.get_label_by_id(label):
                    venn.get_label_by_id(label).set_text(f'{count}')
        
        else:
            # For >3, handled in create_regional_comparison_diagrams with triplet Venns
            ax.axis('off')

class SharedTermsAnalyzer:
    """Analyzes shared terms across different groupings"""
    
    def __init__(self, output_manager: OutputManager):
        self.output_manager = output_manager
    
    def analyze_shared_terms_by_group(self, df: pd.DataFrame) -> Dict[str, Any]:
        """Analyze shared terms across regions for each group"""
        groups = df['group'].unique()
        shared_results = {}
        summary_data = []
        
        for group in groups:
            self.output_manager.log_info(f"Analyzing shared terms for group: {group}")
            
            group_data = df[df['group'] == group]
            regions = group_data['region'].unique()
            
            # Get terms for each region
            region_terms = {}
            for region in regions:
                region_terms[region] = set(group_data[group_data['region'] == region]['Term'])
            
            group_results = {
                'region_terms': region_terms,
                'pairwise_intersections': {},
                'all_shared': set(),
                'unique_to_each': {}
            }
            
            if len(regions) >= 2:
                # Pairwise intersections
                for i, region1 in enumerate(regions):
                    for region2 in regions[i+1:]:
                        shared = region_terms[region1] & region_terms[region2]
                        pair_key = f"{region1} ∩ {region2}"
                        group_results['pairwise_intersections'][pair_key] = shared
                        
                        summary_data.append({
                            'group': group,
                            'analysis_type': 'pairwise_shared',
                            'category': pair_key,
                            'count': len(shared),
                            'terms': '; '.join(sorted(shared)) if shared else 'None'
                        })
                
                # All regions intersection
                if len(regions) > 2:
                    all_shared = set.intersection(*region_terms.values())
                    group_results['all_shared'] = all_shared
                    
                    summary_data.append({
                        'group': group,
                        'analysis_type': 'all_regions_shared',
                        'category': 'all_regions',
                        'count': len(all_shared),
                        'terms': '; '.join(sorted(all_shared)) if all_shared else 'None'
                    })
                
                # Unique terms for each region
                for region in regions:
                    other_regions = [r for r in regions if r != region]
                    if other_regions:
                        other_terms = set.union(*[region_terms[r] for r in other_regions])
                        unique_terms = region_terms[region] - other_terms
                        group_results['unique_to_each'][region] = unique_terms
                        
                        summary_data.append({
                            'group': group,
                            'analysis_type': 'unique_to_region',
                            'category': region,
                            'count': len(unique_terms),
                            'terms': '; '.join(sorted(unique_terms)) if unique_terms else 'None'
                        })
            
            # Add region totals
            for region, terms in region_terms.items():
                summary_data.append({
                    'group': group,
                    'analysis_type': 'region_total',
                    'category': region,
                    'count': len(terms),
                    'terms': '; '.join(sorted(terms))
                })
            
            shared_results[group] = group_results
        
        # Save summary
        summary_df = pd.DataFrame(summary_data)
        self.output_manager.save_csv(summary_df, "shared_terms_analysis_summary")
        
        return shared_results
    
    def create_summary_heatmap(self, df: pd.DataFrame):
        """Create summary heatmap of terms by group and region"""
        summary = df.groupby(['group', 'region'])['Term'].nunique().reset_index()
        pivot_summary = summary.pivot(index='group', columns='region', values='Term')
        
        mask = pivot_summary.isnull()
        
        fig, ax = plt.subplots(figsize=(10, 6))
        sns.heatmap(pivot_summary, annot=True, fmt='.0f', cmap='Blues', 
                    mask=mask, cbar_kws={'label': 'Number of Unique Terms'},
                    linewidths=0.5, ax=ax)
        
        ax.set_title('Number of Unique Terms by Group and Region\n(Gray cells = no data)')
        ax.set_xlabel('Region')
        ax.set_ylabel('Group')
        
        self.output_manager.save_plot(fig, "summary_heatmap_group_region")
        
        # Save the pivot table as CSV
        self.output_manager.save_csv(pivot_summary.reset_index(), "summary_table_group_region")
        
        return pivot_summary


class PathwayAnalyzer:
    """Main class orchestrating the entire pathway analysis workflow"""
    
    def __init__(self, deg_clusters_dir, pathway_results_dir, pathway_plots_dir):
        self.output_manager = OutputManager(deg_clusters_dir, pathway_results_dir, pathway_plots_dir)
        self.enrichr_analyzer = EnrichrAnalyzer(self.output_manager)
        self.venn_analyzer = VennAnalyzer(self.output_manager)
        self.shared_terms_analyzer = SharedTermsAnalyzer(self.output_manager)
        
    def run_complete_analysis(self, input_file: str = None):
        """Run the complete pathway analysis workflow"""
        self.output_manager.log_info("Starting complete pathway analysis workflow")

        # Check if merged pathway results already exist
        merged_pathway_filename = "merged_pathway_results.csv"
        merged_pathway_path = os.path.join(
            self.output_manager.dirs['pathway_results'],
            merged_pathway_filename
        )
        if os.path.exists(merged_pathway_path):
            self.output_manager.log_info(f"Merged pathway results found: {merged_pathway_path}. Loading directly.")
            all_pathways = pd.read_csv(merged_pathway_path)
        else:
            try:
                # Load data
                if input_file is None:
                    input_file = os.path.join(
                        self.output_manager.dirs['deg_clusters'], 
                        'merged_dge_on_clusters.csv'
                    )
                
                df = self._load_and_preprocess_data(input_file)
                
                # Run enrichment analysis
                self.output_manager.log_info("Step 1: Running enrichment analysis")
                enrichment_results = self.enrichr_analyzer.run_enrichment_analysis(df)
                
                if not enrichment_results:
                    self.output_manager.log_info("No significant pathways found. Exiting.")
                    return
                
                # Get merged pathway results
                all_pathways = self.enrichr_analyzer.get_merged_results()
                
                # Save main results
                self.output_manager.save_csv(all_pathways, "merged_pathway_results")

            except Exception as e:
                self.output_manager.log_info(f"Error in analysis: {str(e)}")
                raise
            
        # Log summary statistics
        pathway_subset = all_pathways[['Term', 'group', 'region', 'comparison', 'Genes']]
        combination_counts = pathway_subset.groupby(['group', 'region', 'comparison'])['Term'].nunique().reset_index()
        self.output_manager.save_csv(combination_counts, "pathway_counts_by_combination")
        
        # Step 2: Venn diagram analyses
        self.output_manager.log_info("Step 2: Creating Venn diagram analyses")
        
        # Region-based Venn diagrams
        region_summary = self.venn_analyzer.create_sample_comparison_diagrams(pathway_subset)
        
        # Comparison-based Venn diagrams
        comparison_details = self.venn_analyzer.create_regional_comparison_diagrams(pathway_subset)
        
        # Step 3: Shared terms analysis
        self.output_manager.log_info("Step 3: Analyzing shared terms")
        shared_results = self.shared_terms_analyzer.analyze_shared_terms_by_group(pathway_subset)
        
        # Step 4: Summary visualizations
        self.output_manager.log_info("Step 4: Creating summary visualizations")
        summary_heatmap = self.shared_terms_analyzer.create_summary_heatmap(pathway_subset)
        
        # Final summary
        self._log_final_summary(all_pathways, shared_results)
        
        self.output_manager.log_info("Analysis completed successfully")
            
    
    def _load_and_preprocess_data(self, input_file: str) -> pd.DataFrame:
        """Load and preprocess the DGE data"""
        self.output_manager.log_info(f"Loading data from: {input_file}")
        
        if not os.path.exists(input_file):
            raise FileNotFoundError(f"Input file not found: {input_file}")
        
        df = pd.read_csv(input_file)
        # Detect and rename columns if needed
        colset = set(df.columns)
        # For merged_dge_on_cellTypes.csv
        if {"names", "scores", "pvals", "pvals_adj", "logfoldchanges", "group", "cell_type", "comparison"}.issubset(colset):
            df = df.rename(columns={
                "names": "gene",
                "pvals": "p_val",
                "pvals_adj": "p_val_adj",
                "logfoldchanges": "avg_log2FC",
                "cell_type": "region"  # Treat cell_type as region for downstream compatibility
            })

        df['gene'] = df['gene'].str.replace('.', '-', regex=False)
        
        self.output_manager.log_info(f"Loaded {len(df)} rows of DGE data")
        
        # Log data summary
        summary_stats = {
            'total_rows': len(df),
            'unique_genes': df['gene'].nunique(),
            'unique_groups': df['group'].nunique(),
            'unique_regions': df['region'].nunique(),
            'unique_comparisons': df['comparison'].nunique()
        }
        
        self.output_manager.log_analysis_summary("Data Loading Summary", summary_stats)
        
        return df
    
    def _log_final_summary(self, all_pathways: pd.DataFrame, shared_results: Dict[str, Any]):
        """Log final analysis summary"""
        summary_stats = {
            'total_significant_pathways': len(all_pathways),
            'unique_pathway_terms': all_pathways['Term'].nunique(),
            'groups_analyzed': len(shared_results),
            'total_group_region_combinations': len(all_pathways[['group', 'region']].drop_duplicates())
        }
        
        self.output_manager.log_analysis_summary("Final Analysis Summary", summary_stats)


def main():
    """Main execution function"""
    # Load YAML config for output_dirs
    with open("config/batch_config.yaml", "r") as f:
        config = yaml.safe_load(f)
    output_dirs = config["output_dirs"]

    # List of directories to search
    search_dirs = [
        output_dirs["deg_clusters"], 
        output_dirs["deg_celltypes"], 
        output_dirs.get("deg_combine_celltypes", "")
    ]
    pathway_results_dirs = [
        output_dirs["pathway_analysis_clusters"], 
        output_dirs["pathway_analysis_celltypes"], 
        output_dirs.get("pathway_analysis_combine_celltypes", "")
    ]
    pathway_plots_dirs = [
        output_dirs.get("pathway_analysis_clusters_plots", ""), 
        output_dirs.get("pathway_analysis_celltypes_plots", ""),
        output_dirs.get("pathway_analysis_combine_celltypes_plots", "")
    ]

    for search_dir, pathway_results_dir, pathway_plots_dir in zip(search_dirs, pathway_results_dirs, pathway_plots_dirs):
        pattern = os.path.join(search_dir, "merged_dge_*.csv")
        for input_file in glob.glob(pattern):
            print(f"Processing input file: {input_file}")
            # Ensure output directories exist
            os.makedirs(pathway_results_dir, exist_ok=True)
            os.makedirs(pathway_plots_dir, exist_ok=True)
            # Log the directories being used
            print(f"Using pathway results directory: {pathway_results_dir}")
            print(f"Using pathway plots directory: {pathway_plots_dir}")
            
            # Initialize analyzer with YAML-based output dirs
            analyzer = PathwayAnalyzer(search_dir, pathway_results_dir, pathway_plots_dir)
            analyzer.run_complete_analysis(input_file=input_file)

if __name__ == "__main__":
    main()