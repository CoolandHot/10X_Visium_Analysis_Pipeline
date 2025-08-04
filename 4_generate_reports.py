#!/usr/bin/env python3
"""
Spatial Transcriptomics Analysis Report Generator
Creates standalone HTML5 reports for different analysis sections
"""

from pathlib import Path
import yaml
import csv
import html
import fitz  # PyMuPDF
import os
import shutil

class ReportGenerator:
    """Generates interactive HTML5 reports for spatial transcriptomics analysis"""
    
    def __init__(self, config_path='config/batch_config.yaml'):
        """Initialize with configuration"""
        self.base_dir = Path(__file__).parent
        self.config = self._load_config(config_path)
        self.project_dir = Path(self.config.get('project_dir', self.base_dir))
        # Use cellbrowser_html_output_dir from config if available
        self.cellbrowser_html_output_dir = Path(self.config.get('cellbrowser_html_output_dir', "output/html_reports/cellbrowser"))
        self.report_dir = self.project_dir / "output" / "html_reports"
        self.report_dir.mkdir(exist_ok=True)
        
        # Extract cluster IDs for DGE intro
        dge = self.config.get('Differential_Gene_Analysis', {})
        self.inTumour_clusters = dge.get('inTumour_cluster_nums', [])
        self.outTumour_clusters = dge.get('outTumour_cluster_nums', [])
        self.edgeTumour_clusters = dge.get('edgeTumour_cluster_nums', [])
        
        self.output_dirs = self.config.get('output_dirs', {})
        self.default_dirs = {
            'clustering_umap_spatial': 'output/clustering/clusters_umap_spatial/',
            'clustering_markers': 'output/clustering/markers/',
            'clustering_spatialAware': 'output/clustering/spatial_aware_clustering/',
            'cell_type_spacexr': 'output/cell_type_prediction/spaceXR_results/',
            'cell_type_cell2loc': 'output/cell_type_prediction/cell2location_results/',
            'deg_celltypes': 'output/differential_expression/cell_type_comparisons',
            'deg_combine_celltypes': 'output/differential_expression/combine_cell_type_comparisons',
            'deg_clusters': 'output/differential_expression/cluster_comparisons',
            'pathway_analysis_clusters': 'output/pathway_analysis_clusters/',
            'pathway_analysis_clusters_plots': 'output/pathway_analysis_clusters/plots/',
            'pathway_analysis_celltypes': 'output/pathway_analysis_celltypes/',
            'pathway_analysis_celltypes_plots': 'output/pathway_analysis_celltypes/plots/',
            'pathway_analysis_combine_celltypes': 'output/pathway_analysis_combine_celltypes/',
            'pathway_analysis_combine_celltypes_plots': 'output/pathway_analysis_combine_celltypes/plots/',
            'tfa_activities': 'output/transcription_factor_activity/TF_activities/',
            'tfa_features': 'output/transcription_factor_activity/feature_importance/',
            'tfa_heatmaps': 'output/transcription_factor_activity/heatmaps/',
        }
    
    def _get_dir(self, key):
        """Get output dir from config or default"""
        return self.output_dirs.get(key, self.default_dirs[key])

    def _load_config(self, config_path):
        """Load configuration file"""
        try:
            config_file = Path(config_path)
            if not config_file.is_absolute():
                config_file = self.base_dir / config_path
            with open(config_file, 'r') as f:
                return yaml.safe_load(f)
        except FileNotFoundError:
            print(f"Config file not found: {config_path}")
            return {}
    
    def _copy_visualization_tools(self):
        """Copy visualization tools to html_reports directory if not exists, excluding specific files."""
        viz_tools_src = self.project_dir.parent / "visualisation_tools"
        viz_tools_dest = self.report_dir / "visualisation_tools"
        
        if viz_tools_src.exists() and not viz_tools_dest.exists():
            try:
                # Define files and directories to exclude, relative to viz_tools_src
                exclusions = {
                    'pathway/pathway.html',
                    'pathway/js/pathway.js',
                    'pathway/data'
                }

                def ignore_func(directory, contents):
                    ignored_items = set()
                    for item in contents:
                        # Create a path relative to the source of the copy operation
                        relative_path = (Path(directory).relative_to(viz_tools_src) / item).as_posix()
                        if relative_path in exclusions:
                            ignored_items.add(item)
                    return ignored_items

                shutil.copytree(viz_tools_src, viz_tools_dest, ignore=ignore_func)
                print(f"Copied visualization tools from {viz_tools_src} to {viz_tools_dest}")
            except Exception as e:
                print(f"Warning: Could not copy visualization tools: {e}")
        elif viz_tools_dest.exists():
            print(f"Visualization tools already exist at {viz_tools_dest}")
        else:
            print(f"Warning: Source visualization tools directory not found at {viz_tools_src}")

    def generate_all_reports(self):
        """Generate all reports"""
        print("Generating spatial transcriptomics analysis reports...")
        
        # Copy visualization tools first
        self._copy_visualization_tools()
        
        # Generate main navigation
        self._generate_main_navigation()
        
        # Generate individual section reports
        self._generate_section_reports()
        
        print(f"Reports generated in: {self.report_dir}")
        print(f"Open {self.report_dir}/index.html to view the main navigation")

    def _generate_section_reports(self):
        """Generate all reports"""
        sections = [
            {
                "filename": "clustering_report.html",
                "title": "Spatial Clustering Analysis",
                "intro": """
                Spatial clustering analysis identifies distinct cellular neighborhoods and tissue regions based on gene expression patterns.
                We employ multiple clustering algorithms including Seurat's graph-based clustering, BANKSY for spatial-aware clustering,
                and BayesSpace for enhanced spatial clustering. Each method provides complementary insights into tissue organization.
                """,
                "file_sections": self._clustering_file_sections()
            },
            {
                "filename": "cell_type_report.html",
                "title": "Cell Type Deconvolution Analysis",
                "intro": """
                Cell type deconvolution estimates the cellular composition at each spatial location using reference single-cell datasets.
                We use Cell2Location for probabilistic cell type mapping and SpaceXR separately.
                These methods provide complementary views of cellular heterogeneity across tissue regions.
                """,
                "file_sections": self._cell_type_file_sections()
            },
            {
                "filename": "dge_report.html",
                "title": "Differential Gene Expression Analysis",
                "intro": f"""
                Differential gene expression analysis identifies genes with significant expression changes between conditions.
                We perform both cross-sample comparisons (between treatment groups) and within-sample analysis 
                (comparing high vs low cell type abundance regions). Statistical testing uses appropriate methods 
                for spatial transcriptomics data characteristics.<br>
                
                <b>Cluster ID indications:</b><br>
                <ul>
                  <li><b>inTumour clusters</b>: {self.inTumour_clusters}</li>
                  <li><b>outTumour clusters</b>: {self.outTumour_clusters}</li>
                  <li><b>edgeTumour clusters</b>: {self.edgeTumour_clusters}</li>
                </ul>
                """,
                "file_sections": self._dge_file_sections(),
                "is_dge": True
            },
            {
                "filename": "pathway_report.html",
                "title": "Pathway Enrichment Analysis",
                "intro": """
                Pathway analysis interprets differential expression results by identifying enriched biological pathways and processes.
                We use KEGG database to comprehensively characterize functional changes
                associated with spatial regions and treatment conditions. There are other datasets available, such as GO and Reactome, but we use KEGG for its wide acceptance and avoid creating too much results at a time.
                """,
                "file_sections": self._pathway_file_sections()
            },
            {
                "filename": "tfa_report.html",
                "title": "Transcription Factor Activity Analysis",
                "intro": """
                Transcription factor activity analysis infers regulatory network activity from gene expression data.
                We estimate TF activities using DOROTHEA regulons and identify spatially variable transcription factors
                that may drive regional gene expression differences in brain tumor tissue.
                """,
                "file_sections": self._tfa_file_sections()
            }
        ]
        for section in sections:
            if section.get("is_dge"):
                self._generate_dge_report(
                    section["filename"],
                    section["title"],
                    section["file_sections"],
                    intro_text=section["intro"]
                )
            else:
                self._generate_section_report(
                    section["filename"],
                    section["title"],
                    section["file_sections"],
                    intro_text=section["intro"]
                )

    def _clustering_file_sections(self):
        return [
            {
                "title": "Interactive Visualizations",
                "description": "Interactive HTML plots for cluster exploration",
                "pattern": "**/*.html",
                "subdirs": [self._get_dir('clustering_spatialAware')],
                "exclude_subdirs": [],
                "exclude_patterns": []
            },
            {
                "title": "Spatial-Aware Clustering",
                "description": "Clustering results using BANKSY and BayesSpace methods. For files starting with 'Cell2Location_*', they are clustered based the inferred cell type abundance from Cell2Location, while the rest are clustered based the gene expression.",
                "pattern": "**/*.pdf",
                "subdirs": [self._get_dir('clustering_spatialAware')],
                "exclude_subdirs": [],
                "exclude_patterns": []
            },
            {
                "title": "SNN Clustering (Seurat GNN, traditional leiden method)",
                "description": "2D/3D UMAP plots and spatial cluster using Seurat GNN algorithm visualizations for each sample",
                "pattern": "**/*.pdf",
                "subdirs": [self._get_dir('clustering_umap_spatial')],
                "exclude_subdirs": [],
                "exclude_patterns": []
            },
            {
                "title": "Marker Gene by SNN Clustering", 
                "description": "Top marker genes identified for each cluster",
                "pattern": "**/*markers*.csv",
                "subdirs": [self._get_dir('clustering_markers')],
                "exclude_subdirs": [],
                "exclude_patterns": []
            },
        ]

    def _cell_type_file_sections(self):
        return [
            {
                "title": "Cell2Location Predictions",
                "description": "Cell type abundance maps and spatial distributions from Cell2Location analysis.\nFor the index.html of cellbrowser, you need to open it with a web server backend.",
                "pattern": "**/*.html",
                "subdirs": ['output/differential_expression/cell_type_comparisons', str(self.cellbrowser_html_output_dir)],
                "exclude_subdirs": [],
                "exclude_patterns": []
            },
            {
                "title": "Cell Abundance Data",
                "description": "Quantitative cell type abundance matrices and cluster assignments",
                "pattern": "**/*.csv",
                "subdirs": [
                    self._get_dir('cell_type_cell2loc'),
                    self._get_dir('cell_type_spacexr')
                ],
                "exclude_subdirs": ["output/cell_type_prediction/cell2location_results/nmf_analysis"],
                "exclude_patterns": ["**/n_fact*.csv"]
            },
            {
                "title": "SpaceXR Predictions",
                "description": "Cell type predictions using SpaceXR reference mapping",
                "pattern": "**/*.pdf", 
                "subdirs": [self._get_dir('cell_type_spacexr')],
                "exclude_subdirs": [],
                "exclude_patterns": []
            },
        ]

    def _create_dge_subsections(self, base_dir_key, data_name, 
                                data_within_desc, data_cross_desc,
                                viz_within_desc, viz_cross_desc,
                                is_combined=False):
        base_dir = self._get_dir(base_dir_key)
        
        combined_title_suffix = " (on combined cell types)" if is_combined else ""
        combined_desc_suffix = "\n'Combined cell types' means we combine some of the subtypes into one type, such as 'CD4_T_cell' instead of ['Effector.Memory.CD4.T', 'Naive.CD4.T', 'Stemlike.CD4.T', 'Regulatory.CD4.T', 'Lineage.neg.CD4.T'], etc." if is_combined else ""

        he_viz_context_cell = "We provide the visualizations on H&E images, which are the original images of the tissue sections, to help you understand the spatial context of the cell abundance."
        he_viz_context_gene = "They are  visualizations on H&E images, which are the original images of the tissue sections, to help you understand the spatial context of the cluster regions."
        
        if data_name == "Gene Expression":
            within_viz_he_desc = f"{viz_within_desc} within samples, we also separate the edge regions but they are not included in the comparisons. {he_viz_context_gene}"
            cross_viz_he_desc = f"{viz_cross_desc}, etc. {he_viz_context_gene}"
        else:
            within_viz_he_desc = f"In ADI, {viz_within_desc}, etc. {he_viz_context_cell}"
            cross_viz_he_desc = f"{viz_cross_desc}, etc."

        return [
            # Data
            {
                "title": f"Within-/Cross-Sample {data_name} {'Combined' if is_combined else ''} Data{combined_title_suffix}",
                "description": f"Differential comparison between sample groups ({data_cross_desc}, etc.); within-sample comparisons ({data_within_desc}, etc.){combined_desc_suffix}",
                "pattern": "**/merged_*.csv",
                "subdirs": [base_dir], "exclude_subdirs": [], "exclude_patterns": []
            },
            # Within-Sample
            {
                "title": f"Within-Sample {data_name} {'Combined' if is_combined else ''} Visualizations on H&E images",
                "description": within_viz_he_desc,
                "pattern": "**/*.pdf",
                "subdirs": [os.path.join(base_dir, 'dge_results_within_sample')], "exclude_subdirs": [], "exclude_patterns": []
            },
            {
                "title": f"Within-Sample {data_name} {'Combined' if is_combined else ''} Inspective Visualizations",
                "description": f"In ADI, {viz_within_desc}, etc. Volcano plots, boxplots, and scatter plots to inspect the differential expression results.",
                "pattern": "**/*.pdf",
                "subdirs": [os.path.join(base_dir, 'inspection_plots_pdf', 'within_sample')], "exclude_subdirs": [], "exclude_patterns": []
            },
            {
                "title": f"Within-Sample {data_name} {'Combined' if is_combined else ''} Inspective Visualizations Data",
                "description": f"In ADI, {viz_within_desc}, etc. The data in csv format for the inspective visualizations.",
                "pattern": "**/*.csv",
                "subdirs": [os.path.join(base_dir, 'inspection_plot_data', 'within_sample')], "exclude_subdirs": [], "exclude_patterns": []
            },
            # Cross-Sample
            {
                "title": f"Cross-Sample {data_name} {'Combined' if is_combined else ''} Visualizations on H&E images",
                "description": cross_viz_he_desc,
                "pattern": "**/*.pdf",
                "subdirs": [os.path.join(base_dir, 'dge_results_across_sample')], "exclude_subdirs": [], "exclude_patterns": []
            },
            {
                "title": f"Cross-Sample {data_name} {'Combined' if is_combined else ''} Inspective Visualizations",
                "description": f"{viz_cross_desc}, etc. Volcano plots, boxplots, and scatter plots to inspect the differential expression results.",
                "pattern": "**/*.pdf",
                "subdirs": [os.path.join(base_dir, 'inspection_plots_pdf', 'across_sample')], "exclude_subdirs": [], "exclude_patterns": []
            },
            {
                "title": f"Cross-Sample {data_name} {'Combined' if is_combined else ''} Inspective Visualizations Data",
                "description": f"{viz_cross_desc}, etc. The data in csv format for the inspective visualizations.",
                "pattern": "**/*.csv",
                "subdirs": [os.path.join(base_dir, 'inspection_plot_data', 'across_sample')], "exclude_subdirs": [], "exclude_patterns": []
            }
        ]

    def _dge_file_sections(self):
        
        sections = []
        
        # Cell Abundance
        cell_abundance_args = {
            "base_dir_key": 'deg_celltypes', "data_name": "Cell Abundance",
            "data_within_desc": "regions with high vs low abundance of T cells in ADI",
            "data_cross_desc": "regions that expressed higher population of T cells in ADI vs those regions in SAL",
            "viz_within_desc": "higher quantile of tumour cell abundance density vs lower quantile",
            "viz_cross_desc": "higher quantile of tumour cell abundance density in ADI vs higher quantile of tumour cell abundance density in SAL"
        }
        sections.extend(self._create_dge_subsections(**cell_abundance_args))
        
        # Combined Cell Abundance
        combined_cell_abundance_args = cell_abundance_args.copy()
        combined_cell_abundance_args.update({
            "base_dir_key": 'deg_combine_celltypes',
            "is_combined": True
        })
        sections.extend(self._create_dge_subsections(**combined_cell_abundance_args))
        
        # Gene Expression
        gene_expr_args = {
            "base_dir_key": 'deg_clusters', "data_name": "Gene Expression",
            "data_within_desc": "tumour regions vs non-tumour regions in ADI",
            "data_cross_desc": "tumour regions in ADI vs tumour regions in SAL",
            "viz_within_desc": "Tumour cells vs non-tumour cells comparisons",
            "viz_cross_desc": "The tumour cell clusters/regions in ADI vs those regions in SAL"
        }
        sections.extend(self._create_dge_subsections(**gene_expr_args))
        
        return sections

    def _pathway_file_sections(self):
        section_infos = [
            {
                "title": "Pathway Enrichment Results by Clusters",
                "description": "Pathway enrichment analysis results from differential expression on gene expression of spatial clusters",
                "data_dir": self._get_dir('pathway_analysis_clusters'),
                "plot_dir": self._get_dir('pathway_analysis_clusters_plots')

            },
            {
                "title": "Pathway Enrichment Results by Cell Types",
                "description": "Pathway enrichment analysis results from differential expression on cell abundance of cell types",
                "data_dir": self._get_dir('pathway_analysis_celltypes'),
                "plot_dir": self._get_dir('pathway_analysis_celltypes_plots')
            },
            {
                "title": "Pathway Enrichment Results by Combined Cell Types",
                "description": "Pathway enrichment analysis results from differential expression on cell abundance of combined cell types.\nThis is the combined results of the cell types, such as 'CD4_T_cell' instead of ['Effector.Memory.CD4.T', 'Naive.CD4.T', 'Stemlike.CD4.T', 'Regulatory.CD4.T', 'Lineage.neg.CD4.T'], etc.",
                "data_dir": self._get_dir('pathway_analysis_combine_celltypes'),
                "plot_dir": self._get_dir('pathway_analysis_combine_celltypes_plots')
            }
        ]

        sections = []
        for info in section_infos:
            sections.extend([
                {
                    "title": info["title"],
                    "description": info["description"],
                    "pattern": "**/*.csv",
                    "subdirs": [info["data_dir"]],
                    "exclude_subdirs": [],
                    "exclude_patterns": []
                },
                {
                    "title": f"{info['title']} Visualizations",
                    "description": "Pathway network plots, enrichment plots, and summary visualizations",
                    "pattern": "**/*.pdf",
                    "subdirs": [info["plot_dir"]],
                    "exclude_subdirs": [],
                    "exclude_patterns": []
                }
            ])            

        return sections

    def _tfa_file_sections(self):
        return [
            {
                "title": "TF Activity Matrices",
                "description": "Transcription factor activity scores across spatial locations",
                "pattern": "**/*TF_activities*.csv",
                "subdirs": [self._get_dir('tfa_activities')],
                "exclude_subdirs": [],
                "exclude_patterns": []
            },
            {
                "title": "Important TF Features",
                "description": "Top transcription factors identified by machine learning feature selection",
                "pattern": "**/*top_features*.csv", 
                "subdirs": [self._get_dir('tfa_features')],
                "exclude_subdirs": [],
                "exclude_patterns": []
            },
            {
                "title": "TF Activity Heatmaps",
                "description": "Spatial heatmaps showing transcription factor activity patterns",
                "pattern": "**/*heatmap*.pdf",
                "subdirs": [self._get_dir('tfa_heatmaps')],
                "exclude_subdirs": [],
                "exclude_patterns": []
            }
        ]
    
    def _generate_main_navigation(self):
        """Generate main navigation page"""
        # Get CSS and JS first
        extra_css = self._get_navigation_css()
        extra_js = ""
        
        # Get base template
        html_content = self._get_base_html_template("Spatial Transcriptomics Analysis Reports")
        
        sections = [
            {
                "title": "Spatial Clustering Analysis",
                "description": "Cell clustering results using various algorithms including Seurat, BANKSY, and BayesSpace",
                "file": "clustering_report.html",
                "icon": "🧬"
            },
            {
                "title": "Cell Type Deconvolution", 
                "description": "Cell type predictions using Cell2Location and SpaceXR methods",
                "file": "cell_type_report.html",
                "icon": "🔬"
            },
            {
                "title": "Differential Gene Expression",
                "description": "Cross-sample and within-sample differential expression analysis",
                "file": "dge_report.html", 
                "icon": "📊"
            },
            {
                "title": "Pathway Analysis",
                "description": "Functional pathway enrichment analysis from differential expression results",
                "file": "pathway_report.html",
                "icon": "🧭"
            },
            {
                "title": "Transcription Factor Activity",
                "description": "TF activity analysis and regulatory network inference",
                "file": "tfa_report.html",
                "icon": "🎯"
            },
            {
                "title": "Venn Diagram Tool",
                "description": "Interactive Venn diagram visualization for comparing shared/unique gene sets and pathways",
                "file": "visualisation_tools/venn/venn.html",
                "icon": "⭕"
            },
            {
                "title": "Pathway Visualization Tool",
                "description": "Interactive pathway network visualization on multiple pathways by creating node-link graphs and exploration tool\nPlease input the output/pathway_analysis_*/merged_pathway_results.csv file for inspection",
                "file": "visualisation_tools/pathway/pathway_single_csv.html",
                "icon": "🕸️"
            }
        ]
        
        nav_cards = ""
        for section in sections:
            # Check if file exists for visualization tools
            file_exists = True
            if section['file'].startswith('visualisation_tools/'):
                file_path = self.report_dir / section['file']
                file_exists = file_path.exists()
            
            # Add disabled class if file doesn't exist
            card_class = "nav-card"
            if not file_exists:
                card_class += " disabled"
                section['description'] += " (Tool not available - check if visualisation_tools directory exists)"
            
            # Build onclick string safely
            if file_exists:
                onclick_str = f"window.open('{section['file']}', '_blank')"
            else:
                onclick_str = "alert('Visualization tool not found. Please ensure visualisation_tools directory exists.')"
            
            nav_cards += f"""
            <div class="{card_class}" onclick="{onclick_str}">
                <div class="nav-icon">{section['icon']}</div>
                <h3>{section['title']}</h3>
                <p>{section['description']}</p>
                <div class="nav-arrow">→</div>
            </div>
            """
        
        content = f"""
        <div class="main-header">
            <h1>🧠 Spatial Transcriptomics Analysis</h1>
            <p class="subtitle">Brain Tumor Analysis Pipeline - Interactive Reports</p>
        </div>
        
        <div class="nav-grid">
            {nav_cards}
        </div>
        
        <div class="info-section">
            <h2>About This Analysis</h2>
            <p>This comprehensive spatial transcriptomics analysis examines brain tumor samples using 10X Genomics Visium technology. 
            The pipeline includes spatial clustering, cell type deconvolution, differential expression analysis, pathway enrichment, 
            and transcription factor activity inference across multiple sample groups (ADI, RAD, COMB, SAL).</p>
            
            <h3>Visualization Tools</h3>
            <p>Interactive visualization tools are available for exploring results:
            <ul>
                <li><strong>Venn Diagram Tool:</strong> Compare gene sets and pathways with interactive Venn diagrams</li>
                <li><strong>Pathway Visualization Tool:</strong> Explore pathway networks and relationships interactively</li>
            </ul>
            </p>
        </div>
        """
        
        # Replace placeholders in correct order
        html_content = html_content.replace("{{EXTRA_CSS}}", extra_css)
        html_content = html_content.replace("{{EXTRA_JS}}", extra_js)
        html_content = html_content.replace("{{CONTENT}}", content)
        
        with open(self.report_dir / "index.html", 'w') as f:
            f.write(html_content)
    
    def _generate_section_report(self, filename, title, file_sections, intro_text):
        """Generate a section-specific report"""
        # Get CSS and JS first
        extra_css = self._get_section_css() + self._get_floating_nav_css()
        extra_js = self._get_floating_nav_js() + self._get_view_toggle_js() + self._get_search_js()
        
        # Get base template
        html_content = self._get_base_html_template(title)
        
        # Prepare section anchors for floating nav
        section_anchors = []
        sections_html = ""
        for idx, section in enumerate(file_sections):
            anchor = f"section-{idx+1}"
            section_anchors.append({
                "id": anchor,
                "title": section['title']
            })
            files = self._find_files(
                section["pattern"],
                section.get("subdirs", []),
                section.get("exclude_subdirs", []),
                section.get("exclude_patterns", [])
            )
            file_cards = self._generate_file_cards(files)
            
            sections_html += f"""
            <div class="section" id="{anchor}">
                <div class="section-header">
                    <h2>{section['title']}</h2>
                    <div class="search-container">
                        <input type="text" class="search-input" placeholder="🔍 Search files..." data-section="{anchor}">
                        <span class="file-count">({len(files)} files)</span>
                    </div>
                </div>
                <p class="section-description">{section['description']}</p>
                <div class="file-grid" data-section="{section['title'].lower().replace(' ', '-')}">
                    {file_cards}
                </div>
            </div>
            """

        # Floating navigator HTML (now includes view toggle button)
        floating_nav_html = self._generate_floating_nav_html(section_anchors)

        content = f"""
        <div class="header">
            <button onclick="window.open('index.html', '_blank')" class="back-btn">← Back to Main</button>
            <h1>{title}</h1>
        </div>
        {floating_nav_html}
        <div class="intro">
            {intro_text}
        </div>
        <div class="content">
            {sections_html}
        </div>
        """

        # Replace placeholders in correct order
        html_content = html_content.replace("{{EXTRA_CSS}}", extra_css)
        html_content = html_content.replace("{{EXTRA_JS}}", extra_js)
        html_content = html_content.replace("{{CONTENT}}", content)
        
        with open(self.report_dir / filename, 'w') as f:
            f.write(html_content)

    def _generate_dge_report(self, filename, title, file_sections, intro_text):
        """Generate a redesigned DGE report with tabbed interface and dynamic file filtering"""
        # Get CSS and JS first
        extra_css = self._get_section_css() + self._get_floating_nav_css() + self._get_dge_report_css()
        extra_js = self._get_floating_nav_js() + self._get_view_toggle_js() + self._get_search_js() + self._get_dge_report_js()
        
        # Get base template
        html_content = self._get_base_html_template(title)
        
        # Prepare section anchors for floating nav
        section_anchors = []
        sections_html = ""
        for idx, section in enumerate(file_sections):
            anchor = f"section-{idx+1}"
            section_anchors.append({
                "id": anchor,
                "title": section['title']
            })
            files = self._find_files(
                section["pattern"],
                section.get("subdirs", []),
                section.get("exclude_subdirs", []),
                section.get("exclude_patterns", [])
            )
            file_cards = self._generate_file_cards(files)
            
            sections_html += f"""
            <div class="section" id="{anchor}">
                <div class="section-header">
                    <h2>{section['title']}</h2>
                    <div class="search-container">
                        <input type="text" class="search-input" placeholder="🔍 Search files..." data-section="{anchor}">
                        <span class="file-count">({len(files)} files)</span>
                    </div>
                </div>
                <p class="section-description">{section['description']}</p>
                <div class="file-grid" data-section="{section['title'].lower().replace(' ', '-')}">
                    {file_cards}
                </div>
            </div>
            """

        # Floating navigator HTML 
        floating_nav_html = self._generate_floating_nav_html(section_anchors)
        
        # Get batch names from config
        batch_names = self.config.get('batch_names', [])
        
        # Get cell types from cellType_config.yaml
        cell_types = []
        try:
            celltype_config_path = self.base_dir / "config" / "cellType_config.yaml"
            with open(celltype_config_path, 'r') as f:
                celltype_config = yaml.safe_load(f)
                cell_type_combinations = celltype_config.get('shared', {}).get('cell_type_combinations', {})
                cell_types = celltype_config.get('shared', {}).get('target_cell_types', [])

                if cell_type_combinations:
                    other_cell_types = [ct for ct in cell_types if ct not in [
                        item for sublist in cell_type_combinations.values() for item in sublist]]
                    cell_types = list(cell_type_combinations.keys()) + other_cell_types
                    
        except Exception as e:
            print(f"Warning: Could not load cell type config: {e}")
            cell_types = ["T cells", "Macrophages", "Astrocytes"]  # Default fallback

        content = f"""
        <div class="header">
            <button onclick="window.open('index.html', '_blank')" class="back-btn">← Back to Main</button>
            <h1>{title}</h1>
        </div>
        {floating_nav_html}
        <div class="intro">
            {intro_text}
        </div>
        <div class="content">
            <!-- Tab navigation -->
            <div class="tab-navigation">
                <button class="tab-button active" data-tab="cell-type">Cell Type Comparison</button>
                <button class="tab-button" data-tab="gene-expression">Gene Expression Comparison</button>
            </div>
            
            <!-- Tab content containers -->
            <div class="tab-content active" id="cell-type-content">
                <div class="control-panel">
                    <div class="form-group">
                        <label for="cell-comparison-scope">Comparison Scope:</label>
                        <select id="cell-comparison-scope" onchange="toggleVisibility()">
                            <option value="within">Within-sample</option>
                            <option value="cross">Cross-sample</option>
                        </select>
                    </div>
                    
                    <div class="form-group" id="cell-sample-group">
                        <label for="cell-sample-single">Sample:</label>
                        <select id="cell-sample-single">
                            {''.join([f'<option value="{batch}">{batch}</option>' for batch in batch_names])}
                        </select>
                    </div>
                    
                    <div class="form-group hidden" id="cell-sample-pair-group">
                        <label for="cell-sample-a">Sample A:</label>
                        <select id="cell-sample-a">
                            {''.join([f'<option value="{batch}">{batch}</option>' for batch in batch_names])}
                        </select>
                        <label for="cell-sample-b">Sample B:</label>
                        <select id="cell-sample-b">
                            {''.join([f'<option value="{batch}">{batch}</option>' for batch in batch_names])}
                        </select>
                    </div>
                    
                    <div class="form-row">
                        <div class="form-group">
                            <label for="cell-type-1">Cell Type 1:</label>
                            <select id="cell-type-1">
                                {''.join([f'<option value="{ct.replace(" ", "_")}">{ct}</option>' for ct in cell_types])}
                            </select>
                            <div class="quantile-group">
                                <label>Quantile:</label>
                                <label><input type="radio" name="quantile-1" value="high" checked> High</label>
                                <label><input type="radio" name="quantile-1" value="low"> Low</label>
                            </div>
                        </div>
                        
                        <div class="form-group">
                            <label for="cell-type-2">Cell Type 2:</label>
                            <select id="cell-type-2">
                                {''.join([f'<option value="{ct.replace(" ", "_")}">{ct}</option>' for ct in cell_types])}
                            </select>
                            <div class="quantile-group">
                                <label>Quantile:</label>
                                <label><input type="radio" name="quantile-2" value="high" checked> High</label>
                                <label><input type="radio" name="quantile-2" value="low"> Low</label>
                            </div>
                        </div>
                    </div>
                    
                    <button id="cell-filter-btn" class="filter-button">Filter Files</button>
                </div>
                
                <div class="filtered-files-section">
                    <h3>Filtered Results</h3>
                    <div id="cell-filtered-files" class="file-grid"></div>
                </div>
            </div>
            
            <div class="tab-content" id="gene-expression-content">
                <div class="control-panel">
                    <div class="form-group">
                        <label for="gene-comparison-scope">Comparison Scope:</label>
                        <select id="gene-comparison-scope">
                            <option value="within">Within-sample</option>
                            <option value="cross">Cross-sample</option>
                        </select>
                    </div>
                    
                    <div class="form-group" id="gene-sample-group">
                        <label for="gene-sample-single">Sample:</label>
                        <select id="gene-sample-single">
                            {''.join([f'<option value="{batch}">{batch}</option>' for batch in batch_names])}
                        </select>
                    </div>
                    
                    <div class="form-group hidden" id="gene-sample-pair-group">
                        <label for="gene-sample-a">Sample A:</label>
                        <select id="gene-sample-a">
                            {''.join([f'<option value="{batch}">{batch}</option>' for batch in batch_names])}
                        </select>
                        <label for="gene-sample-b">Sample B:</label>
                        <select id="gene-sample-b">
                            {''.join([f'<option value="{batch}">{batch}</option>' for batch in batch_names])}
                        </select>
                    </div>
                    
                    <button id="gene-filter-btn" class="filter-button">Filter Files</button>
                </div>
                
                <div class="filtered-files-section">
                    <h3>Filtered Results</h3>
                    <div id="gene-filtered-files" class="file-grid"></div>
                </div>
            </div>
            
            <!-- Hidden sections for getAllFiles() -->
            <div style="display: none;">
                {sections_html}
            </div>
        </div>
        """

        # Replace placeholders in correct order
        html_content = html_content.replace("{{EXTRA_CSS}}", extra_css)
        html_content = html_content.replace("{{EXTRA_JS}}", extra_js)
        html_content = html_content.replace("{{CONTENT}}", content)
        
        with open(self.report_dir / filename, 'w') as f:
            f.write(html_content)

    def _get_dge_report_css(self):
        """CSS specific to the DGE report redesign"""
        return """
        .tab-navigation {
            display: flex;
            margin-bottom: 2rem;
            border-bottom: 2px solid #e0e6ed;
        }
        
        .tab-button {
            padding: 1rem 2rem;
            background: #f8f9fa;
            border: none;
            border-radius: 8px 8px 0 0;
            cursor: pointer;
            font-size: 1rem;
            font-weight: 600;
            color: #666;
            transition: all 0.3s ease;
            margin-right: 0.5rem;
        }
        
        .tab-button.active {
            background: #667eea;
            color: white;
        }
        
        .tab-button:hover:not(.active) {
            background: #e9ecef;
        }
        
        .tab-content {
            display: none;
        }
        
        .tab-content.active {
            display: block;
        }
        
        .control-panel {
            background: white;
            border-radius: 15px;
            padding: 1.5rem;
            margin-bottom: 2rem;
            box-shadow: 0 5px 15px rgba(0,0,0,0.08);
        }
        
        .form-group {
            margin-bottom: 1.5rem;
        }
        
        .form-group label {
            display: block;
            margin-bottom: 0.5rem;
            font-weight: 600;
            color: #333;
        }
        
        .form-group select {
            width: 100%;
            padding: 0.75rem;
            border: 2px solid #e0e6ed;
            border-radius: 8px;
            font-size: 1rem;
            background: white;
            transition: border-color 0.3s ease;
        }
        
        .form-group select:focus {
            outline: none;
            border-color: #667eea;
        }
        
        .form-row {
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 1.5rem;
        }
        
        .quantile-group {
            display: flex;
            align-items: center;
            gap: 1rem;
            margin-top: 0.5rem;
        }
        
        .quantile-group label {
            display: flex;
            align-items: center;
            margin: 0;
            font-weight: normal;
            cursor: pointer;
        }
        
        .quantile-group input[type="radio"] {
            margin-right: 0.5rem;
        }
        
        .hidden {
            display: none;
        }
        
        .filter-button {
            padding: 0.75rem 1.5rem;
            background: #28a745;
            color: white;
            border: none;
            border-radius: 8px;
            cursor: pointer;
            font-size: 1rem;
            font-weight: 600;
            transition: background 0.3s ease;
        }
        
        .filter-button:hover {
            background: #218838;
        }
        
        .filtered-files-section {
            background: white;
            border-radius: 15px;
            padding: 1.5rem;
            box-shadow: 0 5px 15px rgba(0,0,0,0.08);
        }
        
        .filtered-files-section h3 {
            margin-bottom: 1rem;
            color: #333;
        }
        
        @media (max-width: 768px) {
            .form-row {
                grid-template-columns: 1fr;
            }
            
            .tab-button {
                padding: 0.75rem 1rem;
                font-size: 0.9rem;
            }
        }
        """

    def _get_dge_report_js(self):
        """JavaScript for the DGE report redesign"""
        return """
        // Tab switching functionality
        function toggleVisibility() {
            var isWithin = document.getElementById('cell-comparison-scope').value === 'within';
            document.getElementById('cell-type-2').parentElement.style.display = isWithin ? 'none' : '';
            document.querySelectorAll('.quantile-group').forEach(function(el) {
            el.style.display = isWithin ? 'none' : 'flex';
            });
        }
        // Initial call to set the visibility based on the default dropdown value on page load
        toggleVisibility();

        document.addEventListener('DOMContentLoaded', function() {
            const tabButtons = document.querySelectorAll('.tab-button');
            const tabContents = document.querySelectorAll('.tab-content');
            
            tabButtons.forEach(button => {
                button.addEventListener('click', () => {
                    const tabId = button.getAttribute('data-tab');
                    
                    // Update active tab button
                    tabButtons.forEach(btn => btn.classList.remove('active'));
                    button.classList.add('active');
                    
                    // Show corresponding tab content
                    tabContents.forEach(content => {
                        content.classList.remove('active');
                        if (content.id === tabId + '-content') {
                            content.classList.add('active');
                        }
                    });
                });
            });
            
            // Handle scope changes for cell type comparison
            document.getElementById('cell-comparison-scope').addEventListener('change', function() {
                const isWithin = this.value === 'within';
                document.getElementById('cell-sample-group').classList.toggle('hidden', !isWithin);
                document.getElementById('cell-sample-pair-group').classList.toggle('hidden', isWithin);
            });
            
            // Handle scope changes for gene expression comparison
            document.getElementById('gene-comparison-scope').addEventListener('change', function() {
                const isWithin = this.value === 'within';
                document.getElementById('gene-sample-group').classList.toggle('hidden', !isWithin);
                document.getElementById('gene-sample-pair-group').classList.toggle('hidden', isWithin);
            });
            
            // Cell type filter button
            document.getElementById('cell-filter-btn').addEventListener('click', filterCellTypeFiles);
            
            // Gene expression filter button
            document.getElementById('gene-filter-btn').addEventListener('click', filterGeneExpressionFiles);
        });
        
        // Filter files for cell type comparison
        function filterCellTypeFiles() {
            const container = document.getElementById('cell-filtered-files');
            container.innerHTML = '';
            const scope = document.getElementById('cell-comparison-scope').value;
            const sampleSingle = document.getElementById('cell-sample-single').value;
            const sampleA = document.getElementById('cell-sample-a').value;
            const sampleB = document.getElementById('cell-sample-b').value;
            const cellType1 = document.getElementById('cell-type-1').value;
            const cellType2 = document.getElementById('cell-type-2').value;
            const quantile1 = document.querySelector('input[name="quantile-1"]:checked').value;
            const quantile2 = document.querySelector('input[name="quantile-2"]:checked').value;
            
            // Get all files from the existing sections
            const allFiles = getAllFiles();
            
            // First filter by path - only include files from cell_type_comparisons or combined_cell_type_comparisons
            const cellTypeFiles = allFiles.filter(file => {
                return file.path.includes('cell_type_comparisons') || file.path.includes('combined_cell_type_comparisons');
            });
            
            let filteredFiles = [];
            
            if (scope === 'within') {
                // Within-sample patterns
                const dgePattern = new RegExp(`dge_${sampleSingle}_${cellType1}_${quantile1}_vs_`);
                const scatter_volcano_Pattern = new RegExp(`volcano_${sampleSingle}_${cellType1}_${quantile1}_vs_`);
                
                filteredFiles = cellTypeFiles.filter(file => {
                    return dgePattern.test(file.path) || 
                           scatter_volcano_Pattern.test(file.path) ||
                           file.path.includes('stacked_barchart') ||
                           file.path.includes('merged_dge_on_');
                });
            } else {
                // Cross-sample patterns
                const dgePattern = new RegExp(`dge_cross_sample_${sampleA}_${quantile1}_vs_${sampleB}_${quantile1}_${cellType1}_`);
                const scatter_volcano_Pattern = new RegExp(`${sampleA}_${quantile1}_.+_${sampleB}_${quantile2}`);
                
                filteredFiles = cellTypeFiles.filter(file => {
                    return dgePattern.test(file.path) || 
                           scatter_volcano_Pattern.test(file.path) ||
                           file.path.includes('stacked_barchart') ||
                           file.path.includes('merged_dge_on_');
                });
            }
            
            displayFilteredFiles(filteredFiles, 'cell-filtered-files');
        }
        
        // Filter files for gene expression comparison
        function filterGeneExpressionFiles() {
            const container = document.getElementById('gene-filtered-files');
            container.innerHTML = '';
            const scope = document.getElementById('gene-comparison-scope').value;
            const sampleSingle = document.getElementById('gene-sample-single').value;
            const sampleA = document.getElementById('gene-sample-a').value;
            const sampleB = document.getElementById('gene-sample-b').value;
            
            // Get all files from the existing sections
            const allFiles = getAllFiles();
            
            // First filter by path - only include files from cluster_comparisons
            const geneExpressionFiles = allFiles.filter(file => {
                return file.path.includes('cluster_comparisons');
            });
            
            let filteredFiles = [];
            
            if (scope === 'within') {
                // Within-sample patterns
                const regionPattern = new RegExp(`(edgeTumour|inTumour|outTumour)_against_(edgeTumour|inTumour|outTumour)_${sampleSingle}`);
                const scatter_volcano_Pattern = new RegExp(`_${sampleSingle}_[^(vs)]+`);
                
                filteredFiles = geneExpressionFiles.filter(file => {
                    return regionPattern.test(file.path) || 
                           scatter_volcano_Pattern.test(file.path) ||
                           file.path.includes('stacked_barchart') ||
                           file.path.includes('merged_dge_on_');
                });
            } else {
                // Cross-sample patterns
                const regionPattern = new RegExp(`(edgeTumour|inTumour|outTumour)_${sampleA}_vs_${sampleB}`);
                const scatter_volcano_Pattern = new RegExp(`_${sampleA}_vs_${sampleB}`);
                
                filteredFiles = geneExpressionFiles.filter(file => {
                    return regionPattern.test(file.path) || 
                           scatter_volcano_Pattern.test(file.path) || 
                           file.path.includes('stacked_barchart') ||
                           file.path.includes('merged_dge_on_');
                });
            }
            
            displayFilteredFiles(filteredFiles, 'gene-filtered-files');
        }
        
        // Get all files from existing sections
        function getAllFiles() {
            const fileCards = document.querySelectorAll('.file-card, .file-link-text');
            const files = [];
            
            fileCards.forEach(card => {
                let fileName = '';
                let filePath = '';
                
                if (card.classList.contains('file-card')) {
                    const fileNameElement = card.querySelector('h4');
                    const filePathElement = card.querySelector('.file-path');
                    if (fileNameElement) fileName = fileNameElement.textContent;
                    if (filePathElement) filePath = filePathElement.textContent;
                } else if (card.classList.contains('file-link-text')) {
                    const linkElement = card.querySelector('a');
                    const filePathElement = card.querySelector('.file-path');
                    if (linkElement) fileName = linkElement.textContent.split(' ')[0];
                    if (filePathElement) filePath = filePathElement.textContent;
                }
                
                files.push({
                    name: fileName,
                    path: filePath,
                    element: card.cloneNode(true)
                });
            });
            
            return files;
        }
        
        // Display filtered files
        function displayFilteredFiles(files, containerId) {
            const container = document.getElementById(containerId);
            
            if (files.length === 0) {
                container.innerHTML = '<div class="no-results"><p>No files match the current filter criteria.</p></div>';
                return;
            }
            
            // Add file elements with proper hyperlinks
            files.forEach(file => {
                // Create a copy of the element to avoid modifying the original
                const fileElement = file.element.cloneNode(true);
                
                // If it's a file card, wrap it in a link
                if (fileElement.classList.contains('file-card')) {
                    const link = document.createElement('a');
                    link.href = `../${file.path}`;
                    link.target = '_blank';
                    link.className = 'file-link';
                    link.appendChild(fileElement);
                    container.appendChild(link);
                } 
                // If it's a file link text, just append it
                else if (fileElement.classList.contains('file-link-text')) {
                    // Update the link href if it exists
                    const link = fileElement.querySelector('a');
                    if (link) {
                        link.href = `../${file.path}`;
                    }
                    container.appendChild(fileElement);
                } 
                // For any other element, just append it
                else {
                    container.appendChild(fileElement);
                }
            });
        }
        """

    def _generate_floating_nav_html(self, section_anchors):
        """Generate HTML for the floating section navigator, with view toggle button next to title"""
        nav_links = ""
        for anchor in section_anchors:
            nav_links += f'<a href="#{anchor["id"]}" class="floating-nav-link">{anchor["title"]}</a>'
        # Button next to title
        nav_header = """
        <div class="floating-nav-header">
            <div class="floating-nav-title">Sections</div>
            <button id="view-toggle-btn" class="floating-view-toggle-btn" onclick="toggleView()" type="button">
                <span class="view-icon">📋</span> List View
            </button>
        </div>
        """
        return f"""
        <nav class="floating-nav" id="floating-nav">
            {nav_header}
            {nav_links}
        </nav>
        """

    def _get_floating_nav_css(self):
        """CSS for the floating section navigator, with nav-header horizontal layout"""
        return """
        .floating-nav {
            position: fixed;
            top: 90px;
            right: 40px;
            z-index: 1000;
            background: rgba(255,255,255,0.97);
            border-radius: 10px;
            box-shadow: 0 4px 18px rgba(0,0,0,0.10);
            padding: 0.7rem 0.7rem 0.7rem 0.7rem;
            min-width: 150px;
            max-width: 220px;
            max-height: 75vh;
            overflow-y: auto;
            font-size: 0.85rem;
            display: flex;
            flex-direction: column;
            gap: 0.3rem;
            transition: box-shadow 0.2s;
            cursor: move;
            user-select: none;
        }
        .floating-nav-header {
            display: flex;
            align-items: center;
            justify-content: space-between;
            gap: 0.5rem;
            margin-bottom: 0.2rem;
        }
        .floating-nav-title {
            font-weight: bold;
            color: #667eea;
            font-size: 0.92rem;
            letter-spacing: 0.01em;
            cursor: move;
        }
        .floating-view-toggle-btn {
            padding: 0.35rem 0.8rem;
            background: #28a745;
            color: white;
            border: none;
            border-radius: 8px;
            cursor: pointer;
            font-size: 0.85rem;
            transition: background 0.3s ease;
            display: flex;
            align-items: center;
            gap: 0.4rem;
        }
        .floating-view-toggle-btn:hover {
            background: #218838;
        }
        .view-icon {
            font-size: 1rem;
        }
        .floating-nav.dragging {
            box-shadow: 0 8px 25px rgba(0,0,0,0.20);
            background: rgba(255,255,255,0.98);
            transition: none;
        }
        .floating-nav-link {
            display: block;
            color: #333;
            text-decoration: none;
            padding: 0.18rem 0.3rem;
            border-radius: 5px;
            transition: background 0.2s, color 0.2s;
            font-size: 0.85rem;
            cursor: pointer;
        }
        .floating-nav-link:hover, .floating-nav-link.active {
            background: #e6eaff;
            color: #3b4890;
        }
        @media (max-width: 900px) {
            .floating-nav {
                display: none;
            }
        }
        """

    def _get_floating_nav_js(self):
        """JS for floating nav: highlight active section, smooth scroll, drag and drop"""
        return """
        // Smooth scroll for floating nav
        document.querySelectorAll('.floating-nav-link').forEach(function(link) {
            link.addEventListener('click', function(e) {
                var targetId = this.getAttribute('href').slice(1);
                var target = document.getElementById(targetId);
                if (target) {
                    e.preventDefault();
                    window.scrollTo({
                        top: target.getBoundingClientRect().top + window.scrollY - 30,
                        behavior: 'smooth'
                    });
                }
            });
        });
        
        // Highlight active section in floating nav
        window.addEventListener('scroll', function() {
            var links = document.querySelectorAll('.floating-nav-link');
            var sections = Array.from(links).map(function(link) {
                var id = link.getAttribute('href').slice(1);
                return document.getElementById(id);
            });
            var scrollPos = window.scrollY + 80;
            var activeIdx = 0;
            for (var i = 0; i < sections.length; i++) {
                if (sections[i] && sections[i].offsetTop <= scrollPos) {
                    activeIdx = i;
                }
            }
            links.forEach(function(link, idx) {
                if (idx === activeIdx) {
                    link.classList.add('active');
                } else {
                    link.classList.remove('active');
                }
            });
        });
        
        // Drag and drop functionality for floating nav
        (function() {
            var floatingNav = document.getElementById('floating-nav');
            if (!floatingNav) return;
            
            var isDragging = false;
            var dragStartX, dragStartY;
            var navStartX, navStartY;
            
            function startDrag(e) {
                // Only start drag if clicking on nav background or title, not on links
                if (e.target.classList.contains('floating-nav-link')) {
                    return;
                }
                
                isDragging = true;
                floatingNav.classList.add('dragging');
                
                dragStartX = e.clientX;
                dragStartY = e.clientY;
                
                var rect = floatingNav.getBoundingClientRect();
                navStartX = rect.left;
                navStartY = rect.top;
                
                e.preventDefault();
            }
            
            function doDrag(e) {
                if (!isDragging) return;
                
                var deltaX = e.clientX - dragStartX;
                var deltaY = e.clientY - dragStartY;
                
                var newX = navStartX + deltaX;
                var newY = navStartY + deltaY;
                
                // Constrain to viewport
                var navRect = floatingNav.getBoundingClientRect();
                var maxX = window.innerWidth - navRect.width;
                var maxY = window.innerHeight - navRect.height;
                
                newX = Math.max(0, Math.min(newX, maxX));
                newY = Math.max(0, Math.min(newY, maxY));
                
                floatingNav.style.left = newX + 'px';
                floatingNav.style.top = newY + 'px';
                floatingNav.style.right = 'auto';
                
                e.preventDefault();
            }
            
            function stopDrag(e) {
                if (!isDragging) return;
                
                isDragging = false;
                floatingNav.classList.remove('dragging');
                e.preventDefault();
            }
            
            // Mouse events
            floatingNav.addEventListener('mousedown', startDrag);
            document.addEventListener('mousemove', doDrag);
            document.addEventListener('mouseup', stopDrag);
            
            // Touch events for mobile (though nav is hidden on mobile)
            floatingNav.addEventListener('touchstart', function(e) {
                if (e.touches.length === 1) {
                    var touch = e.touches[0];
                    startDrag({
                        target: e.target,
                        clientX: touch.clientX,
                        clientY: touch.clientY,
                        preventDefault: function() { e.preventDefault(); }
                    });
                }
            });
            
            document.addEventListener('touchmove', function(e) {
                if (isDragging && e.touches.length === 1) {
                    var touch = e.touches[0];
                    doDrag({
                        clientX: touch.clientX,
                        clientY: touch.clientY,
                        preventDefault: function() { e.preventDefault(); }
                    });
                }
            });
            
            document.addEventListener('touchend', function(e) {
                stopDrag({
                    preventDefault: function() { e.preventDefault(); }
                });
            });
        })();
        """

    def _find_files(self, pattern, subdirs, exclude_subdirs=None, exclude_patterns=None):
        """Find files matching pattern in specified subdirectories, with exclusions"""
        files = []
        exclude_subdirs = exclude_subdirs or []
        exclude_patterns = exclude_patterns or []

        def is_excluded(path):
            for exdir in exclude_subdirs:
                try:
                    if exdir and (self.project_dir / exdir) in path.parents:
                        return True
                except Exception:
                    continue
            for expat in exclude_patterns:
                if expat and path.match(expat):
                    return True
            return False

        if subdirs:
            for subdir in subdirs:
                if isinstance(subdir, str):
                    if not Path(subdir).is_absolute():
                        subdir = self.project_dir / subdir
                    else:
                        subdir = Path(subdir)
                    if subdir.exists():
                        for f in subdir.glob(pattern):
                            if not is_excluded(f):
                                files.append(f)
                else:
                    print(f"Warning: Directory does not exist: {subdir}")

        return sorted(files)
    
    def _generate_file_cards(self, files):
        """Generate HTML cards for file listings"""
        cards = ""
        
        for file_path in files:
            file_type = file_path.suffix.lower()
            file_name = file_path.stem
            file_size = self._get_file_size(file_path)
            
            try:
                relative_path = file_path.relative_to(self.report_dir.parent)
            except ValueError:
                relative_path = file_path
            
            # Direct file access URL
            file_url = f"../{relative_path}"
            
            if file_type == '.pdf':
                # Generate thumbnail if not exists
                thumb_path = self._get_pdf_thumbnail(file_path)
                if thumb_path and thumb_path.exists():
                    thumb_rel = thumb_path.relative_to(self.report_dir.parent)
                    thumb_url = f"../{thumb_rel}"
                    preview_html = f'<img src="{thumb_url}" alt="PDF thumbnail" class="pdf-thumb-img" style="width:100%;height:100%;object-fit:contain;background:#fff;">'
                else:
                    preview_html = """
                        <div class="pdf-thumbnail">
                            <div class="file-icon">📄</div>
                            <div class="file-type-label">PDF</div>
                        </div>
                    """
                cards += f"""
                <a href="{file_url}" target="_blank" class="file-link">
                <div class="file-card" data-file-type="pdf">
                    <div class="file-preview">
                        {preview_html}
                        <div class="file-overlay">
                            <div class="file-type-badge">PDF</div>
                            <div class="file-size-badge">{file_size}</div>
                        </div>
                    </div>
                    <div class="file-info">
                        <h4>{file_name}</h4>
                        <p class="file-path">{relative_path}</p>
                    </div>
                </div>
                </a>
                """
            elif file_type in ['.csv', '.tsv']:
                cards += f"""
                <a href="{file_url}" target="_blank" class="file-link">
                <div class="file-card" data-file-type="data">
                    <div class="file-preview">
                        <div class="csv-preview-static">
                            {self._generate_csv_preview(file_path)}
                        </div>
                        <div class="file-overlay">
                            <div class="file-type-badge">CSV</div>
                            <div class="file-size-badge">{file_size}</div>
                        </div>
                    </div>
                    <div class="file-info">
                        <h4>{file_name}</h4>
                        <p class="file-path">{relative_path}</p>
                    </div>
                </div>
                </a>
                """
            elif file_type == '.html':
                # Use a simple text hyperlink instead of a card
                cards += f"""
                <div class="file-link-text">
                    <a href="{file_url}" target="_blank">{file_name} <span class="file-size">({file_size})</span></a>
                    <span class="file-path">{relative_path}</span>
                </div>
                """
        
        return cards

    def _get_pdf_thumbnail(self, pdf_path):
        """Export first page of PDF as PNG thumbnail, return thumbnail path"""
        try:
            thumb_dir = self.report_dir / "pdf_thumbnails"
            thumb_dir.mkdir(exist_ok=True)
            thumb_path = thumb_dir / (pdf_path.stem + ".png")
            if not thumb_path.exists():
                doc = fitz.open(pdf_path)
                if doc.page_count > 0:
                    page = doc.load_page(0)
                    pix = page.get_pixmap(matrix=fitz.Matrix(0.5, 0.5), alpha=False)
                    pix.save(str(thumb_path))
                doc.close()
            return thumb_path
        except Exception as _:
            # Could not generate thumbnail
            return None

    def _get_file_size(self, file_path):
        """Get human readable file size"""
        try:
            size = file_path.stat().st_size
            for unit in ['B', 'KB', 'MB', 'GB']:
                if size < 1024.0:
                    return f"{size:.1f} {unit}"
                size /= 1024.0
            return f"{size:.1f} TB"
        except Exception as _:
            return "Unknown"

    def _generate_csv_preview(self, file_path):
        """Generate static CSV preview by reading the file"""
        try:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                # Read first few lines to create preview
                reader = csv.reader(f)
                rows = []
                for i, row in enumerate(reader):
                    if i >= 9:  # Read max 9 rows (1 header + 8 data)
                        break
                    rows.append(row)
                
                if not rows:
                    return '<div class="csv-empty">📊 Empty CSV</div>'
                
                # Generate HTML table
                headers = rows[0] if rows else []
                data_rows = rows[1:] if len(rows) > 1 else []
                
                max_cols = min(6, len(headers))
                
                html_content = '<table class="csv-table"><thead><tr>'
                
                # Add headers
                for i in range(max_cols):
                    header = headers[i] if i < len(headers) else ''
                    html_content += f'<th title="{html.escape(header)}">{html.escape(header[:15])}</th>'
                
                if len(headers) > max_cols:
                    html_content += '<th>...</th>'
                
                html_content += '</tr></thead><tbody>'
                
                # Add data rows
                max_rows = min(8, len(data_rows))
                for row_idx in range(max_rows):
                    html_content += '<tr>'
                    row = data_rows[row_idx] if row_idx < len(data_rows) else []
                    
                    for col_idx in range(max_cols):
                        cell = row[col_idx] if col_idx < len(row) else ''
                        cell_display = str(cell)[:15] + ('...' if len(str(cell)) > 15 else '')
                        html_content += f'<td title="{html.escape(str(cell))}">{html.escape(cell_display)}</td>'
                    
                    if len(headers) > max_cols:
                        html_content += '<td>...</td>'
                    
                    html_content += '</tr>'
                
                html_content += '</tbody></table>'
                
                # Add summary info
                total_rows = len(rows)
                total_cols = len(headers)
                if total_rows > 9 or total_cols > max_cols:
                    html_content += f'<div class="csv-more">{total_rows} rows × {total_cols} cols</div>'
                
                return html_content
                
        except Exception as _:
            return '''
                <div class="csv-error">
                    <div class="file-icon">📊</div>
                    <div class="file-type-label">CSV</div>
                    <div class="error-msg">Preview unavailable</div>
                </div>
            '''

    def _get_base_html_template(self, title):
        """Get base HTML template"""
        return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.6;
            color: #333;
            background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
            min-height: 100vh;
        }}
        
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            padding: 20px;
        }}
        
{{{{EXTRA_CSS}}}}
    </style>
</head>
<body>
    <div class="container">
{{{{CONTENT}}}}
    </div>
    
    <script>
{{{{EXTRA_JS}}}}
    </script>
</body>
</html>"""

    def _get_navigation_css(self):
        """Get CSS for main navigation page"""
        return """
        .main-header {
            text-align: center;
            margin-bottom: 3rem;
            padding: 2rem;
            background: white;
            border-radius: 15px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.1);
        }
        
        .main-header h1 {
            font-size: 3rem;
            margin-bottom: 0.5rem;
            background: linear-gradient(45deg, #667eea, #764ba2);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
        }
        
        .subtitle {
            font-size: 1.2rem;
            color: #666;
            margin-bottom: 0.5rem;
        }
        
        .nav-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
            gap: 2rem;
            margin-bottom: 3rem;
        }
        
        .nav-card {
            background: white;
            border-radius: 15px;
            padding: 2rem;
            box-shadow: 0 10px 30px rgba(0,0,0,0.1);
            cursor: pointer;
            transition: transform 0.3s ease, box-shadow 0.3s ease;
            position: relative;
            overflow: hidden;
        }
        
        .nav-card:hover {
            transform: translateY(-10px);
            box-shadow: 0 20px 40px rgba(0,0,0,0.15);
        }
        
        .nav-card.disabled {
            opacity: 0.6;
            cursor: not-allowed;
            background: #f8f9fa;
        }
        
        .nav-card.disabled:hover {
            transform: none;
            box-shadow: 0 10px 30px rgba(0,0,0,0.1);
        }
        
        .nav-icon {
            font-size: 3rem;
            margin-bottom: 1rem;
            text-align: center;
        }
        
        .nav-card h3 {
            font-size: 1.5rem;
            margin-bottom: 1rem;
            color: #333;
        }
        
        .nav-card p {
            color: #666;
            margin-bottom: 1.5rem;
        }
        
        .nav-arrow {
            position: absolute;
            bottom: 1rem;
            right: 1rem;
            font-size: 1.5rem;
            color: #667eea;
            opacity: 0;
            transition: opacity 0.3s ease;
        }
        
        .nav-card:hover .nav-arrow {
            opacity: 1;
        }
        
        .nav-card.disabled .nav-arrow {
            display: none;
        }
        
        .info-section {
            background: white;
            border-radius: 15px;
            padding: 2rem;
            box-shadow: 0 10px 30px rgba(0,0,0,0.1);
        }
        
        .info-section h2, .info-section h3 {
            color: #333;
            margin-bottom: 1rem;
        }
        
        .info-section h3 {
            margin-top: 2rem;
            font-size: 1.3rem;
        }
        
        .info-section ul {
            margin-left: 1.5rem;
            margin-top: 0.5rem;
        }
        
        .info-section li {
            margin-bottom: 0.5rem;
        }
        """
    
    def _get_section_css(self):
        """Get CSS for section pages"""
        return """
        .header {
            display: flex;
            align-items: center;
            gap: 1rem;
            margin-bottom: 2rem;
            padding: 1.5rem;
            background: white;
            border-radius: 15px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.1);
        }
        
        .back-btn {
            padding: 0.5rem 1rem;
            background: #667eea;
            color: white;
            border: none;
            border-radius: 8px;
            cursor: pointer;
            font-size: 0.9rem;
            transition: background 0.3s ease;
        }
        
        .back-btn:hover {
            background: #5a6fd8;
        }
        
        .header h1 {
            font-size: 2.5rem;
            color: #333;
            flex: 1;
        }
        
        .view-controls {
            display: flex;
            gap: 0.5rem;
        }
        
        .view-toggle-btn {
            padding: 0.5rem 1rem;
            background: #28a745;
            color: white;
            border: none;
            border-radius: 8px;
            cursor: pointer;
            font-size: 0.9rem;
            transition: background 0.3s ease;
            display: flex;
            align-items: center;
            gap: 0.5rem;
        }
        
        .view-toggle-btn:hover {
            background: #218838;
        }
        
        .view-icon {
            font-size: 1rem;
        }
        
        .intro {
            background: white;
            border-radius: 15px;
            padding: 2rem;
            margin-bottom: 2rem;
            box-shadow: 0 10px 30px rgba(0,0,0,0.1);
            line-height: 1.8;
        }
        
        .section {
            background: white;
            border-radius: 15px;
            padding: 2rem;
            margin-bottom: 2rem;
            box-shadow: 0 10px 30px rgba(0,0,0,0.1);
        }
        
        .section-header {
            display: flex;
            align-items: center;
            justify-content: space-between;
            margin-bottom: 0.5rem;
            gap: 1rem;
            flex-wrap: wrap;
        }
        
        .section h2 {
            color: #333;
            margin: 0;
            flex-shrink: 0;
        }
        
        .search-container {
            display: flex;
            align-items: center;
            gap: 0.5rem;
            flex-wrap: wrap;
        }
        
        .search-input {
            padding: 0.5rem 1rem;
            border: 2px solid #e0e6ed;
            border-radius: 25px;
            font-size: 0.9rem;
            width: 250px;
            min-width: 200px;
            transition: border-color 0.3s ease, box-shadow 0.3s ease;
            background: #f8f9fa;
        }
        
        .search-input:focus {
            outline: none;
            border-color: #667eea;
            box-shadow: 0 0 0 3px rgba(102, 126, 234, 0.1);
            background: white;
        }
        
        .search-input::placeholder {
            color: #999;
        }
        
        .file-count {
            font-size: 0.85rem;
            color: #666;
            font-weight: 500;
            white-space: nowrap;
        }
        
        .section-description {
            color: #666;
            margin-bottom: 1.5rem;
        }
        
        .no-results {
            grid-column: 1 / -1;
            text-align: center;
            padding: 2rem;
            color: #666;
            font-style: italic;
            background: #f8f9fa;
            border-radius: 10px;
            margin-top: 1rem;
        }
        
        .content.list-view .no-results {
            grid-column: auto;
            margin-bottom: 1rem;
        }
        
        /* Responsive search container */
        @media (max-width: 768px) {
            .section-header {
                flex-direction: column;
                align-items: flex-start;
            }
            
            .search-container {
                width: 100%;
                justify-content: space-between;
            }
            
            .search-input {
                flex-grow: 1;
                min-width: 150px;
            }
        }
        
        /* Card View (Default) */
        .file-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
            gap: 1.5rem;
        }
        
        .file-card {
            background: #f8f9fa;
            border-radius: 10px;
            overflow: hidden;
            box-shadow: 0 5px 15px rgba(0,0,0,0.08);
            transition: transform 0.3s ease, box-shadow 0.3s ease;
        }
        
        .file-card:hover {
            transform: translateY(-5px);
            box-shadow: 0 10px 25px rgba(0,0,0,0.12);
        }
        
        .file-preview {
            height: 180px;
            background: #e9ecef;
            display: flex;
            align-items: center;
            justify-content: center;
            position: relative;
            overflow: hidden;
        }
        
        .file-info {
            padding: 1rem;
        }
        
        .file-info h4 {
            margin-bottom: 0.5rem;
            color: #333;
        }
        
        .file-path {
            font-size: 0.8rem;
            color: #999;
            word-break: break-all;
            margin-bottom: 1rem;
        }
        
        /* List View */
        .content.list-view .file-grid {
            display: block;
        }
        
        .content.list-view .file-card {
            display: flex;
            flex-direction: row;
            margin-bottom: 0.75rem;
            border-radius: 8px;
            overflow: hidden;
            height: auto;
        }
        
        .content.list-view .file-card:hover {
            transform: none;
            box-shadow: 0 3px 10px rgba(0,0,0,0.1);
        }
        
        .content.list-view .file-preview {
            width: 80px;
            height: 60px;
            min-width: 80px;
            border-radius: 0;
        }
        
        .content.list-view .file-preview .file-icon {
            font-size: 1.5rem;
            margin-bottom: 0;
        }
        
        .content.list-view .file-preview .file-type-label {
            font-size: 0.7rem;
        }
        
        .content.list-view .file-preview .csv-table {
            font-size: 0.5rem;
        }
        
        .content.list-view .file-preview .csv-table th,
        .content.list-view .file-preview .csv-table td {
            padding: 1px 2px;
            max-width: 30px;
        }
        
        .content.list-view .file-info {
            padding: 0.75rem;
            flex: 1;
            display: flex;
            flex-direction: column;
            justify-content: center;
        }
        
        .content.list-view .file-info h4 {
            margin-bottom: 0.25rem;
            font-size: 1rem;
        }
        
        .content.list-view .file-path {
            margin-bottom: 0;
            font-size: 0.75rem;
        }
        
        .content.list-view .file-overlay {
            position: static;
            display: flex;
            flex-direction: row;
            gap: 0.25rem;
            margin-left: auto;
            align-items: center;
            padding: 0.5rem;
        }
        
        .content.list-view .file-type-badge,
        .content.list-view .file-size-badge {
            font-size: 0.6rem;
            padding: 1px 4px;
        }
        
        /* Special handling for HTML file links in list view */
        .content.list-view .file-link-text {
            display: flex;
            align-items: center;
            gap: 1rem;
            padding: 0.75rem 1rem;
            background: #f8f9fa;
            border-radius: 8px;
            margin-bottom: 0.5rem;
            box-shadow: 0 2px 5px rgba(0,0,0,0.05);
            transition: box-shadow 0.3s ease;
        }
        
        .content.list-view .file-link-text:hover {
            box-shadow: 0 3px 10px rgba(0,0,0,0.1);
        }
        
        .content.list-view .file-link-text a {
            font-weight: 600;
            color: #007bff;
            text-decoration: none;
            flex: 1;
        }
        
        .content.list-view .file-link-text a:hover {
            text-decoration: underline;
        }
        
        .content.list-view .file-link-text .file-path {
            font-size: 0.75rem;
            color: #999;
            margin: 0;
        }
        
        .file-overlay {
            position: absolute;
            top: 8px;
            right: 8px;
            display: flex;
            flex-direction: column;
            gap: 4px;
        }
        
        .file-type-badge, .file-size-badge {
            background: rgba(0,0,0,0.7);
            color: white;
            padding: 2px 6px;
            border-radius: 4px;
            font-size: 0.7rem;
            font-weight: bold;
        }
        
        .pdf-thumbnail, .html-thumbnail {
            width: 100%;
            height: 100%;
            background: linear-gradient(135deg, #667eea, #764ba2);
            display: flex;
            flex-direction: column;
            align-items: center;
            justify-content: center;
            color: white;
            text-align: center;
        }
        
        .file-icon {
            font-size: 3rem;
            margin-bottom: 0.5rem;
        }
        
        .file-type-label {
            font-size: 1.2rem;
            font-weight: bold;
        }
        
        .csv-preview-static {
            width: 100%;
            height: 100%;
            background: white;
            overflow: hidden;
            position: relative;
        }
        
        .csv-table {
            width: 100%;
            height: 100%;
            font-size: 0.7rem;
            border-collapse: collapse;
        }
        
        .csv-table th,
        .csv-table td {
            border: 1px solid #ddd;
            padding: 2px 4px;
            text-align: left;
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
            max-width: 60px;
        }
        
        .csv-table th {
            background-color: #f8f9fa;
            font-weight: bold;
        }
        
        .csv-more {
            position: absolute;
            bottom: 5px;
            right: 5px;
            background: rgba(0,0,0,0.7);
            color: white;
            padding: 2px 6px;
            border-radius: 3px;
            font-size: 0.6rem;
        }
        
        .csv-empty, .csv-error {
            width: 100%;
            height: 100%;
            display: flex;
            flex-direction: column;
            align-items: center;
            justify-content: center;
            background: linear-gradient(135deg, #667eea, #764ba2);
            color: white;
            text-align: center;
        }
        
        .csv-error .file-icon {
            font-size: 2.5rem;
            margin-bottom: 0.5rem;
        }
        
        .csv-error .file-type-label {
            font-size: 1rem;
            font-weight: bold;
            margin-bottom: 0.3rem;
        }
        
        .error-msg {
            font-size: 0.8rem;
            opacity: 0.8;
        }

        .file-link {
            text-decoration: none;
            color: inherit;
            display: block;
        }
        
        .file-actions {
            display: flex;
            gap: 0.5rem;
            flex-wrap: wrap;
        }
        
        .action-btn {
            padding: 0.4rem 0.8rem;
            border: none;
            border-radius: 5px;
            cursor: pointer;
            font-size: 0.8rem;
            text-decoration: none;
            display: inline-flex;
            align-items: center;
            gap: 0.3rem;
            transition: background 0.3s ease;
            flex: 1;
            justify-content: center;
            min-width: 80px;
        }
        
        .action-btn.primary {
            background: #007bff;
            color: white;
        }
        
        .action-btn.primary:hover {
            background: #0056b3;
        }
        
        .action-btn.secondary {
            background: #6c757d;
            color: white;
        }
        
        .action-btn.secondary:hover {
            background: #545b62;
        }

        .pdf-thumb-img {
            width: 100%;
            height: 100%;
            object-fit: contain;
            background: #fff;
            display: block;
        }
        
        /* File link text styling (for HTML files) */
        .file-link-text {
            margin-bottom: 1rem;
            padding: 1rem;
            background: #f8f9fa;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.05);
        }
        
        .file-link-text a {
            color: #007bff;
            text-decoration: none;
            font-weight: 600;
            font-size: 1rem;
        }
        
        .file-link-text a:hover {
            text-decoration: underline;
        }
        
        .file-link-text .file-size {
            color: #666;
            font-weight: normal;
            font-size: 0.9rem;
        }
        
        .file-link-text .file-path {
            display: block;
            margin-top: 0.5rem;
            font-size: 0.8rem;
            color: #999;
        }
        """

    def _get_view_toggle_js(self):
        """JavaScript for view toggle functionality (selector now inside floating nav)"""
        return """
        // View toggle functionality
        function toggleView() {
            var content = document.querySelector('.content');
            var toggleBtn = document.getElementById('view-toggle-btn');
            var viewIcon = toggleBtn.querySelector('.view-icon');
            
            if (content.classList.contains('list-view')) {
                // Switch to card view
                content.classList.remove('list-view');
                toggleBtn.innerHTML = '<span class="view-icon">📋</span> List View';
                localStorage.setItem('fileViewMode', 'card');
            } else {
                // Switch to list view
                content.classList.add('list-view');
                toggleBtn.innerHTML = '<span class="view-icon">🔲</span> Card View';
                localStorage.setItem('fileViewMode', 'list');
            }
        }
        
        // Initialize view based on saved preference
        document.addEventListener('DOMContentLoaded', function() {
            var savedView = localStorage.getItem('fileViewMode');
            var content = document.querySelector('.content');
            var toggleBtn = document.getElementById('view-toggle-btn');
            
            if (toggleBtn) {
                if (savedView === 'list') {
                    content.classList.add('list-view');
                    toggleBtn.innerHTML = '<span class="view-icon">🔲</span> Card View';
                } else {
                    // Default to card view
                    toggleBtn.innerHTML = '<span class="view-icon">📋</span> List View';
                }
            }
        });
        """

    def _get_search_js(self):
        """JavaScript for search functionality"""
        return r"""
        // Search functionality
        document.addEventListener('DOMContentLoaded', function() {
            var searchInputs = document.querySelectorAll('.search-input');
            
            searchInputs.forEach(function(input) {
                input.addEventListener('input', function() {
                    var sectionId = this.getAttribute('data-section');
                    var section = document.getElementById(sectionId);
                    var searchTerm = this.value.toLowerCase().trim();
                    var fileGrid = section.querySelector('.file-grid');
                    var fileCards = fileGrid.querySelectorAll('.file-card, .file-link-text');
                    var fileCount = section.querySelector('.file-count');
                    var visibleCount = 0;
                    
                    // Split search term into keywords for AND operation
                    var keywords = searchTerm.split(/\s+/).filter(function(keyword) {
                        return keyword.length > 0;
                    });
                    
                    fileCards.forEach(function(card) {
                        var fileName = '';
                        
                        // Get filename from different card types
                        if (card.classList.contains('file-card')) {
                            var fileNameElement = card.querySelector('h4');
                            if (fileNameElement) {
                                fileName = fileNameElement.textContent.toLowerCase();
                            }
                        } else if (card.classList.contains('file-link-text')) {
                            var linkElement = card.querySelector('a');
                            if (linkElement) {
                                fileName = linkElement.textContent.toLowerCase();
                            }
                        }
                        
                        // Also search in file path
                        var filePathElement = card.querySelector('.file-path');
                        var filePath = filePathElement ? filePathElement.textContent.toLowerCase() : '';
                        
                        // Combine filename and path for searching
                        var searchableText = fileName + ' ' + filePath;
                        
                        // Check if ALL keywords are present (AND operation)
                        var matchesAllKeywords = keywords.length === 0 || keywords.every(function(keyword) {
                            return searchableText.includes(keyword);
                        });
                        
                        // Show/hide based on search term
                        if (matchesAllKeywords) {
                            card.style.display = '';
                            visibleCount++;
                        } else {
                            card.style.display = 'none';
                        }
                    });
                    
                    // Update file count
                    fileCount.textContent = '(' + visibleCount + ' of ' + fileCards.length + ' files)';
                    
                    // Show/hide "no results" message
                    var noResults = section.querySelector('.no-results');
                    if (visibleCount === 0 && searchTerm !== '') {
                        if (!noResults) {
                            noResults = document.createElement('div');
                            noResults.className = 'no-results';
                            var keywordText = keywords.length > 1 ? 
                                'keywords: "' + keywords.join('", "') + '"' : 
                                '"' + searchTerm + '"';
                            noResults.innerHTML = '<p>No files found matching ' + keywordText + '</p>';
                            fileGrid.appendChild(noResults);
                        } else {
                            var keywordText = keywords.length > 1 ? 
                                'keywords: "' + keywords.join('", "') + '"' : 
                                '"' + searchTerm + '"';
                            noResults.innerHTML = '<p>No files found matching ' + keywordText + '</p>';
                        }
                        noResults.style.display = 'block';
                    } else if (noResults) {
                        noResults.style.display = 'none';
                    }
                });
                
                // Clear search on Escape key
                input.addEventListener('keydown', function(e) {
                    if (e.key === 'Escape') {
                        this.value = '';
                        this.dispatchEvent(new Event('input'));
                        this.blur();
                    }
                });
            });
        });
        """


if __name__ == "__main__":
    generator = ReportGenerator()
    generator.generate_all_reports()