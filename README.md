# 10X Visium Analysis Pipeline

This repository contains a comprehensive analysis pipeline for 10X Genomics Visium and VisiumHD spatial transcriptomics data. The pipeline includes preprocessing, clustering, cell type annotation, differential gene expression analysis, and visualization components.

## Overview

The pipeline is designed to process spatial transcriptomics data from 10X Genomics Visium and VisiumHD platforms. It provides tools for:

- Data preprocessing and quality control
- Spatial clustering and cell type annotation
- Differential gene expression (DGE) analysis
- Transcription factor activity analysis
- Visualization of spatial gene expression patterns
- Generation of interactive web reports

## Key Components

### Python Scripts
- `1_VisiumHD_optimized.py`: Optimized processing pipeline for VisiumHD data
- `3_cellbrowser_export_web.py`: Export data for web-based cell browser visualization
- `4_generate_reports.py`: Generate comprehensive analysis reports

### R Scripts
- `1_VisiumHD_read_merge_cluster.r`: Data loading, merging, and initial clustering
- `2_visual_clusters.r`: Visualization of clustering results
- `3_wrap_analysis_to_Loupe.r`: Prepare data for Loupe browser
- Cell type annotation workflows:
  - CARD-based annotation
  - Cell2location-based deconvolution
  - SpaceXR-based annotation

### Analysis Modules

#### Spatial Clustering
- BayesSpace
- BANKSY
- DPNBVI

#### Cell Type Prediction
- Cell2location for cell type deconvolution
- CARD for cell annotation
- SpaceXR for spatial context-aware annotation

#### Differential Gene Expression (DGE)
- Cluster-based DGE analysis
- Cell type-based DGE analysis
- Pathway analysis of DGE results

#### Transcription Factor Analysis (TFA)
- Transcription factor activity prediction
- Important gene identification
- TF activity visualization

## Requirements

- Python 3.8+
- R 4.0+
- Scanpy
- Seurat
- Squidpy
- SpatialData
- Cell2location
- SpaceXR
- CARD

## Usage

1. Configure input data paths in the scripts
2. Run the preprocessing pipeline:
   ```
   Rscript 1_VisiumHD_read_merge_cluster.r
   ```
3. Perform clustering and annotation:
   ```
   Rscript 2_visual_clusters.r
   ```
4. Run downstream analysis:
   ```
   python 1_VisiumHD_optimized.py
   ```

## Directory Structure

```
.
├── DGE/                 # Differential Gene Expression analysis
├── TFA/                 # Transcription Factor Analysis
├── cellType/            # Cell type prediction and annotation
├── cluster_labelling/   # Manual cluster labeling
├── config/              # Configuration files
├── spatial_clustering/  # Spatial clustering algorithms
├── rds_data/            # Processed R data objects
├── output/              # Analysis outputs
└── util_headers.r       # Utility functions
```

## Contributing

This is an active research project. Contributions are welcome through issues and pull requests.

## License

This project is licensed under the MIT License - see the LICENSE file for details.