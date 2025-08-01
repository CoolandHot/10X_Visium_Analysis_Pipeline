import scanpy as sc
import squidpy as sq
import os
import yaml
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# ---- User input ----
with open('config/cellType_config.yaml', 'r') as f:
    config = yaml.safe_load(f)
h5ad_path = config['paths']['spatial_data']
output_dir = config['paths']['results_folder']
os.makedirs(output_dir, exist_ok=True)
sc.settings.figdir = output_dir

# ---- Read AnnData ----
adata = sc.read_h5ad(h5ad_path)

# Use squidpy to plot all batches with correct image placement
if "library_id" in adata.obs.columns:
    # sq.pl.spatial_scatter(
    #     adata,
    #     library_key="library_id",
    #     color="orig.ident",
    #     img=True,
    #     save="_validate_all_batches.pdf"
    # )
    with PdfPages("spatial_scatter_plots.pdf") as pdf:
        for subset_id in adata.obs["library_id"].unique():
            subset = adata[adata.obs["library_id"] == subset_id].copy()
            for k in list(subset.uns['spatial'].keys()):
                if k != subset_id:
                    del subset.uns['spatial'][k]
            sq.pl.spatial_scatter(
                subset,
                library_key="library_id",
                color="orig.ident",
                img=True,
                img_res_key="hires",
                title=subset_id
            )
            fig = plt.gcf()
            pdf.savefig(fig)
            plt.close(fig)
    print("Plotted spatial scatter for all libraries using squidpy.")
else:
    print("No 'library_id' column found in adata.obs. Cannot plot with squidpy.")

print("Validation and plotting complete. Figures saved to:", output_dir)
