"""
WebSlides Presentation Generator - creates HTML presentations from PDF reports
Uses pip-installed webslides package to create interactive presentations from DGE analysis PDFs

When assigning {"results": {"highlights", "title", "body", "footer"}} to a category, The category depth of the assigned content should be less than or equal to 2.
"""

import os
import sys
import yaml
from pathlib import Path
import webslides as ws
import fitz  # PyMuPDF
from pdf2image import convert_from_path
from cellType.pred_cell2location_utils import load_configuration as load_config

# --- Shared modal functionality constants ---
MODAL_CSS_JS = """
/* Modal popup styles */
.image-modal {
    display: none !important;
    position: fixed !important;
    z-index: 9999 !important;
    left: 0 !important;
    top: 0 !important;
    width: 100% !important;
    height: 100% !important;
    background-color: rgba(0,0,0,0.8) !important;
    cursor: pointer !important;
}
.image-modal.show {
    display: flex !important;
    align-items: center !important;
    justify-content: center !important;
}
.image-modal img {
    max-width: 95% !important;
    max-height: 95% !important;
    object-fit: contain !important;
    border-radius: 8px !important;
    box-shadow: 0 4px 20px rgba(0,0,0,0.5) !important;
    cursor: default !important;
}
.image-modal .close-text {
    position: absolute !important;
    top: 20px !important;
    right: 30px !important;
    color: white !important;
    font-size: 24px !important;
    font-weight: bold !important;
    cursor: pointer !important;
}
</style>
<script>
// Image modal functionality
document.addEventListener('DOMContentLoaded', function() {
    // Create modal element
    const modal = document.createElement('div');
    modal.className = 'image-modal';
    modal.innerHTML = '<span class="close-text">&times;</span><img>';
    document.body.appendChild(modal);
    
    const modalImg = modal.querySelector('img');
    const closeBtn = modal.querySelector('.close-text');
    
    // Add click listeners to all images
    document.addEventListener('click', function(e) {
        if (e.target.tagName === 'IMG' && (e.target.closest('.image-container, [style*="grid"]') || e.target.closest('div[style*="display: flex"]'))) {
            modal.classList.add('show');
            modalImg.src = e.target.src;
            modalImg.alt = e.target.alt;
            document.body.style.overflow = 'hidden';
        }
    });
    
    // Close modal functions
    function closeModal() {
        modal.classList.remove('show');
        document.body.style.overflow = '';
    }
    
    // Close on background click
    modal.addEventListener('click', function(e) {
        if (e.target === modal) {
            closeModal();
        }
    });
    
    // Close on X button click
    closeBtn.addEventListener('click', closeModal);
    
    // Close on Escape key
    document.addEventListener('keydown', function(e) {
        if (e.key === 'Escape' && modal.classList.contains('show')) {
            closeModal();
        }
    });
});
</script>
<style>
"""

DGE_SPECIFIC_CSS = """
.image-grid {
    display: grid !important;
    grid-template-columns: repeat(2, 1fr) !important;
    grid-template-rows: repeat(2, 1fr) !important;
    gap: 20px !important;
    height: 80vh !important;
}
.image-container {
    position: relative !important;
    border: 1px solid #ddd !important;
    border-radius: 5px !important;
    overflow: hidden !important;
    cursor: pointer !important;
    transition: transform 0.2s !important;
}
.image-container:hover {
    transform: scale(1.02) !important;
    box-shadow: 0 4px 12px rgba(0,0,0,0.15) !important;
}
.image-container h4 {
    position: absolute !important;
    top: 10px !important;
    left: 10px !important;
    background: rgba(255,255,255,0.9) !important;
    padding: 5px !important;
    border-radius: 3px !important;
    margin: 0 !important;
    z-index: 1 !important;
    font-size: 14px !important;
    pointer-events: none !important;
}
.image-container img {
    cursor: pointer !important;
}
"""

def get_combined_css(extra_css=""):
    """Get combined CSS with modal functionality and optional extra styles."""
    return MODAL_CSS_JS + extra_css

def create_clickable_image_container(rel_path, title, container_style="", img_style="", title_style=""):
    """Create a standardized clickable image container with hover effects."""
    default_container_style = "cursor: pointer; transition: transform 0.2s;"
    default_img_style = "cursor: pointer;"
    default_title_style = "pointer-events: none;"
    
    return f'''
    <div style="{container_style} {default_container_style}"
         onmouseover="this.style.transform='scale(1.02)';this.style.boxShadow='0 4px 12px rgba(0,0,0,0.15)'"
         onmouseout="this.style.transform='scale(1)';this.style.boxShadow='none'">
        <h4 style="{title_style} {default_title_style}">{title}</h4>
        <img src="{rel_path}" style="{img_style} {default_img_style}" alt="{title}">
    </div>'''

def create_placeholder_container(message, container_style=""):
    """Create a placeholder container for missing images."""
    return f'''
    <div style="background:#f0f0f0;border:2px dashed #ccc;
                display:flex;align-items:center;justify-content:center;
                border-radius:5px;{container_style}">
        <span style="color:#666;text-align:center;">{message}</span>
    </div>'''

def normalize_name(name):
    return name.lower().replace(" ", "_")

def extract_subplots_by_title(
    pdf_path,
    output_folder,
    lookup_keyword,
    vertical_offset_percentage=0.08,
    horizontal_blank_percentage=0.18,
):
    Path(output_folder).mkdir(parents=True, exist_ok=True)
    
    existing_files = [
        f for f in os.listdir(output_folder)
        if f.endswith('.png') and normalize_name(os.path.splitext(f)[0]) == normalize_name(lookup_keyword)
    ]
    if existing_files:
        print(f"    ✓ Found existing subplot(s) for '{lookup_keyword}': {existing_files}")
        return
    
    # Open the PDF file with PyMuPDF
    doc = fitz.open(pdf_path)

    # Convert PDF to a list of images
    images = convert_from_path(pdf_path)

    # Iterate over each page
    for page_num in range(len(doc)):
        page = doc.load_page(page_num)  # Load the current page
        blocks = page.get_text("dict")["blocks"]  # Extract text blocks

        # Find titles in the text blocks
        titles = []
        for block in blocks:
            if block["type"] == 0:  # Text block
                for line in block["lines"]:
                    for span in line["spans"]:
                        text = span["text"].strip()
                        if (
                            text and text[0].isupper()
                        ):  # Simple heuristic to detect titles
                            titles.append((text, span["bbox"]))

        # Crop and save subplots based on titles
        if page_num < len(images):
            image = images[page_num]
            width, height = image.size

            # Scale PDF bbox coords to image pixel coords
            page_width, page_height = page.rect.width, page.rect.height
            factor_x = width / page_width
            factor_y = height / page_height
            scaled_titles = []
            for title, bbox in titles:
                scaled_bbox = (
                    bbox[0] * factor_x,
                    bbox[1] * factor_y,
                    bbox[2] * factor_x,
                    bbox[3] * factor_y,
                )
                scaled_titles.append((title, scaled_bbox))

            # Define the subplot crop coordinates and print them
            # dynamically choose 1×1, 1×2 or 2×2 grid based on detected titles
            num_plots = len(titles)-1
            if num_plots == 1:
                rows, cols = 1, 1
            elif num_plots == 2:
                rows, cols = 1, 2
            else:
                rows, cols = 2, 2

            vertical_offset = int(vertical_offset_percentage * height)
            horizontal_blank = int(horizontal_blank_percentage * width)
            blank_x = horizontal_blank if cols > 1 else 0
            # height below the top offset divided by rows
            plot_height = (height - vertical_offset) // rows
            # width minus blank_x divided by cols
            plot_width = (width - blank_x) // cols

            # assign human-readable labels
            if rows == 1 and cols == 1:
                labels = ["Single"]
            elif rows == 1 and cols == 2:
                labels = ["Left", "Right"]
            else:
                labels = ["Top-Left", "Top-Right", "Bottom-Left", "Bottom-Right"]

            coords_labels = []
            idx = 0
            for r in range(rows):
                for c in range(cols):
                    x1 = c * (plot_width + blank_x)
                    x2 = x1 + plot_width if c < cols - 1 else width
                    y1 = vertical_offset + r * plot_height
                    y2 = y1 + plot_height if r < rows - 1 else height
                    coords_labels.append(((x1, y1, x2, y2), labels[idx]))
                    idx += 1

            subplots = [(image.crop(coord), lbl) for coord, lbl in coords_labels]

            for i, (subplot, position) in enumerate(subplots):
                for title, bbox in scaled_titles:
                    # Check if the title's bounding box is within the subplot's area
                    if normalize_name(lookup_keyword) == normalize_name(title):
                        if (
                            bbox[0] >= (i % cols) * (plot_width + horizontal_blank)
                            and bbox[2] <= (i % cols) * (plot_width + horizontal_blank) + plot_width
                            and bbox[1] >= vertical_offset + (i // cols) * plot_height
                            and bbox[3] <= vertical_offset + (i // cols + 1) * plot_height
                        ):
                            filename = f"{title}.png"
                            save_path = os.path.join(output_folder, filename)
                            subplot.save(save_path)
                            print(f"    ✓ Saved subplot: {filename}")



def load_config(config_file="config/cellType_config.yaml"):
    """Loads configuration from a YAML file."""
    try:
        with open(config_file, "r") as f:
            config_data = yaml.safe_load(f)

        return config_data

    except FileNotFoundError:
        print(f"Error: Configuration file '{config_file}' not found.")
        return None
    except yaml.YAMLError as e:
        print(f"Error parsing YAML file: {e}")
        return None


def extract_pdf_first_page(pdf_path, output_dir, base_name, config):
    """Extract the first page of a PDF as an image."""
    try:
        # Import here to provide better error message if not available
        import pdf2image

        image_path = os.path.join(
            output_dir, f"{base_name}.{config['image_format']}"
        )
        # Check if the image already exists
        if os.path.exists(image_path):
            print(f"Image already exists: {image_path}")
            return image_path
        
        # Convert first page to image
        images = pdf2image.convert_from_path(
            pdf_path, first_page=1, last_page=1, dpi=config["image_dpi"]
        )

        if images:
            images[0].save(image_path, config["image_format"].upper())
            print(f"Extracted first page: {image_path}")
            return image_path
        else:
            print(f"No images extracted from {pdf_path}")
            return None

    except ImportError:
        print("Error: pdf2image library not found. Install with: pip install pdf2image")
        return None
    except Exception as e:
        print(f"Error extracting PDF {pdf_path}: {e}")
        return None


def process_pdf_to_image(pdf_path, images_dir, base_name, config):
    """Process a single PDF to extract its first page as an image."""
    image_format = config['shared']['image_settings']['format']
    image_dpi = config['shared']['image_settings']['dpi']
    image_path = os.path.join(images_dir, f"{base_name}.{image_format}")

    if os.path.exists(pdf_path):
        print("    ✓ PDF found")
        # Only extract if image doesn't exist
        if not os.path.exists(image_path):
            image_path = extract_pdf_first_page(pdf_path, images_dir, base_name, 
                                              {"image_dpi": image_dpi, "image_format": image_format})
        else:
            print(f"    ✓ Image already exists: {image_path}")
        return image_path
    else:
        print(f"    ✗ PDF not found: {pdf_path}")
        return None

def create_slide_content(images_data, config):
    """Create HTML content for a 2x2 image grid."""
    rows = config['local_report']['dge_presentations']['layout_rows']
    cols = config['local_report']['dge_presentations']['layout_cols']

    html_content = f'''
    <div style="display: grid; grid-template-columns: repeat({cols}, 1fr); grid-template-rows: repeat({rows}, 1fr); gap: 20px; height: 80vh;">
    '''

    for img_data in images_data:
        if img_data and img_data["path"] and os.path.exists(img_data["path"]):
            rel_path = os.path.relpath(img_data["path"], config['paths']['dge_report'])
            html_content += create_clickable_image_container(
                rel_path=rel_path,
                title=img_data["cell_type"],
                container_style="position: relative; border: 1px solid #ddd; border-radius: 5px; overflow: hidden;",
                img_style="width: 100%; height: 100%; object-fit: contain;",
                title_style="position: absolute; top: 10px; left: 10px; padding: 5px; border-radius: 3px; margin: 0; z-index: 1; font-size: 14px; background: rgba(255,255,255,0.9);"
            )
        else:
            html_content += create_placeholder_container(
                f"Image not available<br>{img_data['cell_type'] if img_data else ''}"
            )

    html_content += "</div>"
    return html_content

def create_slide_for_comparison(
    pdf_dir, filename_pattern, comparison_info, config, cell_type_group
):
    """Create a slide for a comparison with 2x2 image grid for a specific cell type group."""
    images_data = []

    print(
        f"\nProcessing {comparison_info['description']} - Group: {cell_type_group['group_name']}"
    )

    # Process each cell type in the group
    for cell_type in cell_type_group["cell_types"]:
        # Convert cell type name to underscore format for file operations
        cell_type_file = cell_type.replace(" ", "_")
        
        # Generate PDF filename using the pattern with underscore format
        pdf_filename = filename_pattern.format(cell_type=cell_type_file, **comparison_info)
        pdf_path = os.path.join(pdf_dir, pdf_filename)

        print(f"  Looking for: {pdf_filename}")

        # Generate base name for image using underscore format
        base_name = comparison_info["base_name"].format(cell_type=cell_type_file)

        # Process PDF to image
        image_path = process_pdf_to_image(
            pdf_path, os.path.join(config['paths']['dge_report'], "images"), base_name, config
        )

        images_data.append(
            {
                "path": image_path,
                "cell_type": cell_type,  # Keep original format with spaces for display
                "sheet_name": comparison_info["sheet_name"].format(cell_type=cell_type_file),
            }
        )

    # Create slide content with 2x2 image grid
    slide_body = create_slide_content(images_data, config)

    return {
        "title": f"{comparison_info['slide_title']} - {cell_type_group['group_name']}",
        "highlights": comparison_info["highlights"]
        + [f"- Cell type group: {cell_type_group['group_name']}"],
        "body": slide_body,
        "footer": [f"- Generated from PDF reports in {pdf_dir}"],
    }

# --- new functions for cell‐population presentation ---
def create_cell_population_slide(batch, config):
    """Create one slide showing combined + per‐cell plots for a batch."""
    settings = config["local_report"]["cell_population"]
    out_dir = config['paths']['cell_population_reports']
    out_images = os.path.join(out_dir, "images")
    Path(out_images).mkdir(parents=True, exist_ok=True)

    images_data = []
    # 1) combined plot first page
    sel_pdf = settings["selected_pdf_pattern"].format(
        output_dir=config['paths']['cell_abundance_heatmap_path'],
        batch=batch
    )
    img1 = extract_pdf_first_page(
        sel_pdf,
        out_images,
        f"{batch}_selected_cell_types",
        {
            "image_dpi": config['shared']['image_settings']['dpi'], 
            "image_format": config['shared']['image_settings']['format']
        },
    )
    images_data.append({"path": img1, "cell_type": "Combined Plot"})

    # 2) per‐cell subplots
    all_pdf = settings["all_cell_abundance_pdf_pattern"].format(
        output_dir=config['paths']['cell_abundance_heatmap_path'],
        batch=batch
    )
    sub_fld = os.path.join(out_images, f"{batch}_subplots")
    Path(sub_fld).mkdir(parents=True, exist_ok=True)
    for ct in settings["cell_types_for_population"]:
        extract_subplots_by_title(all_pdf, sub_fld, lookup_keyword=ct)
        # pick first match
        matches = [
            f for f in os.listdir(sub_fld) if normalize_name(ct) == normalize_name(os.path.splitext(f)[0])
        ]
        img_path = os.path.join(sub_fld, matches[0]) if matches else None
        images_data.append({"path": img_path, "cell_type": ct})

    # build HTML grid
    rows, cols = config['local_report']['cell_population']['layout_rows'], config['local_report']['cell_population']['layout_cols']
    html = f'<div style="display:grid;grid-template-columns:repeat({cols},1fr);grid-template-rows:repeat({rows},1fr);gap:5px;height:85vh;width:100%;">'
    
    for img in images_data:
        if img["path"] and os.path.exists(img["path"]):
            rel = os.path.relpath(img["path"], out_dir)
            html += create_clickable_image_container(
                rel_path=rel,
                title=img["cell_type"],
                container_style="position:relative;border:1px solid #ddd;border-radius:5px;overflow:hidden;display:flex;flex-direction:column;min-height:0;",
                img_style="width:100%;height:100%;object-fit:contain;flex:1;",
                title_style="position:absolute;top:5px;left:5px;background:rgba(255,255,255,0.9);padding:3px 6px;margin:0;font-size:12px;z-index:10;border-radius:3px;"
            )
        else:
            html += create_placeholder_container(
                f"Image not available<br>{img['cell_type']}",
                "min-height:0;"
            )
    html += "</div>"

    return {
        "title": f"Cell Population - {batch}",
        "highlights": [],
        "body": html,
        "footer": [f"- Combined plot: {sel_pdf}", f"- Abundance details: {all_pdf}"],
    }

def generate_cell_population_presentations(config):
    """Generate a WebSlides presentation for cell population plots."""
    out_dir = config['paths']['cell_population_reports']
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    # wrap slides under one top‐level category
    content = {"Cell Population": {}}
    for batch in config['shared']['within_sample_groups']:
        content["Cell Population"][batch] = {
            "results": create_cell_population_slide(batch, config)
        }

    pres_path = os.path.join(out_dir, "cell_population_presentation.html")
    title_page = {
        "title": "Cell Population Analysis",
        "summary": {
            "Context": (
                "Cell2Location estimates cell type populations at the spot level, allowing each spot to contain multiple cell types. "
                "The presented heatmaps visualize the distribution of cell type fractions across entire tissue sections, with color scales standardized across all samples. "
                "These visualizations help reveal spatial migration patterns and abundance differences for each cell type."
            ),
            "Author": "Huan",
        }
    }
    ws.create(
        content=content,
        title_page=title_page,
        fname=pres_path,
        open_in_browser=False,
        show_index_page=True,
        show_topcat=True,
        show_subcat=True,
        custom_css=get_combined_css(),
    )
    print(f"Cell‐population presentation available at: {pres_path}")
    return True

# --- new functions for NMF analysis presentation ---
def create_nmf_heatmap_slide(n_factor, config):
    """Create one slide showing only the cell type fractions heatmap."""
    nmf_settings = config["paths"]["nmf_analysis_path"]
    out_dir = config['paths']['nmf_reports']
    out_images = os.path.join(out_dir, "images")
    Path(out_images).mkdir(parents=True, exist_ok=True)

    # Cell type fractions heatmap
    heatmap_pdf = os.path.join(
        nmf_settings,
        "cell_type_fractions_heatmap",
        f"n_fact{n_factor}.pdf"
    )
    heatmap_img = extract_pdf_first_page(
        heatmap_pdf,
        out_images,
        f"n_fact{n_factor}_cell_type_fractions",
        {
            "image_dpi": config['shared']['image_settings']['dpi'], 
            "image_format": config['shared']['image_settings']['format']
        },
    )

    # Build HTML for full-size heatmap
    if heatmap_img and os.path.exists(heatmap_img):
        rel_path = os.path.relpath(heatmap_img, out_dir)
        html = f'''
        <div style="display: flex; justify-content: center; align-items: center; height: 85vh;">
            {create_clickable_image_container(
                rel_path=rel_path,
                title=f"Cell Type Fractions - Factor {n_factor}",
                container_style="width: 90%; height: 90%; border: 1px solid #ddd; border-radius: 5px; overflow: hidden; position: relative;",
                img_style="width: 100%; height: 100%; object-fit: contain;",
                title_style="position: absolute; top: 15px; left: 15px; background: rgba(255,255,255,0.9); padding: 8px 12px; margin: 0; font-size: 16px; z-index: 10; border-radius: 5px;"
            )}
        </div>'''
    else:
        html = create_placeholder_container(
            f"Cell Type Fractions heatmap not available<br>Factor {n_factor}",
            "padding: 40px; font-size: 18px; margin: auto; width: fit-content; height: fit-content;"
        )

    return {
        "title": f"NMF Factor {n_factor} - Cell Type Fractions Dotplot",
        "highlights": [
            f"- NMF Factor: {n_factor}",
            "- Cell type composition across factors",
            "- Dotplot showing the weights of each cell type contribution",
        ],
        "body": html,
        "footer": [f"- Source: {heatmap_pdf}"],
    }

def create_nmf_spatial_slide(n_factor, sample, config):
    """Create one slide showing spatial density plot for a specific sample."""
    nmf_settings = config["paths"]["nmf_analysis_path"]
    out_dir = config['paths']['nmf_reports']
    out_images = os.path.join(out_dir, "images")
    Path(out_images).mkdir(parents=True, exist_ok=True)

    # Spatial density plot for the sample
    spatial_pdf = os.path.join(
        nmf_settings, 
        "spatial", 
        f"showcell_density_mean_n_fact{n_factor}_s{sample}_p99.2.pdf"
    )
    spatial_img = extract_pdf_first_page(
        spatial_pdf,
        out_images,
        f"{sample}_n_fact{n_factor}_spatial_density",
        {
            "image_dpi": config['shared']['image_settings']['dpi'], 
            "image_format": config['shared']['image_settings']['format']
        },
    )

    # Build HTML for full-size spatial plot
    if spatial_img and os.path.exists(spatial_img):
        rel_path = os.path.relpath(spatial_img, out_dir)
        html = f'''
        <div style="display: flex; justify-content: center; align-items: center; height: 85vh;">
            {create_clickable_image_container(
                rel_path=rel_path,
                title=f"Spatial Density - {sample} (Factor {n_factor})",
                container_style="width: 90%; height: 90%; border: 1px solid #ddd; border-radius: 5px; overflow: hidden; position: relative;",
                img_style="width: 100%; height: 100%; object-fit: contain;",
                title_style="position: absolute; top: 15px; left: 15px; background: rgba(255,255,255,0.9); padding: 8px 12px; margin: 0; font-size: 16px; z-index: 10; border-radius: 5px;"
            )}
        </div>'''
    else:
        html = create_placeholder_container(
            f"Spatial density plot not available<br>{sample} - Factor {n_factor}",
            "padding: 40px; font-size: 18px; margin: auto; width: fit-content; height: fit-content;"
        )

    return {
        "title": f"NMF Factor {n_factor} - Spatial Density ({sample})",
        "highlights": [
            f"- NMF Factor: {n_factor}",
            f"- Sample: {sample}",
            "- Spatial distribution of cellular patterns",
            "- Density visualization at p99.2 threshold"
        ],
        "body": html,
        "footer": [f"- Source: {spatial_pdf}"],
    }

def generate_nmf_presentations(config):
    """Generate a WebSlides presentation for NMF analysis."""
    out_dir = config['paths']['nmf_reports']
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    print("\n=== NMF Analysis Presentation ===")
    
    # Get factor range from config
    n_fact_start, n_fact_end, n_fact_step = config['nmf_analysis']['n_fact_range']
    factors = list(range(n_fact_start, n_fact_end, n_fact_step))
    
    print(f"Processing factors: {factors}")
    print(f"Processing samples: {config['shared']['within_sample_groups']}")

    # Create content structure with separate categories for heatmaps and spatial plots
    content = {}
    
    # Process heatmap slides (one per factor)
    for n_factor in factors:
        factor_key = f"Factor_{n_factor}"
        print(f"  Processing Heatmap for Factor {n_factor}")
        
        content[factor_key] = {}

        heatmap_slide = create_nmf_heatmap_slide(n_factor, config)
        content[factor_key]["NMF_Heatmaps"] = {"results": heatmap_slide}
        
        for sample in config['shared']['within_sample_groups']:
            sample_key = f"{sample}"
            print(f"  Processing Spatial plot for Factor {n_factor}, Sample {sample}")
            
            spatial_slide = create_nmf_spatial_slide(n_factor, sample, config)
            content[factor_key][sample_key] = {"results": spatial_slide}

    # Generate presentation
    pres_path = os.path.join(out_dir, "nmf_analysis_presentation.html")
    title_page = {
        "title": "NMF Analysis - Cellular Colocalization Patterns",
        "summary": {
            "Context": (
                "This presentation summarizes the results of Non-negative Matrix Factorization (NMF) analysis applied to spatial transcriptomics data. "
                "NMF is a powerful dimensionality reduction technique that decomposes complex gene expression matrices into a set of interpretable factors, "
                "each representing a distinct pattern of cell type colocalization or spatial organization within the tissue."
            ),
            "Purpose": (
                "By analyzing these factors, we can uncover hidden cellular niches, co-localization of cell types, and spatially restricted biological processes. "
                "NMF helps to identify regions where specific cell types tend to co-occur, which may correspond to functional microenvironments or disease-relevant structures."
            ),
            "Factors": f"Analyzed factors: {factors}",
            "Samples": f"Samples: {', '.join(config['shared']['within_sample_groups'])}",
            "Structure": (
                "For each factor, we first display a dotplot showing the weights of each cell type, highlighting which cell types are most associated with that factor. "
                "This is followed by a spatial density heatmap of the factor's activity (mean_nUMI_factors), visualizing how the cellular pattern is distributed across the tissue section. "
                "Together, these plots allow us to interpret both the cellular composition and the spatial localization of each discovered pattern."
            ),
            "Note": (
                "The full-resolution heatmaps (cell type fractions) and spatial density plots are available in the analysis output directory: "
                f"{config['paths']['nmf_analysis_path']}"
            ),
            "Author": "Huan",
        }
    }
    
    ws.create(
        content=content,
        title_page=title_page,
        fname=pres_path,
        open_in_browser=False,
        show_index_page=True,
        show_topcat=True,
        show_subcat=True,
        custom_css=get_combined_css(".nmf-full-image { width: 90% !important; height: 90% !important; }"),
    )
    
    # Calculate total slides
    heatmap_slides = len(factors)
    spatial_slides = len(factors) * len(config['shared']['within_sample_groups'])
    total_slides = heatmap_slides + spatial_slides
    
    print("NMF presentation created:")
    print(f"  - Heatmap slides: {heatmap_slides}")
    print(f"  - Spatial slides: {spatial_slides}")
    print(f"  - Total slides: {total_slides}")
    print(f"NMF presentation available at: {pres_path}")
    return True

# Post-process: replace all absolute paths with relative ones in generated HTML files
def replace_absolute_paths_in_html(html_path, old_prefix, new_prefix):
    if not os.path.exists(html_path):
        print(f"File not found: {html_path}")
        return
    with open(html_path, "r", encoding="utf-8") as f:
        content = f.read()
    new_content = content.replace(old_prefix, new_prefix)
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(new_content)
    print(f"Replaced '{old_prefix}' with '{new_prefix}' in {html_path}")


if __name__ == "__main__":
    print("--- Starting WebSlides Presentation Generation ---")

    # Load configuration
    config = load_config()
    if not config:
        sys.exit(1)

    # Generate all presentations
    generate_cell_population_presentations(config)
    generate_nmf_presentations(config)

    # Replace local absolute paths in generated HTML files
    base_path = str(Path(__file__).parent.parent.parent.absolute())
    presentation_files = [
        os.path.join(config['paths']['cell_population_reports'], "cell_population_presentation.html"),
        os.path.join(config['paths']['nmf_reports'], "nmf_analysis_presentation.html"),
    ]
    
    for html_file in presentation_files:
        replace_absolute_paths_in_html(html_file, base_path, ".")

    print("--- All presentations generated successfully ---")

