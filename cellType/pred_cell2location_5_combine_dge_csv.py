# %%
"""
Local CSV report generator - creates Excel reports from CSV files
"""

import pandas as pd
import os
import re
import sys
from pathlib import Path
from pred_cell2location_utils import load_configuration as load_configuration


# def create_excel_report(csv_folder_path, output_excel_path, sheet_name_pattern):
#     """
#     Creates an Excel report from all CSV files in a folder.
#     Each CSV becomes a separate sheet in the Excel file.
#     """
#     print(f"Creating Excel report: {output_excel_path}")
    
#     # Create output directory if it doesn't exist
#     Path(output_excel_path).parent.mkdir(parents=True, exist_ok=True)
    
#     file_count = 0
#     file_list = []
    
#     # First pass: collect all file information
#     for filename in sorted(os.listdir(csv_folder_path)):
#         if filename.endswith(".csv"):
#             file_path = os.path.join(csv_folder_path, filename)
            
#             # Extract sheet name from filename using custom logic
#             sheet_name = extract_sheet_name(filename, sheet_name_pattern)
            
#             # Excel sheet names have a 31 character limit and can't contain certain characters
#             sheet_name = re.sub(r'[\\/*?:\[\]]', '_', sheet_name)[:31]
            
#             try:
#                 # Read CSV file to get dimensions
#                 df = pd.read_csv(file_path)
                
#                 file_count += 1
#                 file_list.append({
#                     "filename": filename, 
#                     "sheet_name": sheet_name, 
#                     "rows": len(df), 
#                     "cols": len(df.columns)
#                 })
                
#             except Exception as e:
#                 print(f"Error processing {filename}: {e}")
#                 error_sheet_name = f"ERROR_{sheet_name}"[:31]
#                 file_list.append({
#                     "filename": filename, 
#                     "sheet_name": error_sheet_name, 
#                     "rows": 0, 
#                     "cols": 0
#                 })
    
#     # Create Excel writer object
#     with pd.ExcelWriter(output_excel_path, engine='openpyxl') as writer:
#         # Create summary sheet first with hyperlinks
#         if file_list:
#             summary_df = pd.DataFrame(file_list)
#             # Add hyperlink formula to sheet_name column
#             summary_df['sheet_name'] = summary_df['sheet_name'].apply(
#                 lambda x: f'=HYPERLINK("#\'{x}\'!A1","{x}")'
#             )
#             summary_df.to_excel(writer, sheet_name='_Summary', index=False)
        
#         # Second pass: create actual data sheets
#         for filename in sorted(os.listdir(csv_folder_path)):
#             if filename.endswith(".csv"):
#                 file_path = os.path.join(csv_folder_path, filename)
                
#                 # Extract sheet name from filename using custom logic
#                 sheet_name = extract_sheet_name(filename, sheet_name_pattern)
                
#                 # Excel sheet names have a 31 character limit and can't contain certain characters
#                 sheet_name = re.sub(r'[\\/*?:\[\]]', '_', sheet_name)[:31]
                
#                 print(f"Processing {filename} -> Sheet: {sheet_name}")
                
#                 try:
#                     # Read CSV file
#                     df = pd.read_csv(file_path)
                    
#                     # Write to Excel sheet
#                     df.to_excel(writer, sheet_name=sheet_name, index=False)
                    
#                 except Exception as e:
#                     print(f"Error processing {filename}: {e}")
#                     # Create an error sheet
#                     error_df = pd.DataFrame({
#                         'Error': [f"Failed to process {filename}"],
#                         'Details': [str(e)]
#                     })
#                     error_sheet_name = f"ERROR_{sheet_name}"[:31]
#                     error_df.to_excel(writer, sheet_name=error_sheet_name, index=False)
    
#     print(f"Excel report created: {output_excel_path}")
#     print(f"Total files processed: {file_count}")
    
#     return file_count, file_list


# def extract_sheet_name(filename, pattern_type):
#     """
#     Extract sheet name from filename based on pattern type.
#     """
#     if "cross_sample" in pattern_type:
#         # For cross-sample files like:
#         # dge_cross_sample_ADI_high_vs_SAL_high_NK_cells_high_abundance_top1000.csv
#         # Extracts to: high_ADI_vs_SAL_NK_cells
#         # This regex uses (high|low) for specific abundance capture and 
#         # backreferences (\2) to ensure these keywords are consistent in the filename.
#         match = re.search(r'dge_cross_sample_(\w+)_(high|low)_vs_(\w+)_\2_(.+?)_\2_abundance_top1000\.csv', filename)
#         if match:
#             # Groups captured:
#             # 1: sample1 (e.g., 'ADI')
#             # 2: abundance (e.g., 'high' or 'low') - this is used for all abundance keywords via \2
#             # 3: sample2 (e.g., 'SAL')
#             # 4: cell_type (e.g., 'NK_cells' or 'Macrophages')
#             sample1, abundance, sample2, cell_type = match.groups()
            
#             # The 'abundance' variable is now directly from the second captured group.
#             # The regex structure ensures it's 'high' or 'low' and consistent.
#             return f"{abundance}_{sample1}_vs_{sample2}_{cell_type}"
#     else:
#         # For within-sample files: dge_ADI_high_vs_low_Monocytes_top1000.csv  
#         # Extract: high_vs_low_ADI_Monocytes (assuming this is the actual pattern)
#         # But based on your example, it seems like: dge_cross_sample_ADI_low_vs_SAL_low_Monocytes_low_abundance_top1000.csv
#         # Let me handle the pattern you provided:
#         match = re.search(r'dge_(.+?)_high_vs_low_top1000\.csv', filename)
#         if match:
#             part = match.group(1)
#             # Split by underscore and reconstruct
#             parts = part.split('_')
#             if len(parts) >= 2:
#                 sample = parts[0]
#                 cell_type = '_'.join(parts[1:])
#                 return f"high_vs_low_{sample}_{cell_type}"
    
#     # Fallback to original filename without extension
#     return os.path.splitext(filename)[0]

def extract_within_sample_info(filename):
    match = re.match(r'(\w+)_(.+)_(high_vs_low)\.csv', filename)
    if match:
        sample, cell_type, comparison = match.groups()
        return cell_type, f"{sample}_{comparison}"
    return None, None

def extract_across_sample_info(filename):
    match = re.match(r'cross_sample_(\w+)_(low|high)_vs_(\w+)_(low|high)_(.+)_(low|high)_abundance\.csv', filename)
    if match:
        s1, a1, s2, a2, cell_type, _ = match.groups()
        return cell_type, f"{s1}_{a1}_vs_{s2}_{a2}"
    return None, None

def merge_csvs_to_one(csv_folder_path, group, pattern_func, merged_rows):
    for filename in sorted(os.listdir(csv_folder_path)):
        if filename.endswith(".csv"):
            cell_type, comparison = pattern_func(filename)
            if cell_type is None or comparison is None:
                print(f"Skipping {filename}: pattern not matched")
                continue
            file_path = os.path.join(csv_folder_path, filename)
            try:
                df = pd.read_csv(file_path)
                df['group'] = group
                df['cell_type'] = cell_type
                df['comparison'] = comparison
                merged_rows.append(df)
            except Exception as e:
                print(f"Error processing {filename}: {e}")

if __name__ == "__main__":
    print("--- Starting Local CSV Merge ---")

    config = load_configuration()
    if not config:
        sys.exit(1)

    # --- choose correct output folders based on combine_cell_types ---
    combine_cell_types = config.get('shared', {}).get('combine_cell_types', False)
    if combine_cell_types:
        csv_folder_path_cross = config['paths']['dge_results_across_sample_combined']
        csv_folder_path_within = config['paths']['dge_results_within_sample_combined']
        output_csv_path = f"{config['paths']['dge_csv_combine_combined']}/merged_dge_on_cellTypes.csv"
    else:
        csv_folder_path_cross = config['paths']['dge_results_across_sample']
        csv_folder_path_within = config['paths']['dge_results_within_sample']
        output_csv_path = f"{config['paths']['dge_csv_combine']}/merged_dge_on_cellTypes.csv"

    merged_rows = []

    try:
        # Across-sample
        merge_csvs_to_one(
            csv_folder_path_cross,
            group="across_sample",
            pattern_func=extract_across_sample_info,
            merged_rows=merged_rows
        )

        # Within-sample
        merge_csvs_to_one(
            csv_folder_path_within,
            group="within_sample",
            pattern_func=extract_within_sample_info,
            merged_rows=merged_rows
        )

        # Concatenate and save
        if merged_rows:
            merged_df = pd.concat(merged_rows, ignore_index=True)
            Path(output_csv_path).parent.mkdir(parents=True, exist_ok=True)
            merged_df.to_csv(output_csv_path, index=False)
            print(f"Merged CSV saved to: {output_csv_path}")
        else:
            print("No data merged.")

    except KeyError as e:
        print(f"Error: Missing configuration key in config/cellType_config.yaml: {e}")
        sys.exit(1)

    print("--- Local CSV merge finished ---")

    # Optionally delete original CSVs
    for folder in [csv_folder_path_within, csv_folder_path_cross]:
        for filename in os.listdir(folder):
            if filename.endswith(".csv"):
                file_path = os.path.join(folder, filename)
                try:
                    os.remove(file_path)
                    print(f"Deleted CSV file: {file_path}")
                except Exception as e:
                    print(f"Error deleting {file_path}: {e}")
# %%
