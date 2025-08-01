#!/bin/bash
cd cellType

# Parse output_dir and cell_type_cell2loc from config/batch_config.yaml, trim whitespace
output_dir=$(grep '^output_dir:' ../config/batch_config.yaml | sed 's/output_dir:[[:space:]]*["'\'']\?\([^"'\'']*\)["'\'']\?/\1/' | sed 's|/$||' | xargs)
cell_type_cell2loc=$(grep 'cell_type_cell2loc:' ../config/batch_config.yaml | sed 's/cell_type_cell2loc:[[:space:]]*["'\'']\?\([^"'\'']*\)["'\'']\?/\1/' | sed 's|/$||' | xargs)

# Parse spatial_output from config/cellType_config.yaml, trim whitespace
spatial_output=$(grep 'spatial_output:' ../config/cellType_config.yaml | sed 's/spatial_output:[[:space:]]*["'\'']\?\([^"'\'']*\)["'\'']\?/\1/' | sed 's|/$||' | xargs)

# Compose destination and source directories
dest_dir="../${output_dir}/${cell_type_cell2loc}"
src_dir="${spatial_output}"

echo "Destination directory: ${dest_dir}"
echo "Source directory: ${src_dir}"

# Create target directory if it doesn't exist
mkdir -p "${dest_dir}"

cell_abund_csv="cell_abundances_and_clusters.csv"

# Check if cell_abundances_and_clusters.csv exists in spatial_output
if [ -f "${src_dir}/${cell_abund_csv}" ]; then
    echo "Found ${cell_abund_csv}, moving to output directory..."
    mv "${src_dir}/${cell_abund_csv}" "${dest_dir}/${cell_abund_csv}"
    echo "File moved successfully."
# If not found, check for merged file
elif [ -f "${src_dir}/merged_results/merged_${cell_abund_csv}" ]; then
    echo "Found merged_${cell_abund_csv}, moving to output directory..."
    mv "${src_dir}/merged_results/merged_${cell_abund_csv}" "${dest_dir}/${cell_abund_csv}"
    echo "Merged file moved successfully."
else
    echo "Neither ${cell_abund_csv} nor merged_${cell_abund_csv} found."
    exit 1
fi
