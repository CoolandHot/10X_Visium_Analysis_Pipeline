#!/bin/bash

# bash create_imgLinks_to_spatials.sh unlink
# bash create_imgLinks_to_spatials.sh create

parent_folders=(
    "/vol/research/brainTumorST/10X_Genomics_KD/Run 1 - 22-08-24"
    "/vol/research/brainTumorST/10X_Genomics_KD/Run 2 - 11-10-24"
    "/vol/research/brainTumorST/10X_Genomics_KD/Run 3 - 14-11-24"
)

# Function to create symbolic links
create_links() {
    local parent_folder0="$1"

    for sample in $(ls -1 "$parent_folder0" | grep -E '^[0-9]'); do
        parent_folder1="${parent_folder0}/${sample}"
        parent_folder="${parent_folder1}/outs"
        source_dir="${parent_folder}/spatial"

        target_dirs=(
            "${parent_folder}/binned_outputs/square_002um/spatial"
            "${parent_folder}/binned_outputs/square_008um/spatial"
            "${parent_folder}/binned_outputs/square_016um/spatial"
        )

        files=(
            "tissue_lowres_image.png"
            "tissue_hires_image.png"
        )

        for target_dir in "${target_dirs[@]}"; do
            for file in "${files[@]}"; do
                # Calculate the relative path from the target directory to the source file
                relative_path=$(realpath --relative-to="${target_dir}" "${source_dir}")
                ln -sf "${relative_path}/${file}" "${target_dir}/${file}"
            done
        done
    done

    echo "Symbolic links created successfully."
}

# Function to remove symbolic links
unlink_files() {
    local parent_folder0="$1"

    for sample in $(ls -1 "$parent_folder0" | grep -E '^[0-9]'); do
        parent_folder1="${parent_folder0}/${sample}"
        parent_folder="${parent_folder1}/outs"

        target_dirs=(
            "${parent_folder}/binned_outputs/square_002um/spatial"
            "${parent_folder}/binned_outputs/square_008um/spatial"
            "${parent_folder}/binned_outputs/square_016um/spatial"
        )

        files=(
            "tissue_lowres_image.png"
            "tissue_hires_image.png"
        )

        for target_dir in "${target_dirs[@]}"; do
            for file in "${files[@]}"; do
                if [ -L "${target_dir}/${file}" ]; then
                    rm "${target_dir}/${file}"
                    echo "Removed symbolic link: ${target_dir}/${file}"
                fi
            done
        done
    done

    echo "Symbolic links removed successfully."
}

# Main execution
case "${1:-create}" in
    "create")
        for parent_folder0 in "${parent_folders[@]}"; do
            create_links "$parent_folder0"
        done
        ;;
    "unlink")
        for parent_folder0 in "${parent_folders[@]}"; do
            unlink_files "$parent_folder0"
        done
        ;;
    *)
        echo "Usage: $0 [create|unlink]"
        echo "  create - Create symbolic links (default)"
        echo "  unlink - Remove symbolic links"
        exit 1
        ;;
esac
