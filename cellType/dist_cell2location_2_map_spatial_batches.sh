#!/bin/bash

# Script to distribute cell2location batch processing across otter08 to otter20 hosts

# Configuration
SOURCE_SCRIPT="/scratch/hh01116/10X_KD/cellType/pred_cell2location_2_map_spatial.py"
UTIL_SCRIPT="/scratch/hh01116/10X_KD/cellType/pred_cell2location_utils.py"
SOURCE_CONFIG="/scratch/hh01116/10X_KD/config"

DEST_BASE_DIR="/vol/research/scratch1/NOBACKUP/hh01116/10X_KD_infer"
SPLIT_DATASET_DIR="/scratch/hh01116/10X_KD/cellType/gbm_analysis/split_datasets"
TMUX_SESSION_NAME="cell2loc_batch"

cp "$SOURCE_SCRIPT" "$DEST_BASE_DIR/cellType"
cp "$UTIL_SCRIPT" "$DEST_BASE_DIR/cellType"
cp -r "$SOURCE_CONFIG" "$DEST_BASE_DIR"
cp -r "$SPLIT_DATASET_DIR" "$DEST_BASE_DIR/cellType/gbm_analysis/"


process_batch() {
    local BATCH_NUM="$1"
    local HOSTNAME="$2"

    # Validate batch file exists
    local BATCH_FILE="${DEST_BASE_DIR}/cellType/gbm_analysis/split_datasets/batch_${BATCH_NUM}.h5ad"
    if [ ! -f "$BATCH_FILE" ]; then
        echo "Error: Batch file $BATCH_FILE does not exist"
        return 1
    fi

    echo "Targeting $HOSTNAME, processing batch $BATCH_NUM"

    # Execute commands on remote host
    echo "Starting processing on $HOSTNAME"
    ssh "$HOSTNAME" "cd '$DEST_BASE_DIR' && tmux new-session -d -s '${TMUX_SESSION_NAME}_${BATCH_NUM}'"
    # ssh "$HOSTNAME" "tmux send-keys -t '${TMUX_SESSION_NAME}_${BATCH_NUM}' 'cd $DEST_BASE_DIR' Enter"
    ssh "$HOSTNAME" "tmux send-keys -t '${TMUX_SESSION_NAME}_${BATCH_NUM}' 'apptainer exec --nv --bind /usr/lib/locale:/usr/lib/locale --bind /vol/research/brainTumorST:/vol/research/brainTumorST --bind /vol/research/scratch1/NOBACKUP/hh01116/:/vol/research/scratch1/NOBACKUP/hh01116/ --env R_LIBS_USER=\"\" --env R_ENVIRON_USER=\"\" /vol/research/scratch1/NOBACKUP/hh01116/docker_imgs/cell2location.sif bash' Enter"
    ssh "$HOSTNAME" "tmux send-keys -t '${TMUX_SESSION_NAME}_${BATCH_NUM}' 'python cellType/pred_cell2location_2_map_spatial.py --batch $BATCH_NUM' Enter"

    echo "Started batch $BATCH_NUM processing on $HOSTNAME in tmux session ${TMUX_SESSION_NAME}_${BATCH_NUM}"
    echo "To monitor progress, attach to the session with:"
    echo "ssh -t $HOSTNAME 'tmux attach'"
}

declare -A batch_host_dict=( ["0"]="otter08" ["1"]="otter12"  ["2"]="otter21"  ["3"]="otter49"  ["4"]="otter39"  ["5"]="otter60"  ["6"]="otter04"  ["7"]="otter66"  ["8"]="otter79"  ["9"]="otter50"  ["10"]="otter48"  ["11"]="otter47"  ["12"]="otter51"  ["13"]="otter65"  ["14"]="otter64" )
for batch in "${!batch_host_dict[@]}"; do
    process_batch "$batch" "${batch_host_dict[$batch]}"
done
