#!/bin/bash

# Script to distribute cell2location batch processing across otter08 to otter20 hosts

# Configuration
SOURCE_SCRIPT="/scratch/hh01116/10X_KD/cellType/pred_cell2location_2_map_spatial.py"
UTIL_SCRIPT="/scratch/hh01116/10X_KD/cellType/pred_cell2location_utils.py"
SOURCE_CONFIG="/scratch/hh01116/10X_KD/config"
DEST_BASE_DIR="/scratch/hh01116/10X_KD_infer"
SPLIT_DATASET_DIR="/scratch/hh01116/10X_KD/cellType/gbm_analysis/split_datasets"
TMUX_SESSION_NAME="cell2loc_batch"

process_batch() {
    local BATCH_NUM="$1"
    local HOSTNAME="$2"

    # Validate batch file exists
    local BATCH_FILE="${SPLIT_DATASET_DIR}/batch_${BATCH_NUM}.h5ad"
    if [ ! -f "$BATCH_FILE" ]; then
        echo "Error: Batch file $BATCH_FILE does not exist"
        return 1
    fi

    echo "Targeting $HOSTNAME, processing batch $BATCH_NUM"

    # Create destination directory on remote host
    echo "Creating destination directory on $HOSTNAME: $DEST_BASE_DIR"
    ssh "$HOSTNAME" "mkdir -p '$DEST_BASE_DIR/cellType/gbm_analysis/split_datasets'"
    ssh "$HOSTNAME" "mkdir -p '$DEST_BASE_DIR/config'"

    # Copy files to remote host
    echo "Copying $BATCH_FILE to $HOSTNAME"
    ssh "$HOSTNAME" "[ -f '$DEST_BASE_DIR/cellType/gbm_analysis/split_datasets/batch_${BATCH_NUM}.h5ad' ] || exit 1"
    if [ $? -ne 0 ]; then
        scp "$BATCH_FILE" "$HOSTNAME:$DEST_BASE_DIR/cellType/gbm_analysis/split_datasets/"
    else
        echo "$BATCH_FILE already exists on $HOSTNAME, skipping copy."
    fi
    scp "${SPLIT_DATASET_DIR}/batch_info.json" "$HOSTNAME:$DEST_BASE_DIR/cellType/gbm_analysis/split_datasets/"

    echo "Copying $SOURCE_SCRIPT & $UTIL_SCRIPT to $HOSTNAME"
    scp "$SOURCE_SCRIPT" "$HOSTNAME:$DEST_BASE_DIR/cellType/"
    scp "$UTIL_SCRIPT" "$HOSTNAME:$DEST_BASE_DIR/cellType/"

    echo "Copying $SOURCE_CONFIG to $HOSTNAME"
    scp -r "$SOURCE_CONFIG"/* "$HOSTNAME:$DEST_BASE_DIR/config/"

    # Execute commands on remote host
    echo "Starting processing on $HOSTNAME"
    ssh "$HOSTNAME" "cd '$DEST_BASE_DIR' && tmux new-session -d -s '${TMUX_SESSION_NAME}_${BATCH_NUM}'"
    ssh "$HOSTNAME" "tmux send-keys -t '${TMUX_SESSION_NAME}_${BATCH_NUM}' 'cd $DEST_BASE_DIR' Enter"
    ssh "$HOSTNAME" "tmux send-keys -t '${TMUX_SESSION_NAME}_${BATCH_NUM}' 'apptainer exec --nv --bind /usr/lib/locale:/usr/lib/locale --bind /vol/research/brainTumorST:/vol/research/brainTumorST --bind /vol/research/scratch1/NOBACKUP/hh01116/:/vol/research/scratch1/NOBACKUP/hh01116/ --bind /scratch/hh01116/:/scratch/hh01116/ --env R_LIBS_USER=\"\" --env R_ENVIRON_USER=\"\" /vol/research/scratch1/NOBACKUP/hh01116/docker_imgs/cell2location.sif bash' Enter"
    ssh "$HOSTNAME" "tmux send-keys -t '${TMUX_SESSION_NAME}_${BATCH_NUM}' 'python cellType/pred_cell2location_2_map_spatial.py --batch $BATCH_NUM' Enter"

    echo "Started batch $BATCH_NUM processing on $HOSTNAME in tmux session ${TMUX_SESSION_NAME}_${BATCH_NUM}"
    echo "To monitor progress, attach to the session with:"
    echo "ssh -t $HOSTNAME 'tmux attach'"
}

declare -A batch_host_dict=( ["0"]="otter08" ["1"]="otter12"  ["2"]="otter21"  ["3"]="otter49"  ["4"]="otter39"  ["5"]="otter60"  ["6"]="otter04"  ["7"]="otter66"  ["8"]="otter79"  ["9"]="otter50"  ["10"]="otter48"  ["11"]="otter47"  ["12"]="otter51"  ["13"]="otter65"  ["14"]="otter64" )
for batch in "${!batch_host_dict[@]}"; do
    process_batch "$batch" "${batch_host_dict[$batch]}"
done

batch_host_dict_finish=( ["0"]="otter08" ["1"]="otter12"  ["2"]="otter21"  ["3"]="otter49"  ["4"]="otter39"  ["5"]="otter60"  ["7"]="otter66"  ["9"]="otter50"  ["10"]="otter48"  ["11"]="otter47"  ["12"]="otter51"  ["13"]="otter65"  ["14"]="otter64" )

ssh -t otter79 'tmux attach'
ssh -t otter04 'tmux attach'
exit

# retrieve the results
echo "Retrieving results from remote hosts..."
for batch in "${!batch_host_dict_finish[@]}"; do
    echo "Copying results from ${batch_host_dict_finish[$batch]}..."
    scp -r "${batch_host_dict_finish[$batch]}:$DEST_BASE_DIR/cellType/gbm_analysis/cell2location_map/batch_$batch" "/scratch/hh01116/10X_KD/cellType/gbm_analysis/cell2location_map/"
done