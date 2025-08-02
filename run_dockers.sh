podman run --rm -it \
    --name seurat-R \
    -v /vol/research/scratch1/NOBACKUP/hh01116/:/vol/research/scratch1/NOBACKUP/hh01116/ \
    -v /vol/research/brainTumorST:/vol/research/brainTumorST \
    -v /scratch/hh01116/:/scratch/hh01116/ \
    -w $(pwd) \
    container-registry.surrey.ac.uk/shared-containers/bioinformatics-r-seurat5-with-recommended-pkgs

    container-registry.surrey.ac.uk/shared-containers/bioinformatics-r-seurat-5-with-banksy \
    bash -c "Rscript /app/data/3_spatial_domain_clustering.r"
    container-registry.surrey.ac.uk/shared-containers/bioinformatics-r-seurat-5-for-transcription-factor-analysis



apptainer exec --nv --bind /usr/lib/locale:/usr/lib/locale --bind /vol/research/brainTumorST:/vol/research/brainTumorST --bind /vol/research/scratch1/NOBACKUP/hh01116/:/vol/research/scratch1/NOBACKUP/hh01116/ --bind /scratch/hh01116/:/scratch/hh01116/ --env R_LIBS_USER="" --env R_ENVIRON_USER="" \
    /vol/research/scratch1/NOBACKUP/hh01116/docker_imgs/R_seurat5_with_pkgs.sif bash

    /vol/research/scratch1/NOBACKUP/hh01116/docker_imgs/cell2location.sif bash
    
    /vol/research/scratch1/NOBACKUP/hh01116/docker_imgs/R_seurat5_with_banksy.sif radian



export TMPDIR=/scratch/hh01116/tmpdir && \
apptainer pull $VOL_NOBACKUP/docker_imgs/R_seurat5_with_pkgs.sif \
    docker://container-registry.surrey.ac.uk/shared-containers/bioinformatics-r-seurat5-with-recommended-pkgs:latest

export TMPDIR=/scratch/hh01116/tmpdir && \
    apptainer pull $VOL_NOBACKUP/docker_imgs/cell2location.sif \
    oras://container-registry.surrey.ac.uk/shared-containers/apptainer-bioinformatics-cell2location



ssh otter21 'cd /vol/research/scratch1/NOBACKUP/hh01116/Jake_project/cellType && tmux new-session -d "apptainer exec --nv --bind /usr/lib/locale:/usr/lib/locale --bind /vol/research/brainTumorST:/vol/research/brainTumorST --bind /vol/research/scratch1/NOBACKUP/hh01116/:/vol/research/scratch1/NOBACKUP/hh01116/ --bind /scratch/hh01116/:/scratch/hh01116/ --env R_LIBS_USER=\"\" --env R_ENVIRON_USER=\"\" /vol/research/scratch1/NOBACKUP/hh01116/docker_imgs/cell2location.sif bash"'
# Then attach to it
ssh -t otter21 tmux attach


find . -type f \( -name "*.r" -o -name "*.py" -o -name "*.yml" -o -name "*.yaml" -o -name "*.sh" \) \
    -not -path "*/.ipynb_checkpoints/*" \
    -exec sh -c '
    dest_dir=/scratch/hh01116/10X_KD/$(dirname "{}" | sed "s|^\./||");
    ssh otter47 "mkdir -p \"$dest_dir\""
    scp "{}" "otter47:$dest_dir/"
' \;