# Utility functions for cell2location analysis pipelines.


import os
import pickle
import yaml
import scanpy as sc


def load_configuration(config_path='config/cellType_config.yaml'):
    """Loads configuration from a YAML file."""
    with open(config_path, 'r') as f:
        config_data = yaml.safe_load(f)
    return config_data


def extract_adata_additions(adata_before, adata_after, output_path):
    """
    Extract adata additions from the after- AnnData object and save to file.
    """
    print("\n--- Extracting additions ---")

    obs_before = set(adata_before.obs.columns)
    obs_after = set(adata_after.obs.columns)
    adata_obs_cols = sorted(obs_after - obs_before)

    uns_before = set(adata_before.uns.keys())
    uns_after = set(adata_after.uns.keys())
    adata_uns_keys = sorted(uns_after - uns_before)

    obsm_before = set(adata_before.obsm.keys())
    obsm_after = set(adata_after.obsm.keys())
    adata_obsm_keys = sorted(obsm_after - obsm_before)

    varm_before = set(adata_before.varm.keys())
    varm_after = set(adata_after.varm.keys())
    adata_varm_keys = sorted(varm_after - varm_before)

    layers_before = set(adata_before.layers.keys()) if adata_before.layers else set()
    layers_after = set(adata_after.layers.keys()) if adata_after.layers else set()
    adata_layers_keys = sorted(layers_after - layers_before)

    adata_additions = {
        'obs_data': {},
        'uns_data': {},
        'obsm_data': {},
        'varm_data': {},
        'layers_data': {},
        'metadata': {
            'obs_columns': adata_obs_cols,
            'uns_keys': adata_uns_keys,
            'obsm_keys': adata_obsm_keys,
            'varm_keys': adata_varm_keys,
            'layers_keys': adata_layers_keys,
            'obs_names': adata_after.obs_names.tolist(),
            'var_names': adata_after.var_names.tolist()
        }
    }

    for col in adata_obs_cols:
        adata_additions['obs_data'][col] = adata_after.obs[col].values
        print(f"  Extracted .obs column: {col}")

    for key in adata_uns_keys:
        adata_additions['uns_data'][key] = adata_after.uns[key]
        print(f"  Extracted .uns key: {key}")

    for key in adata_obsm_keys:
        adata_additions['obsm_data'][key] = adata_after.obsm[key]
        print(f"  Extracted .obsm key: {key}")

    for key in adata_varm_keys:
        adata_additions['varm_data'][key] = adata_after.varm[key]
        print(f"  Extracted .varm key: {key}")

    for key in adata_layers_keys:
        adata_additions['layers_data'][key] = adata_after.layers[key]
        print(f"  Extracted .layers key: {key}")

    with open(output_path, 'wb') as f:
        pickle.dump(adata_additions, f)

    print(f"Additions saved to: {output_path}")
    print(f"Summary: {len(adata_obs_cols)} .obs columns, {len(adata_uns_keys)} .uns keys, "
          f"{len(adata_obsm_keys)} .obsm keys, {len(adata_varm_keys)} .varm keys, "
          f"{len(adata_layers_keys)} .layers keys")

    return adata_additions


def restore_adata_additions(adata_base, adata_additions_path):
    """
    Read adata additions from file and merge them into a base AnnData object.
    """
    print(f"\n--- Restoring additions from {adata_additions_path} ---")

    with open(adata_additions_path, 'rb') as f:
        adata_additional = pickle.load(f)

    adata_restored = adata_base.copy()

    if (adata_restored.obs_names.tolist() != adata_additional['metadata']['obs_names'] or
            adata_restored.var_names.tolist() != adata_additional['metadata']['var_names']):
        print("Warning: obs_names or var_names don't match between base AnnData and saved additions.")
        print("This may cause issues. Proceeding anyway...")

    for col in adata_additional['metadata']['obs_columns']:
        adata_restored.obs[col] = adata_additional['obs_data'][col]
        print(f"  Restored .obs column: {col}")

    for key in adata_additional['metadata']['uns_keys']:
        adata_restored.uns[key] = adata_additional['uns_data'][key]
        print(f"  Restored .uns key: {key}")

    for key in adata_additional['metadata']['obsm_keys']:
        adata_restored.obsm[key] = adata_additional['obsm_data'][key]
        print(f"  Restored .obsm key: {key}")

    for key in adata_additional['metadata']['varm_keys']:
        adata_restored.varm[key] = adata_additional['varm_data'][key]
        print(f"  Restored .varm key: {key}")

    for key in adata_additional['metadata']['layers_keys']:
        adata_restored.layers[key] = adata_additional['layers_data'][key]
        print(f"  Restored .layers key: {key}")

    print("Successfully restored additions to AnnData object")
    print(f"Summary: {len(adata_additional['metadata']['obs_columns'])} .obs columns, "
          f"{len(adata_additional['metadata']['uns_keys'])} .uns keys, "
          f"{len(adata_additional['metadata']['obsm_keys'])} .obsm keys, "
          f"{len(adata_additional['metadata']['varm_keys'])} .varm keys, "
          f"{len(adata_additional['metadata']['layers_keys'])} .layers keys")

    return adata_restored


def load_incremental_adata(run_name_path, base_adata_filename="sp.h5ad", model_increments_filename="sp_model_increments.pkl", final_increments_filename="sp_final_increments.pkl"):
    """
    Loads an AnnData object from a base file and applies model and final increments if present.
    Returns the fully restored AnnData object, or None if files are missing.
    """
    import scanpy as sc
    import os

    base_adata_path = os.path.join(run_name_path, base_adata_filename)
    model_increments_path = os.path.join(run_name_path, model_increments_filename)
    final_increments_path = os.path.join(run_name_path, final_increments_filename)

    if os.path.exists(final_increments_path) and os.path.exists(base_adata_path):
        print("Found incremental models.")
        adata_base = sc.read_h5ad(base_adata_path)
        if os.path.exists(model_increments_path):
            adata_base = restore_adata_additions(adata_base, model_increments_path)
        adata = restore_adata_additions(adata_base, final_increments_path)
        return adata
    else:
        return None


def load_processed_spatial_data(run_name_path):
    """
    Loads the final processed AnnData object.
    """
    adata = load_incremental_adata(run_name_path)
    if adata is not None:
        pass  # already loaded
    else:
        print("No incremental models found, loading directly from base AnnData.")
        adata_path = [os.path.join(run_name_path, "sp_final.h5ad"),
                      os.path.join(run_name_path, "merged_results", "merged_adata.h5ad")]
        for path in adata_path:
            if os.path.exists(path):
                print(f"Loading processed spatial data from legacy format: {path}")
                adata = sc.read_h5ad(path)
                break
        if adata is None:
            raise FileNotFoundError(f"Processed AnnData file not found at both path: {adata_path}. "
                                    "Please ensure pred_cell2location_2_map_spatial.py has been run successfully.")

    if adata.raw is None:
        print("adata.raw is not set. DGE will be performed on adata.X. "
              "For DGE on raw counts, ensure adata.raw is appropriately populated.")

    return adata
