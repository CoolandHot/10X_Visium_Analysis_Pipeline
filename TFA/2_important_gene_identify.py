import pandas as pd
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
import yaml

config = yaml.load(open("config/batch_config.yaml"), Loader=yaml.FullLoader)
batch_names = config["batch_names"]

cell_types = pd.read_csv("./TFA/output/ordered_cell_types.csv", index_col=0)
for batch_id in batch_names:
    TFA = pd.read_csv(f"./TFA/output/{batch_id}_TF_activities.csv", index_col=0)
    TFA["cell_types"] = TFA["cell_types"].map(cell_types.x)

    top_features_by_cell_type_gb = {}
    top_features_by_cell_type_rf = {}
    unique_cell_types = TFA["cell_types"].unique()

    for cell_type in unique_cell_types:
        if pd.isna(cell_type):
            continue
        # Create binary labels: 1 for the current cell type, 0 for all others
        TFA["binary_target"] = (TFA["cell_types"] == cell_type).astype(int)

        X = TFA.drop(columns=["cell_types", "binary_target"])
        y = TFA["binary_target"]

        rf = RandomForestClassifier(random_state=42)
        rf.fit(X, y)
        gbc = GradientBoostingClassifier(random_state=42)
        gbc.fit(X, y)

        top_features_by_cell_type_gb[cell_type] = pd.Series(
            gbc.feature_importances_, index=X.columns
        ).sort_values(ascending=False)
        top_features_by_cell_type_rf[cell_type] = pd.Series(
            rf.feature_importances_, index=X.columns
        ).sort_values(ascending=False)

    # export to csv
    pd.DataFrame(
        [
            {"cell_type": cell_type, "feature": feature, "importance": importance}
            for cell_type, features in top_features_by_cell_type_gb.items()
            for feature, importance in features.items()
        ]
    ).to_csv(
        f"./TFA/output/{batch_id}_top_features_by_cell_type_GradientBoosting.csv",
        index=False,
    )

    pd.DataFrame(
        [
            {"cell_type": cell_type, "feature": feature, "importance": importance}
            for cell_type, features in top_features_by_cell_type_rf.items()
            for feature, importance in features.items()
        ]
    ).to_csv(
        f"./TFA/output/{batch_id}_top_features_by_cell_type_RandomForest.csv",
        index=False,
    )
