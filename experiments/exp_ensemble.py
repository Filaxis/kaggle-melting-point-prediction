import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import Pipeline

from features.feature_config import FP_COLS, DESC_COLS
from features.scaffold_features import encode_scaffolds
from modeling.preprocessing import make_preprocessor
from modeling.xgb_models import make_xgb
from modeling.ensemble import run_seed_ensemble_cv

# -------------------------
# Data loading
# -------------------------

def load_data():
    data_dir = Path("data/03_processed")
    train = pd.read_csv(data_dir / "train_merged_fp_desc.csv")
    test = pd.read_csv(data_dir / "test_fp_desc_scaffold.csv")
    return train, test

# -------------------------
# Main experiment
# -------------------------

def main():
    df_train, df_test = load_data()

    # Encode scaffolds
    df_train, df_test = encode_scaffolds(df_train, df_test)

    X_train = df_train[FP_COLS + DESC_COLS]
    y_train = df_train["Tm"]
    X_test = df_test[FP_COLS + DESC_COLS]

    # Sample weights (external data)
    sample_weight = df_train["is_external"].map({
        0: 1.0,
        1: 0.3,
    }).values

    preprocessor = make_preprocessor(FP_COLS, DESC_COLS)

    seeds = [42, 1337, 2023, 7, 99]
    gkf = GroupKFold(n_splits=5)

    mean_mae, std_mae = run_seed_ensemble_cv(
        X=X_train,
        y=y_train,
        groups=df_train["scaffold_id"],
        preprocessor=preprocessor,
        model_factory=make_xgb,
        seeds=seeds,
        cv=gkf,
        sample_weight=sample_weight,
    )

    print(f"Scaffold CV MAE: {mean_mae:.2f} ± {std_mae:.2f}")

    # -------------------------
    # Train ensemble + predict
    # -------------------------

    preds = np.zeros(len(X_test))

    for seed in seeds:
        model = make_xgb(seed)
        pipe = Pipeline([
            ("features", preprocessor),
            ("model", model),
        ])

        pipe.fit(
            X_train,
            y_train,
            model__sample_weight=sample_weight,
        )

        preds += pipe.predict(X_test)

    preds /= len(seeds)

    submission = pd.DataFrame({
        "id": df_test["id"],
        "Tm": preds,
    })

    submission.to_csv("my_submission.csv", index=False)
    print("Submission saved.")


if __name__ == "__main__":
    main()