import numpy as np
from sklearn.pipeline import Pipeline
from sklearn.metrics import mean_absolute_error

def run_seed_ensemble_cv(
    X,
    y,
    groups,
    preprocessor,
    model_factory,
    seeds,
    cv,
    sample_weight=None,
):
    """
    Runs GroupKFold CV for a seed-based ensemble.
    Returns mean and std of MAE across seeds.
    """

    seed_maes = []

    for seed in seeds:
        fold_maes = []

        for train_idx, val_idx in cv.split(X, y, groups):
            X_tr, X_val = X.iloc[train_idx], X.iloc[val_idx]
            y_tr, y_val = y.iloc[train_idx], y.iloc[val_idx]

            if sample_weight is not None:
                sw_tr = sample_weight[train_idx]
            else:
                sw_tr = None

            model = model_factory(seed)
            pipe = Pipeline([
                ("features", preprocessor),
                ("model", model),
            ])

            if sw_tr is not None:
                pipe.fit(X_tr, y_tr, model__sample_weight=sw_tr)
            else:
                pipe.fit(X_tr, y_tr)

            preds = pipe.predict(X_val)
            fold_maes.append(mean_absolute_error(y_val, preds))

        seed_maes.append(np.mean(fold_maes))

    return float(np.mean(seed_maes)), float(np.std(seed_maes))