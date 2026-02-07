from xgboost import XGBRegressor


BASE_XGB_PARAMS = dict(
    n_estimators=800,
    max_depth=6,
    learning_rate=0.05,
    min_child_weight=10,
    subsample=0.8,
    colsample_bytree=0.8,
    reg_alpha=0.1,
    reg_lambda=1.0,
    objective="reg:absoluteerror",
    n_jobs=-1,
)


def make_xgb(seed):
    return XGBRegressor(**BASE_XGB_PARAMS, random_state=seed)