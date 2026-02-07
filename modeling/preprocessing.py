from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import QuantileTransformer


def make_preprocessor(fp_cols, desc_cols):
    return ColumnTransformer(
        transformers=[
            ("fp", "passthrough", fp_cols),
            (
                "descriptors",
                QuantileTransformer(
                    output_distribution="normal",
                    n_quantiles=1000,
                    random_state=42,
                ),
                desc_cols,
            ),
        ],
        remainder="drop",
    )