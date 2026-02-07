
import pandas as pd

# Load CSV
df_train = pd.read_csv("melt_train_morgan_fps_512.csv")
df_ext = pd.read_csv("external_morgan_fps_512.csv")

# Sanity check
assert list(df_train.columns) == list(df_ext.columns)