from sklearn.preprocessing import LabelEncoder


def encode_scaffolds(train_df, test_df, col="scaffold"):
    """
    Fit scaffold encoder on training data only.
    """
    le = LabelEncoder()
    train_df["scaffold_id"] = le.fit_transform(train_df[col])

    mapping = dict(zip(le.classes_, le.transform(le.classes_)))

    test_df["scaffold_id"] = (
        test_df[col]
        .map(mapping)
        .fillna(-1)
        .astype(int)
    )

    return train_df, test_df