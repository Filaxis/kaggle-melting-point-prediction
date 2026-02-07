import argparse
import pandas as pd
from rdkit import Chem


def canonical_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, canonical=True)


def main(args):
    # Load datasets
    df_train = pd.read_csv(args.train)
    df_test = pd.read_csv(args.test)
    df_ext = pd.read_csv(args.external)

    # Canonicalize SMILES
    df_train["canon_smiles"] = df_train["SMILES"].apply(canonical_smiles)
    df_test["canon_smiles"] = df_test["SMILES"].apply(canonical_smiles)
    df_ext["canon_smiles"] = df_ext["SMILES"].apply(canonical_smiles)

    # Drop invalid external SMILES
    df_ext = df_ext[df_ext["canon_smiles"].notnull()]

    # Remove external duplicates vs train
    train_smiles = set(df_train["canon_smiles"].dropna())
    df_ext = df_ext[~df_ext["canon_smiles"].isin(train_smiles)]

    # Remove duplicates within external
    df_ext = df_ext.drop_duplicates(subset="canon_smiles")

    # Remove external duplicates vs test
    test_smiles = set(df_test["canon_smiles"].dropna())
    df_ext = df_ext[~df_ext["canon_smiles"].isin(test_smiles)]

    print(f"External dataset size after deduplication: {len(df_ext)}")

    # Drop helper column
    df_ext = df_ext.drop(columns=["canon_smiles"])

    # Save
    df_ext.to_csv(args.output, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Deduplicate external dataset against train and test using canonical SMILES"
    )

    parser.add_argument("--train", required=True, help="Processed Kaggle train CSV")
    parser.add_argument("--test", required=True, help="Processed Kaggle test CSV")
    parser.add_argument("--external", required=True, help="Processed external CSV")
    parser.add_argument("--output", required=True, help="Deduplicated external output CSV")

    args = parser.parse_args()
    main(args)