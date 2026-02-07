import argparse
from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd

descriptor_fns = [
    ("MolWt", Descriptors.MolWt),
    ("MolLogP", Descriptors.MolLogP),
    ("TPSA", Descriptors.TPSA),
    ("NumHDonors", Descriptors.NumHDonors),
    ("NumHAcceptors", Descriptors.NumHAcceptors),
    ("NumRotatableBonds", Descriptors.NumRotatableBonds),
    ("RingCount", Descriptors.RingCount),
    ("FractionCSP3", Descriptors.FractionCSP3),
]

def compute_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return [None] * len(descriptor_fns)
    return [fn(mol) for _, fn in descriptor_fns]


def main(args):
    df = pd.read_csv(args.input)

    desc_df = df["SMILES"].apply(compute_descriptors).apply(pd.Series)
    desc_df.columns = [name for name, _ in descriptor_fns]

    df_out = pd.concat([df, desc_df], axis=1)
    df_out.to_csv(args.output, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)

    args = parser.parse_args()
    main(args)