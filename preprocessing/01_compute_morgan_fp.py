import argparse
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import pandas as pd
import numpy as np

RADIUS = 2
N_BITS = 512

def smiles_to_morgan_counts(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.zeros(N_BITS, dtype=int)

    fp = rdMolDescriptors.GetHashedMorganFingerprint(
        mol, radius=RADIUS, nBits=N_BITS
    )

    arr = np.zeros(N_BITS, dtype=int)
    for idx, value in fp.GetNonzeroElements().items():
        arr[idx] = value
    return arr


def main(args):
    df = pd.read_csv(args.input)

    # Optional ID creation (external only)
    if args.add_id:
        if "id" in df.columns:
            raise ValueError("ID column already exists")
        df.insert(0, "id", range(args.start_id, args.start_id + len(df)))

    fps = df["SMILES"].apply(smiles_to_morgan_counts)
    fp_df = pd.DataFrame(
        fps.tolist(),
        columns=[f"fp_{i}" for i in range(N_BITS)]
    )

    base_cols = [c for c in ["id", "SMILES", "Tm"] if c in df.columns]
    df_out = pd.concat([df[base_cols], fp_df], axis=1)

    df_out.to_csv(args.output, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--add-id", action="store_true")
    parser.add_argument("--start-id", type=int, default=3330)

    args = parser.parse_args()
    main(args)