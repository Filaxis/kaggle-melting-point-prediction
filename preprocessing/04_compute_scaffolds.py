import argparse
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

# -----------------------------
# Murcko scaffold helper
# -----------------------------
def compute_scaffold(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "INVALID"      # invalid or unparsable SMILES

    scaffold = MurckoScaffold.MurckoScaffoldSmiles(
        mol=mol,
        includeChirality=False  # group by 2D scaffold only
    )

    if scaffold == "":
        return "NO_SCAFFOLD"   # valid molecule, no Murcko core

    return scaffold

def main(args):
    df = pd.read_csv(args.input)
    df["scaffold"] = df["SMILES"].apply(compute_scaffold)
    df.to_csv(args.output, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    main(parser.parse_args())