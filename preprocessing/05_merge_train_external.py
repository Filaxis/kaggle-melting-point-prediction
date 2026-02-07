import argparse
import pandas as pd

def main(args):
    df_train = pd.read_csv(args.train)
    df_ext   = pd.read_csv(args.external)

    df_train["is_external"] = 0
    df_ext["is_external"]   = 1

    df_merged = pd.concat([df_train, df_ext], ignore_index=True)
    df_merged.to_csv(args.output, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--train", required=True, help="Processed Kaggle train CSV")
    parser.add_argument("--external", required=True, help="Processed external CSV")
    parser.add_argument("--output", required=True, help="Merged train and external output CSV")
    main(parser.parse_args())