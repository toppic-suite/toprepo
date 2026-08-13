import os
import pandas as pd
import numpy as np
import argparse

def read_file_remove_comments(file_path):
    with open(file_path, "r") as f:
        for i, line in enumerate(f):
            if "Data file name" in line:
                header_line = i
                break

    df = pd.read_csv(
        file_path,
        sep="\t",
        skiprows=header_line
    )
    return df


def main():
    parser = argparse.ArgumentParser(
        description="Read a TSV file and skip comment lines before the header."
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input TSV file"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV file"
    )
    args = parser.parse_args()

    df = read_file_remove_comments(args.input)
    df.to_csv(args.output, sep="\t", index=False)

    print(f"Read {len(df):,} rows")
    print(f"Saved to: {args.output}")

if __name__ == "__main__":
    main()
