import pandas as pd
from pathlib import Path
import argparse

def combine_tsv_tables(base_path, output_filename):
    base_path = Path(base_path)
    tsv_files = [f for f in base_path.glob("*_anno_table.tsv")]
    dataframes = []
    for file in tsv_files:
        df = pd.read_csv(file, sep="\t")
        dataframes.append(df)
        
    combined_df_all = pd.concat(dataframes, ignore_index=True)
    combined_df_all.to_csv(output_filename, sep='\t', index=False)     
    print(f"Combined file saved to: {output_filename}")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="combine four type of activation TSV tables into a TSV file.")
    parser.add_argument("--input", required=True, help="Input folder containing TSV files")
    parser.add_argument("--output", required=True, help="Output TSV file path", default="combine_anno_table.tsv")

    args = parser.parse_args()
    combine_tsv_tables(args.input, args.output)

