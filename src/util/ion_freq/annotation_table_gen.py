import pandas as pd
from pathlib import Path
import argparse


def count_frequency(df_termins, anno_df):
    df_termins = df_termins.copy()
    count_map = dict(
        zip(
            anno_df['Ion'],
            anno_df['Count']
        )
    )
    freq_map = dict(
        zip(
            anno_df['Ion'],
            anno_df['Coverage']
        )
    )
    df_termins['count'] = (
       df_termins['label'].map(count_map)
    )
    df_termins['coverage'] = (
       df_termins['label'].map(freq_map)
    )
    return df_termins


def anno_table_gen(input_tsv, table_tsv, out_prefix, activation_type):
    df = pd.read_csv(table_tsv, sep='\t')

    if activation_type.lower() == "etd":
        # remove rows labeled 'other' 
        df = df[df['note'] != 'other']
    else:
        # Remove rows labeled 'etd'
        df = df[df['note'] != 'etd']
        
    
    anno_df = pd.read_csv(input_tsv, sep='\t')
    
    nc_results = count_frequency(df, anno_df)
    nc_results = nc_results.drop(columns='note')
    nc_results['activation'] = activation_type.upper()
  
    nc_results = nc_results[['activation'] + [col for col in nc_results.columns if col != 'activation']]
    # Save outputs
    out_file = Path(out_prefix) / f"{activation_type}_anno_table.tsv"
    out_file.parent.mkdir(parents=True, exist_ok=True)
    
    nc_results.to_csv(out_file, sep='\t', index=False)
    
    print(f"Saved: {out_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate annotation table for N and C-terminal ions")
    parser.add_argument("--input", required=True, type=str, help="Input TSV with ion counts")
    parser.add_argument("--table", required=True, type=str, help="A TSV table with selected ion types", default="basic_annotation_table.tsv")
    parser.add_argument("--out", required=True, type=str, help="Output folder for annotated TSVs")
    parser.add_argument("--activation", required=True, type=str, help="Activation type (hcd, cid, etd, etc.)")
    args = parser.parse_args()

    anno_table_gen(args.input, args.table, args.out, args.activation)
