import os
import pandas as pd
import sys


def feature_merge(msalign_file, feature_file, out_filename=None, wfile=False):

    msalign_df = pd.read_csv(msalign_file, sep="\t")
    msalign_df["FILE_NAME"] = msalign_df["FILE_NAME"].astype(str)
    msalign_df["MS2_SCANS"] = msalign_df["MS2_SCANS"].astype(str)
    msalign_df["DATASET_ID"] = msalign_df["DATASET_ID"].astype(str)

    feature_df = pd.read_csv(feature_file, sep="\t")
    # rename 
    rename_map = {
            "File_name": "FILE_NAME",
            "Scans": "MS2_SCANS",
            "DATASET_id": "DATASET_ID"
        }
    feature_df = feature_df.rename(columns=rename_map)

    feature_df["FILE_NAME"] = feature_df["FILE_NAME"].astype(str)
    feature_df["MS2_SCANS"] = feature_df["MS2_SCANS"].astype(str)
    feature_df["DATASET_ID"] = feature_df["DATASET_ID"].astype(str)

    feature_df_merged = msalign_df.merge(
        feature_df,
        on=["DATASET_ID", "FILE_NAME", "MS2_SCANS"],
        how="left"
    )
    # Save merged file
    if wfile and out_filename:
        feature_df_merged.to_csv(out_filename, sep='\t', index=False)
        print(f"Merged file saved to: {out_filename}")
    return feature_df_merged

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <msalign_info_filename> <feature_info_filename> <output_tsv_filename>")
    else:
        msalign_info_tsv = sys.argv[1]
        feature_info_tsv = sys.argv[2]
        out_filename = sys.argv[3]
        df = feature_merge(msalign_info_tsv, feature_info_tsv, out_filename, wfile=True)
