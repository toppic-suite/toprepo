import os
import pandas as pd
import sys


def add_project_id(file_name, feature_df):
    fname = os.path.basename(file_name)
    project_id = fname.split("_ms2")[0].split("_")[0]
    feature_df["PROJECT_ID"] = project_id
    return feature_df


def feature_merge(msalign_file, feature_file, out_filename):

    msalign_df = pd.read_csv(msalign_file, sep="\t")
    msalign_df["FILE_NAME"] = msalign_df["FILE_NAME"].astype(str)
    msalign_df["MS2_SCANS"] = msalign_df["MS2_SCANS"].astype(str)
    msalign_df["PROJECT_ID"] = msalign_df["PROJECT_ID"].astype(str)

    feature_df = pd.read_csv(feature_file, sep="\t")
    feature_df = add_project_id(feature_file, feature_df) # check if it has project_id column
    # rename 
    rename_map = {
            "File_name": "FILE_NAME",
            "Scans": "MS2_SCANS",
        }
    feature_df = feature_df.rename(columns=rename_map)

    # normalize
    feature_df["FILE_NAME"] = feature_df["FILE_NAME"].astype(str)
    feature_df["MS2_SCANS"] = feature_df["MS2_SCANS"].astype(str)
    feature_df["PROJECT_ID"] = feature_df["PROJECT_ID"].astype(str)

    feature_df_merged = msalign_df.merge(
        feature_df,
        on=["PROJECT_ID", "FILE_NAME", "MS2_SCANS"],
        how="left"
    )
    # Save merged file
    feature_df_merged.to_csv(out_filename, sep='\t', index=False)
    print(f"Merged file saved to: {out_filename}")
   

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <msalign_info_filename> <feature_info_filename> <output_tsv_filename>")
    else:
        msalign_info_tsv = sys.argv[1]
        feature_info_tsv = sys.argv[2]
        out_filename = sys.argv[3]
        feature_merge(msalign_info_tsv, feature_info_tsv, out_filename)
