import os
import pandas as pd
import sys

def process_feature_file(feature_file):
    temp_df = pd.read_csv(feature_file, sep="\t")
    temp_df['File_name'] = temp_df['File_name'].apply(lambda x: os.path.basename(str(x)))

    temp_df = temp_df.sort_values(
        by=['File_name', 'Scans', 'Precursor_intensity'],
        ascending=[True, True, False]
    )

    def agg_func(group):
        return pd.Series({
            "feature_id": ":".join(group['Fraction_feature_ID'].astype(str)),

            "feature_intensity": ":".join(
                group['Fraction_feature_intensity'].map(lambda x: f"{x:.2f}")
            ),

            "precursor_intensity": ":".join(
                group['Precursor_intensity'].map(lambda x: f"{x:.2f}")
            ),

            "feature_score": ":".join(
                group['Fraction_feature_score'].map(lambda x: f"{x:.5f}")
            ),

            "feature_apex_time": ":".join(
                group['Fraction_feature_apex_time'].map(lambda x: f"{x:.2f}")
            )
        })

    result_df = (
        temp_df
        .groupby(['File_name', 'Scans'], as_index=False)
        .apply(agg_func, include_groups=False)
        .reset_index(drop=True)
    )

    return result_df


def write_feature_file(dataset_id, feature_file, output_filename):
    df = process_feature_file(feature_file)
    df["DATASET_id"] = dataset_id
    # save the file
    df.to_csv(output_filename, sep="\t", index=False)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <dataset_id> <feature_filename> <output_tsv_filename>")
    else:
        dataset_id = sys.argv[1]
        feature_file = sys.argv[2]
        out_filename = sys.argv[3]
        write_feature_file(dataset_id, feature_file, out_filename)
