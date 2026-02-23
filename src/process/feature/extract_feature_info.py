import os
import pandas as pd

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


def process_all_feature_files(feature_folder, output_file):
    
    feature_files = [f for f in os.listdir(feature_folder) if f.endswith(".feature")]

    for i, fname in enumerate(feature_files, 1):
        fpath = os.path.join(feature_folder, fname)
        print(f"[{i}/{len(feature_files)}] Processing {fname}")
        project_id = fname.split("_")[0]

        df = process_feature_file(fpath, project_id)
        # add project id 
        df["PROJECT_ID"] = project_id
        # save to each file
        output_filename = fname.split('.')[0] + '_feature.tsv'
        output_file = os.path.join(feature_folder, output_filename)
        df.to_csv(output_file, sep="\t", index=False)
       
       
    # print(f"Saved to: {output_file}")
