import pandas as pd
import os
import sys


def add_project_id(input_folder, output_folder, filename):
    os.makedirs(output_folder, exist_ok=True)

    project_df = pd.read_csv(filename, sep='\t')

    tsv_files = [f for f in os.listdir(input_folder) if f.endswith("_combined_info.tsv")]

    total_rows = 0
    fcount = 0

    for i, fname in enumerate(tsv_files, 1):
        fpath = os.path.join(input_folder, fname)
        print(f"[{i}/{len(tsv_files)}] {fname}")

        df = pd.read_csv(fpath, sep="\t", low_memory=False)

        mzml_fname = df['MZML_file_name'].iloc[0]
        dataset_id = df['DATASET_id'].iloc[0]

        match = project_df[
            (project_df['MZML_file_name'] == mzml_fname) &
            (project_df['Dataset_id'] == dataset_id)
        ]

        project_id = match['Project_id'].iloc[0]
        df.insert(0, 'PROJECT_id', project_id)

        out_path = os.path.join(output_folder, fname)
        df.to_csv(out_path, sep="\t", index=False)

        fcount += 1
        total_rows += len(df)

    print(f"\nDone! Total {fcount} files generated with total {total_rows} rows.")


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python script.py <input_folder> <output_folder> <project_table>")
        sys.exit(1)

    input_folder = sys.argv[1]
    output_folder = sys.argv[2]
    tsv_filename = sys.argv[3]

    add_project_id(input_folder, output_folder, tsv_filename)

