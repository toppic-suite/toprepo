import sys
import os
import re
import pandas as pd


def msalign_meta_extract(dataset_id, msalign_path):
    records = []
    current = None
    peak_count = 0

    with open(msalign_path, 'r') as f:
        # current = {}
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith("BEGIN IONS"):
                current = {}
                peak_count  = 0
            elif line.startswith("END IONS"):
                # Save only if all required fields are found
                if current:
                    mzml_filename = os.path.basename(current.get("FILE_NAME"))
                    msalign_fullname = os.path.basename(msalign_path)
                    if dataset_id in msalign_fullname:
                        msalign_filename_extract = msalign_fullname.replace(f"{dataset_id}_", "", 1)
                    else:
                        msalign_filename_extract = msalign_fullname
                    records.append({
                        "DATASET_ID": dataset_id,
                        "FILE_NAME": mzml_filename,
                        "SPECTRUM_ID": current.get("SPECTRUM_ID"),
                        "MS2_SCANS": current.get("SCANS"),
                        "MS_ONE_ID": current.get("MS_ONE_ID"),
                        "MS_ONE_SCAN": current.get("MS_ONE_SCAN"),
                        "MSALIGN_FILE_NAME": msalign_filename_extract,
                        "ACTIVATION": current.get("ACTIVATION"),
                        "PRECURSOR_WINDOW_BEGIN": current.get("PRECURSOR_WINDOW_BEGIN"),
                        "PRECURSOR_WINDOW_END": current.get("PRECURSOR_WINDOW_END"),
                        "PRECURSOR_MZ": current.get("PRECURSOR_MZ"),
                        "PRECURSOR_CHARGE": current.get("PRECURSOR_CHARGE"),
                        "PRECURSOR_MASS": current.get("PRECURSOR_MASS"),
                        "PRECURSOR_INTENSITY": current.get("PRECURSOR_INTENSITY"),
                        "PRECURSOR_FEATURE_ID": current.get("PRECURSOR_FEATURE_ID"),
                        "MSALIGN_number_of_fragment_ions": peak_count  
                    })
                    current = None
            elif "=" in line and current is not None:
                key, val = line.split("=", 1)
                current[key] = val
            elif current is not None:
                parts = line.split() 
                if len(parts) >= 4:
                    peak_count += 1

    return records    


def process_msalign_folder(dataset_id: str, msalign_filename: str, output_filename: str):
    print(f"Extracting metadata from: {msalign_filename} ({dataset_id})")
        
        
    result = msalign_meta_extract(dataset_id, msalign_filename)
        
    df = pd.DataFrame.from_records(result)
    df.to_csv(output_filename, sep='\t', index=False)
    print(f"\nMetadata extracted for {len(df)} MS2 spectra from {msalign_filename}.")
    print(f"Saved to: {output_filename}")



if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <dataset_id> <input_msalign_filename> <output_tsv_filename>")
        sys.exit(1)

    dataset_id = sys.argv[1]
    input_filename = sys.argv[2]
    output_tsv_filename = sys.argv[3]
    process_msalign_folder(dataset_id, input_filename, output_tsv_filename)
