"""
remove duplicate proteoform
"""

import pandas as pd
import sys


def calculate_overlap_percentage(r1_start, r1_end, r2_start, r2_end):
    try:
        r1_start, r1_end = int(r1_start), int(r1_end)
        r2_start, r2_end = int(r2_start), int(r2_end)

        overlap_start = max(r1_start, r2_start)
        overlap_end   = min(r1_end, r2_end)

        if overlap_start > overlap_end:
            return 0.0

        overlap_len = overlap_end - overlap_start + 1
        len1 = r1_end - r1_start + 1
        len2 = r2_end - r2_start + 1
        shorter_len = min(len1, len2)
        return overlap_len / shorter_len if shorter_len > 0 else 0
    except:
        return 0.0


def is_duplicate(a, b):
    """Return True if b is duplicate of a and should be removed."""
    
    if a["MSALIGN_file_name"] == b["MSALIGN_file_name"]:
        return False

    if a["TOPPIC_protein_accession"] != b["TOPPIC_protein_accession"]:
        return False

    if abs(float(a['MSALIGN_precursor_mass']) -
           float(b['MSALIGN_precursor_mass'])) > 2.2:
        return False
    
    ov = calculate_overlap_percentage(
        a["TOPPIC_first_residue_position"], a["TOPPIC_last_residue_position"],
        b["TOPPIC_first_residue_position"], b["TOPPIC_last_residue_position"]
    )
    if ov < 0.70:
        return False

    return True


def deduplicate_group(df):
    df = df.copy()
    df["TOPPIC_e-value"] = df["TOPPIC_e-value"].astype(float)

    files = []
    for f, sub in df.groupby("MSALIGN_file_name"):
        sub = sub.sort_values("TOPPIC_e-value")  
        files.append(sub)

    survivors = []

    for file_df in files:
        for _, row in file_df.iterrows():
            replace_idx = None
            is_dup = False

            for i, s in enumerate(survivors):
                if is_duplicate(s, row):
                    is_dup = True

                    if row["TOPPIC_e-value"] < s["TOPPIC_e-value"]:
                        replace_idx = i
                    break

            if replace_idx is not None:
                survivors[replace_idx] = row

            elif not is_dup:
                survivors.append(row)

    return pd.DataFrame(survivors)


def main(input_file, output_file):
    df = pd.read_csv(
        input_file,
        sep="\t",
        low_memory=False
    )
    print(f"{len(df)}")
    df['MSALIGN_precursor_mass'] = df['MSALIGN_precursor_monoisotopic_mass'].apply(lambda x: x.split(':')[0])
    groups = []
    for acc, group in df.groupby(
        ["TOPPIC_protein_accession"]
    ):
        groups.append(deduplicate_group(group))
    
    deduped_chunk = pd.concat(groups, ignore_index=True)
    print(f"Producing {len(deduped_chunk)} deduplicated rows")
    
    deduped_chunk.drop(columns=["MSALIGN_precursor_mass"], axis=1, inplace=True)
    num_ptms = deduped_chunk['TOPPIC_proteoform_id'].notna().sum()
    
    deduped_chunk.to_csv(
        output_file,
        sep="\t",
        index=False)


    print(f"Total spectra (rows): {len(df)}, total remove duplicates: {len(deduped_chunk)}")
    print(f"Rows with proteoform ID: {num_ptms}")
    print("Done.")
    

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    main(input_file, output_file)

