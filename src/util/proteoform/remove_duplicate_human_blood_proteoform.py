"""
remove duplicate proteoform
"""

import pandas as pd
import sys


HUMAN_BLOOD_DATASETS = list(set([
    'PXD026123','PXD026124','PXD026126','PXD026127','PXD026128','PXD026129',
    'PXD026130','PXD026131','PXD026132','PXD026133','PXD026134','PXD026135',
    'PXD026136','PXD026137','PXD026138','PXD026139','PXD026140','PXD026141',
    'PXD026142','PXD026143','PXD026144','PXD026145','PXD026146','PXD026147',
    'PXD026148','PXD026149','PXD026150','PXD026151','PXD026152','PXD026153',
    'PXD026154','PXD026155','PXD026156','PXD026157','PXD026158','PXD026159',
    'PXD026160','PXD026161','PXD026162','PXD026163','PXD026164','PXD026165',
    'PXD026166','PXD026167','PXD026168','PXD026169','PXD026170','PXD026171',
    'PXD026172','PXD026173','PXD026174','PXD026175','PXD026176','PXD026177',
    'PXD026178'
]))


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


def process_human_blood(df):
    human_blood_df = df[df['DATASET_id'].isin(HUMAN_BLOOD_DATASETS)].copy()
    rest_df = df[~df['DATASET_id'].isin(HUMAN_BLOOD_DATASETS)].copy()

    return human_blood_df, rest_df


def main(input_file, output_file):
    df = pd.read_csv(
        input_file,
        sep="\t",
        low_memory=False
    )
    print(f"{len(df)}")  
    human_df, rest_df = process_human_blood(df)
    print(f"Human blood rows: {len(human_df)}")
    print(f"Rest rows: {len(rest_df)}")

    # --- Deduplicate human blood ---
    human_df['MSALIGN_precursor_mass'] = human_df['MSALIGN_precursor_monoisotopic_mass'].apply(lambda x: x.split(':')[0])
    groups = []
    for acc, group in human_df.groupby(
        ["TOPPIC_protein_accession"]
    ):
        groups.append(deduplicate_group(group))
    
    dedup_human_df = pd.concat(groups, ignore_index=True)
    print(f"Producing {len(dedup_human_df)} deduplicated rows")
    
    dedup_human_df.drop(columns=["MSALIGN_precursor_mass"], axis=1, inplace=True)
    final_df = pd.concat([rest_df, dedup_human_df], ignore_index=True)
    
    num_ptms = final_df['TOPPIC_proteoform_id'].notna().sum()
    final_df.to_csv(
        output_file,
        sep="\t",
        index=False)


    print(f"Total spectra (rows): {len(final_df)}")
    print(f"Rows with proteoform ID: {num_ptms}")
    print("Done.")
    

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    main(input_file, output_file)

