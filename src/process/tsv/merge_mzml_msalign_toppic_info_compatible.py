"""
Merge mzML, msalign, feature, and TopPIC information.

Compatible with TopPIC output formats from versions 1.7.11 and 1.8.1.
The TopPIC format is detected automatically.
"""

import os
import pandas as pd
import sys
import merge_msalign_feature_info as mf
import merge_mzml_msalign_info as mm


def toppic_file_preprocess(top_filename):

    top_df = pd.read_csv(
        top_filename,
        sep="\t",
        low_memory=False,
        dtype=str
    )

    top_df = top_df.drop(
        columns=[
            'Spectrum ID',
            'Charge',
            'Precursor mass',
            'Fragmentation',
            'Feature ID',
            'Retention time',
            '#peaks',
            'Feature intensity',
            'Feature score',
            'Feature apex time'
        ],
        errors='ignore'
    )

    # basename of data file
    top_df["Data file name"] = (
        top_df["Data file name"]
        .apply(os.path.basename)
    )

    # TopPIC 1.8.1 already has these columns
    if (
        "Previous amino acid" in top_df.columns
        and "Next amino acid" in top_df.columns
    ):
        print("Detected newer TopPIC TSV format")

    else:
        print("Detected older TopPIC TSV format")

        col_idx = top_df.columns.get_loc("Proteoform")

        if top_df.shape[0] > 0: 
            parts = top_df["Proteoform"].str.split(".", n=1, expand=True)
            prev_residue = parts[0]
            rest = parts[1].str.rsplit(".", n=1, expand=True)
            proteoform_seq = rest[0]
            next_residue = rest[1]
            top_df["Proteoform"] = proteoform_seq
            top_df.insert(col_idx, "Previous amino acid", prev_residue)
            top_df.insert(col_idx + 2, "Next amino acid", next_residue)
        else:
            top_df.insert(col_idx, "Previous amino acid", "")
            top_df.insert(col_idx + 2, "Next amino acid", "")

    rename_map = {
        'Data file name':
            "MSALIGN file name",
        'Scan(s)':
            "MZML MS2 scan",
        'Prsm ID':
            "TOPPIC PrSM ID",
        'Adjusted precursor mass':
            "TOPPIC adjusted precursor mass",
        'Proteoform ID':
            "TOPPIC proteoform ID",
        'Proteoform intensity':
            "TOPPIC proteoform intensity",
        # -------- v1.7.11 + v1.8.1 --------
        '#Protein hits':
            "TOPPIC number of protein hits",

        '# Protein hits':
            "TOPPIC number of protein hits",

        'Protein accession':
            "TOPPIC protein accession",
        'Protein description':
            "TOPPIC protein description",
        'First residue':
            "TOPPIC first residue position",
        'Last residue':
            "TOPPIC last residue position",
        'Special amino acids':
            "TOPPIC special amino acids",
        'Database protein sequence':
            "TOPPIC database protein sequence",
        'Proteoform mass':
            "TOPPIC proteoform mass",
        'Protein N-terminal form':
            "TOPPIC protein N-terminal form",
        'Fixed PTMs':
            "TOPPIC fixed PTMs",

        # unexpected modifications
        '#unexpected modifications':
            "TOPPIC number of unexpected modifications",

        '# Unexpected modifications':
            "TOPPIC number of unexpected modifications",

        'unexpected modifications':
            "TOPPIC unexpected modifications",

        'Unexpected modifications':
            "TOPPIC unexpected modifications",

        # variable PTMs
        '#variable PTMs':
            "TOPPIC number of variable PTMs",

        '# Variable PTMs':
            "TOPPIC number of variable PTMs",

        'variable PTMs':
            "TOPPIC variable PTMs",

        'Variable PTMs':
            "TOPPIC variable PTMs",

        'MIScore':
            "TOPPIC MIScore",

        # matched peaks/masses
        '#matched peaks':
            "TOPPIC number of matched experimental fragment ions",

        '# Matched masses':
            "TOPPIC number of matched experimental fragment ions",

        '#matched fragment ions':
            "TOPPIC number of matched theoretical fragment masses",

        '# Matched fragments':
            "TOPPIC number of matched theoretical fragment masses",

        'E-value':
            "TOPPIC E-value",

        'Spectrum-level Q-value':
            "TOPPIC spectrum-level Q-value",

        'Proteoform-level Q-value':
            "TOPPIC proteoform-level Q-value",

        'Proteoform':
            "TOPPIC proteoform",

        'Previous amino acid':
            "TOPPIC previous residue",

        'Next amino acid':
            "TOPPIC next residue"
    }

    top_df = top_df.rename(columns=rename_map)
    return top_df


def get_mzml_info_filename(msalign_meta_filename):
    df = pd.read_csv(msalign_meta_filename, sep="\t")
    mzml_filename = df["FILE_NAME"].iloc[0]
    dataset_id = str(df["DATASET_ID"].iloc[0])

    mzml_info_filename = mzml_filename.replace(".mzML", "_mzml_info.tsv")
    # same directory as msalign info file
    folder = os.path.dirname(msalign_meta_filename)    
    # without dataset ID prefix
    filename1 = os.path.join(folder, mzml_info_filename)
    # with dataset ID prefix
    filename2 = os.path.join(folder, f"{dataset_id}_{mzml_info_filename}")
    # check both
    if os.path.isfile(filename1):
        return filename1

    elif os.path.isfile(filename2):
        return filename2

    else:
        raise FileNotFoundError(
            f"Cannot find mzML info file:\n"
            f"{filename1}\n"
            f"or\n"
            f"{filename2}"
        )


def header_str_gen():
        header_str_basic = (
            "DATASET_id\tMZML_file_name\tMZML_instrument\tMZML_ms1_scan\t"
            "MZML_ms1_scan_window_lower_limit\tMZML_ms1_scan_window_upper_limit\t"
            "MZML_ms1_retention_time\tMZML_ms1_total_ion_current\t"
            "MZML_ms1_mass_resolving_power\tMZML_ms1_ion_injection_time\t"
            "MZML_ms1_lowest_observed_mz\tMZML_ms1_highest_observed_mz\t"
            "MZML_ms2_scan\tMZML_ms2_scan_window_lower_limit\t"
            "MZML_ms2_scan_window_upper_limit\tMZML_ms2_retention_time\t"
            "MZML_ms2_total_ion_current\tMZML_ms2_mass_resolving_power\t"
            "MZML_ms2_ion_injection_time\tMZML_ms2_lowest_observed_mz\t"
            "MZML_ms2_highest_observed_mz\tMZML_isolation_window_target_mz\t"
            "MZML_isolation_window_lower_offset\tMZML_isolation_window_upper_offset\t"
            "MZML_selected_ion_mz\tMZML_selected_ion_peak_intensity\t"
            "MZML_selected_ion_charge\tMZML_activation\tMZML_collision_energy\t"
            "MSALIGN_file_name\tMSALIGN_ms1_id\tMSALIGN_ms2_id\t"
            "MSALIGN_precursor_charge\tMSALIGN_precursor_monoisotopic_mass\t"
            "MSALIGN_precursor_intensity\tMSALIGN_feature_id\t"
            "MSALIGN_feature_intensity\tMSALIGN_feature_score\t"
            "MSALIGN_feature_apex_time\tMSALIGN_number_of_fragment_ions\t"
            "TOPPIC_prsm_id\tTOPPIC_adjusted_precursor_mass\tTOPPIC_proteoform_id\t"
            "TOPPIC_proteoform_intensity\tTOPPIC_number_of_protein_hits\t"
            "TOPPIC_protein_accession\tTOPPIC_protein_description\t"
            "TOPPIC_first_residue_position\tTOPPIC_last_residue_position\t"
            "TOPPIC_special_amino_acids\tTOPPIC_database_sequence\t"
            "TOPPIC_proteoform_mass\tTOPPIC_protein_n-terminal_form\t"
            "TOPPIC_fixed_modifications\tTOPPIC_number_of_unexpected_modifications\t"
            "TOPPIC_unexpected_modifications\tTOPPIC_number_of_variable_modifications\t"
            "TOPPIC_variable_modifications\tTOPPIC_miscore\t"
            "TOPPIC_number_of_matched_experimental_fragment_ions\t"
            "TOPPIC_number_of_matched_theoretical_fragment_masses\t"
            "TOPPIC_e-value\tTOPPIC_spectrum-level_q_value\t"
            "TOPPIC_proteoform-level_q_value\tTOPPIC_proteoform\t"
            "TOPPIC_previous_residue\tTOPPIC_next_residue"
        )
        header_str_with_project = (
            "PROJECT_id\tSUBDATASET_id\t" + header_str_basic
        )
        return header_str_with_project, header_str_basic 

def rename_cols_orders(df):     
    # Reorder 
    cols = [
            "DATASET ID", "MZML file name", "MZML instrument", "MZML MS1 scan", 
            "MZML MS1 scan window lower limit", "MZML MS1 scan window upper limit", "MZML MS1 retention time",
            "MZML MS1 total ion current", "MZML MS1 mass resolving power", "MZML MS1 ion injection time",
            "MZML MS1 lowest observed mz", "MZML MS1 highest observed mz",
            "MZML MS2 scan", "MZML MS2 scan window lower limit", "MZML MS2 scan window upper limit",
            "MZML MS2 retention time", "MZML MS2 total ion current",
            "MZML MS2 mass resolving power", "MZML MS2 ion injection time", "MZML MS2 lowest observed mz", "MZML MS2 highest observed mz",
            "MZML isolation window target mz", "MZML isolation window lower offset", "MZML isolation window upper offset",
            "MZML selected ion mz", "MZML selected ion peak intensity", "MZML selected ion charge",
            "MZML activation", "MZML collision energy", "MSALIGN file name", "MSALIGN MS1 ID", "MSALIGN MS2 ID", 
            "MSALIGN precursor charge", "MSALIGN precursor monoisotopic mass",
            "MSALIGN precursor intensity", 
            "MSALIGN feature ID", "MSALIGN feature intensity", "MSALIGN feature score", "MSALIGN feature apex time",
            "MSALIGN number of fragment ions", "TOPPIC PrSM ID", "TOPPIC adjusted precursor mass", "TOPPIC proteoform ID",
            "TOPPIC proteoform intensity", "TOPPIC number of protein hits", 
            "TOPPIC protein accession", "TOPPIC protein description",
            "TOPPIC first residue position", "TOPPIC last residue position",
            "TOPPIC special amino acids","TOPPIC database protein sequence", 
            "TOPPIC proteoform mass","TOPPIC protein N-terminal form", 
            "TOPPIC fixed PTMs", "TOPPIC number of unexpected modifications",
            "TOPPIC unexpected modifications", "TOPPIC number of variable PTMs", 
            "TOPPIC variable PTMs","TOPPIC MIScore", 
            "TOPPIC number of matched experimental fragment ions",
            "TOPPIC number of matched theoretical fragment masses", 
            "TOPPIC E-value","TOPPIC spectrum-level Q-value","TOPPIC proteoform-level Q-value",
            "TOPPIC proteoform", "TOPPIC previous residue", "TOPPIC next residue"
            ]
  
    df = df[cols]
    return df


def rename_cols_with_header_str(df, header_str):
    """Rename DataFrame columns using tab-separated names from header_str.

    The number of tab-separated names in header_str must match the number of
    columns in df (i.e. the DataFrame should already be reordered by
    rename_cols_orders before calling this function).
    """
    new_cols = header_str.split('\t')
    old_cols = list(df.columns)
    if len(new_cols) != len(old_cols):
        raise ValueError(
            f"header_str has {len(new_cols)} columns but DataFrame has {len(old_cols)}"
        )
    return df.rename(columns=dict(zip(old_cols, new_cols)))


def info_merge(msalign_meta_filename, feature_meta_filename, top_filename, output_file, file_info_filename=None):
    """
    msalign_meta_filename: extracted meta info from msalign file
    feature_meta_filename: extracted meta info from feature file 
    top_filename: toppic output tsv file, 
    output_file: file name of the output in tsv format 
    """
    header_str_with_project, header_str_basic = header_str_gen()

    top_df = toppic_file_preprocess(top_filename)
    # merge msalign + feature file
    ms2_feature_df = mf.feature_merge(msalign_meta_filename, feature_meta_filename, out_filename=None, wfile=False)

    # find corresponding mzML info file
    mzml_meta_filename = get_mzml_info_filename(msalign_meta_filename)
    # merge mzml + msalign + feature file
    meta_df = mm.meta_merge(mzml_meta_filename, ms2_feature_df, output_file=None, wfile=False)
    meta_df = meta_df.drop(columns=['title'], errors='ignore')

    meta_df["MSALIGN file name"] = meta_df["MSALIGN file name"].astype(str)
    meta_df["MZML MS2 scan"] = meta_df["MZML MS2 scan"].astype(str)
    meta_df["DATASET ID"] = meta_df["DATASET ID"].astype(str)

    # merge mzml + msalign + feature + toppic 
    # merge: keep all metadata rows, even if no match in TopPIC
    top_df_merged = meta_df.merge(
        top_df,
        on = ['DATASET ID', 'MSALIGN file name', 'MZML MS2 scan'],
        how='left',
        suffixes=('', '_new')
    )
    #print(top_df_merged.columns)
    # rename column names and change orders
    top_df_merged = rename_cols_orders(top_df_merged)

    # inner join with file info TSV if provided
    if file_info_filename is not None:
        file_info_df = pd.read_csv(file_info_filename, sep="\t", dtype=str)
        # Make either form work
        file_info_df = file_info_df.rename(
            columns={
                "DATASET id": "DATASET ID",
                "DATASET_id": "DATASET ID"
            }
        )
        required_cols = ["DATASET ID", "MSALIGN file name", "PROJECT id", "SUBDATASET id"]
        missing_cols = [
            c for c in required_cols
            if c not in file_info_df.columns
        ]

        if missing_cols:
            raise ValueError(
                f"Missing required columns in file info TSV: {missing_cols}"
            )

        file_info_join = file_info_df[required_cols]
        top_df_merged = top_df_merged.merge(
            file_info_join,
            on=["DATASET ID", "MSALIGN file name"],
            how="inner"
        )

        cols = top_df_merged.columns.tolist()
        new_cols = (
            ["PROJECT id", "SUBDATASET id"]
            + [
                c for c in cols
                if c not in ["PROJECT id", "SUBDATASET id"]
            ]
        )

        top_df_merged = top_df_merged[new_cols]

        header_str = header_str_with_project
    else:
        header_str = header_str_basic

    if header_str is not None:
        top_df_merged = rename_cols_with_header_str(top_df_merged, header_str)
    # Save merged file
    id_cols = [
        "TOPPIC_first_residue_position",
        "TOPPIC_last_residue_position"
    ]
    for c in id_cols:
        if c in top_df_merged.columns:
            top_df_merged[c] = top_df_merged[c].fillna("").astype(str)

    top_df_merged.to_csv(output_file, sep='\t', index=False)
    print(f"Merged file saved to: {output_file}")
    print(f"Total rows: {len(top_df_merged)}, matched: {top_df_merged['TOPPIC_proteoform_id'].notna().sum()}")
    


if __name__ == "__main__":
    if len(sys.argv) not in [5, 6]:
        print("Usage: python script.py <msalign_info_filename> <feature_info_filename> <toppic_info_filename> [file_info_filename] <output_tsv_filename>")
        sys.exit(1)

    msalign_info_filename = sys.argv[1]
    feature_info_filename = sys.argv[2]
    toppic_info_filename = sys.argv[3]

    # without file_info
    if len(sys.argv) == 5:
        file_info_filename = None
        out_filename = sys.argv[4]

    # with file_info
    else:
        file_info_filename = sys.argv[4]
        out_filename = sys.argv[5]

    info_merge(msalign_info_filename, feature_info_filename, toppic_info_filename, out_filename, file_info_filename)
    
