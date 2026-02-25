import os
import pandas as pd
import sys


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


def info_merge(top_filename, mzml_meta_filename, output_file, header_str=None):
    """
    top_filename: toppic output tsv file, 
    mzml_meta_filename: extracted metadata from mzml file
    output_file: file name of the output in tsv format 
    """
    top_df = pd.read_csv(top_filename, sep='\t',low_memory=False, dtype=str)
    top_df = top_df.drop(columns=['Spectrum ID', 'Charge', 'Precursor mass', 'Fragmentation', 'Feature ID', 'Retention time', '#peaks','Feature intensity','Feature score', 'Feature apex time'], errors='ignore')
    top_df["Data file name"] = top_df["Data file name"].apply(os.path.basename)

    # Split Proteoform into three parts on first and last "."
    parts = top_df["Proteoform"].str.split(".", n=1, expand=True)
    prev_residue = parts[0]
    rest = parts[1].str.rsplit(".", n=1, expand=True)
    proteoform_seq = rest[0]
    next_residue = rest[1]

    # Insert new columns, replacing original Proteoform
    col_idx = top_df.columns.get_loc("Proteoform")
    top_df.insert(col_idx, "Previous residue", prev_residue)
    top_df["Proteoform"] = proteoform_seq
    top_df.insert(col_idx + 2, "Next residue", next_residue)
    #print(top_df.columns)  
    top_df = top_df.rename(columns={
        'Data file name': "MSALIGN file name",
        'Scan(s)': "MZML MS2 scan", 
        'Prsm ID': "TOPPIC PrSM ID",
        'Adjusted precursor mass': "TOPPIC adjusted precursor mass",
        'Proteoform ID': "TOPPIC proteoform ID", 
        'Proteoform intensity': "TOPPIC proteoform intensity",
        '#Protein hits': "TOPPIC number of protein hits", 
        'Protein accession': "TOPPIC protein accession", 
        'Protein description': "TOPPIC protein description",
        'First residue': "TOPPIC first residue position",
        'Last residue': "TOPPIC last residue position",
        'Special amino acids': "TOPPIC special amino acids",
        'Database protein sequence': "TOPPIC database protein sequence", 
        'Proteoform mass': "TOPPIC proteoform mass",
        'Protein N-terminal form': "TOPPIC protein N-terminal form", 
        'Fixed PTMs': "TOPPIC fixed PTMs", 
        '#unexpected modifications': "TOPPIC number of unexpected modifications",
        'unexpected modifications': "TOPPIC unexpected modifications", 
        '#variable PTMs': "TOPPIC number of variable PTMs", 
        'variable PTMs': "TOPPIC variable PTMs",
        'MIScore': "TOPPIC MIScore", 
        '#matched peaks': "TOPPIC number of matched experimental fragment ions",
        '#matched fragment ions': "TOPPIC number of matched theoretical fragment masses", 
        'E-value': "TOPPIC E-value",
        'Spectrum-level Q-value': "TOPPIC spectrum-level Q-value", 
        'Proteoform-level Q-value': "TOPPIC proteoform-level Q-value",
        'Proteoform': "TOPPIC proteoform", 
        'Previous residue': "TOPPIC previous residue", 
        'Next residue': "TOPPIC next residue",
       })
    
    meta_df = pd.read_csv(mzml_meta_filename, sep='\t', dtype=str)
    meta_df = meta_df.drop(columns=['title'])
    
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
    if header_str is not None:
        top_df_merged = rename_cols_with_header_str(top_df_merged, header_str)
    # Save merged file
    top_df_merged.to_csv(output_file, sep='\t', index=False)
    print(f"Merged file saved to: {output_file}")
    print(f"Total rows: {len(top_df_merged)}, matched: {top_df_merged['TOPPIC_proteoform_id'].notna().sum()}")
    

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <toppic_info_filename> <mzml_msalign_info_filename> <output_tsv_filename>")
    else:
        toppic_info_filename = sys.argv[1]
        mzml_msalign_info_filename = sys.argv[2]
        out_filename = sys.argv[3]
        header_str ="DATASET_id	MZML_file_name	MZML_instrument	MZML_ms1_scan	MZML_ms1_scan_window_lower_limit	MZML_ms1_scan_window_upper_limit	MZML_ms1_retention_time	MZML_ms1_total_ion_current	MZML_ms1_mass_resolving_power	MZML_ms1_ion_injection_time	MZML_ms1_lowest_observed_mz	MZML_ms1_highest_observed_mz	MZML_ms2_scan	MZML_ms2_scan_window_lower_limit	MZML_ms2_scan_window_upper_limit	MZML_ms2_retention_time	MZML_ms2_total_ion_current	MZML_ms2_mass_resolving_power	MZML_ms2_ion_injection_time	MZML_ms2_lowest_observed_mz	MZML_ms2_highest_observed_mz	MZML_isolation_window_target_mz	MZML_isolation_window_lower_offset	MZML_isolation_window_upper_offset	MZML_selected_ion_mz	MZML_selected_ion_peak_intensity	MZML_selected_ion_charge	MZML_activation	MZML_collision_energy	MSALIGN_file_name	MSALIGN_ms1_id	MSALIGN_ms2_id	MSALIGN_precursor_charge	MSALIGN_precursor_monoisotopic_mass	MSALIGN_precursor_intensity	MSALIGN_feature_id	MSALIGN_feature_intensity	MSALIGN_feature_score	MSALIGN_feature_apex_time	MSALIGN_number_of_fragment_ions	TOPPIC_prsm_id	TOPPIC_adjusted_precursor_mass	TOPPIC_proteoform_id	TOPPIC_proteoform_intensity	TOPPIC_number_of_protein_hits	TOPPIC_protein_accession	TOPPIC_protein_description	TOPPIC_first_residue_position	TOPPIC_last_residue_position	TOPPIC_special_amino_acids	TOPPIC_database_sequence	TOPPIC_proteoform_mass	TOPPIC_protein_n-terminal_form	TOPPIC_fixed_modifications	TOPPIC_number_of_unexpected_modifications	TOPPIC_unexpected_modifications	TOPPIC_number_of_variable_modifications	TOPPIC_variable_modifications	TOPPIC_miscore	TOPPIC_number_of_matched_experimental_fragment_ions	TOPPIC_number_of_matched_theoretical_fragment_masses	TOPPIC_e-value	TOPPIC_spectrum-level_q-value	TOPPIC_proteoform-level_q-value	TOPPIC_proteoform	TOPPIC_previous_residue	TOPPIC_next_residue"
        info_merge(toppic_info_filename, mzml_msalign_info_filename, out_filename, header_str)
        
        
        
