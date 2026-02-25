import sys
# import os
# import re
import pandas as pd


def meta_merge(mzml_meta_filename, msalign_meta_filename, output_file):
    """
    mzml_meta_filename: meta extracted from mzml file in tsv, 
    msalign_meta_filename: meta extracted from msalign file in tsv
    output_file: file name of the output in tsv format 
    """
    meta_df1 = pd.read_csv(mzml_meta_filename, sep='\t',low_memory=False)
    meta_df1 = meta_df1.rename(columns={
        'dataset_id': "DATASET ID",
        'instrument': "MZML instrument",
        'file_name': "MZML file name",
        'ms2_scan_id': "MZML MS2 scan", 
        'ms2_scan_begin': "MZML MS2 scan window lower limit", 
        'ms2_scan_end': "MZML MS2 scan window upper limit", 
        'ms2_retention_time': "MZML MS2 retention time", 
        'pepmass_mz': "MZML selected ion mz",
        'selected_ion_charge': "MZML selected ion charge", 
        'peak_intensity': "MZML selected ion peak intensity", 
        'collision_energy': "MZML collision energy",
        'ms2_total_ion_current': "MZML MS2 total ion current", 
        'ms2_lower_obsevered_mz': "MZML MS2 lowest observed mz",
        'ms2_highest_obsevered_mz': "MZML MS2 highest observed mz", 
        'ms2_injection_time': "MZML MS2 ion injection time", 
        'ms2_resolution': "MZML MS2 mass resolving power",
        'isolation_window_mz': "MZML isolation window target mz", 
        'isolation_window_lower_offset': "MZML isolation window lower offset",
        'isolation_window_upper_offset': "MZML isolation window upper offset", 
        'ms1_scan_id': "MZML MS1 scan", 
        'ms1_scan_begin': "MZML MS1 scan window lower limit",
        'ms1_scan_end': "MZML MS1 scan window upper limit", 
        'ms1_retention_time': "MZML MS1 retention time", 
        'ms1_total_ion_current': "MZML MS1 total ion current",
        'ms1_injection_time': "MZML MS1 ion injection time", 
        'ms1_resolution': "MZML MS1 mass resolving power", 
        'ms1_lower_obsevered_mz': "MZML MS1 lowest observed mz",
        'ms1_highest_obsevered_mz': "MZML MS1 highest observed mz"
        })
    
    meta_df2 = pd.read_csv(msalign_meta_filename, sep='\t')
    meta_df2 = meta_df2.drop(columns=['feature_id','precursor_intensity'], errors='ignore')
    meta_df2 = meta_df2.rename(columns={
        'FILE_NAME': 'MZML file name',
        'DATASET_ID': 'DATASET ID',
        'MS_ONE_ID': 'MSALIGN MS1 ID',
        'MS_ONE_SCAN': 'MZML MS1 scan',
        'MSALIGN_FILE_NAME': 'MSALIGN file name',
        'MS2_SCANS': 'MZML MS2 scan',
        'SPECTRUM_ID': 'MSALIGN MS2 ID', 
        'ACTIVATION': 'MZML activation',
        'PRECURSOR_WINDOW_BEGIN': 'MSALIGN precursor window begin',
        'PRECURSOR_WINDOW_END': 'MSALIGN precursor window end',
        'PRECURSOR_MZ': 'MSALIGN precursor mz',
        'PRECURSOR_CHARGE': 'MSALIGN precursor charge',
        'PRECURSOR_MASS': 'MSALIGN precursor monoisotopic mass',
        'PRECURSOR_INTENSITY': 'MSALIGN precursor intensity',
        'PRECURSOR_FEATURE_ID': 'MSALIGN feature ID',
        'MSALIGN_number_of_fragment_ions': 'MSALIGN number of fragment ions',
        'feature_intensity': 'MSALIGN feature intensity',
        'feature_score': 'MSALIGN feature score',
        'feature_apex_time': 'MSALIGN feature apex time'
    })

    merged_meta_df = meta_df1.merge(
        meta_df2,
        on=['DATASET ID', 'MZML file name', 'MZML MS2 scan'],
        how='left',
        suffixes=('', '_new')
    )

    # Save merged file
    merged_meta_df.to_csv(output_file, sep='\t', index=False)
    print(f"Merged file saved to: {output_file}")
    print(f"Total rows: {len(merged_meta_df)}")
    

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <mzml_info_filename> <msalign_info_filename> <output_tsv_filename>")
    else:
        mzml_info_tsv = sys.argv[1]
        msalign_info_tsv = sys.argv[2]
        out_filename = sys.argv[3]
        meta_merge(mzml_info_tsv, msalign_info_tsv, out_filename)