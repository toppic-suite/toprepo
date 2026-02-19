#import pandas as pd
import numpy as np
import argparse
import time
from process.msalign import msalign_reader
from process.msalign import msalign_writer


def merge_msalign_prsm(input_msalign_file, input_tsv_file, output_msalign_file, input_option="raw", format="msalign"):
    is_msalign = (format.lower() == "msalign")
    input_f = open(input_tsv_file, "r", buffering=1024*1024*1024)  # 1G buffer
    header = input_f.readline().strip().split("\t")
    tsv_array = []
    count = 0
    start_time = time.time()
    for line in input_f:
        tsv_array.append(line)
        count += 1
        if count % 100000 == 0:
            elapsed_time = time.time() - start_time
            print(f"\rLoaded {count} rows. Elapsed time: {elapsed_time:.2f} seconds.", end="", flush=True)
    input_f.close()
    print(f"\nLoaded {len(tsv_array)} rows from {input_tsv_file}")

    # Build map using PROJECT_ID, MSALIGN_file_name, MZML_MS2_scan as index
    project_id_idx = header.index("DATASET_id")
    msalign_file_name_idx = header.index("MSALIGN_file_name")
    mzml_ms2_scan_idx = header.index("MZML_ms2_scan")
    tsv_dict = {}
    count = 0
    start_time = time.time()
    for idx, row in enumerate(tsv_array):
        field_values = row.strip().split("\t")
        key = (field_values[project_id_idx], field_values[msalign_file_name_idx], field_values[mzml_ms2_scan_idx])
        tsv_dict[key] = idx
        count += 1
        if count % 100000 == 0:
            elapsed_time = time.time() - start_time
            print(f"\rProcessed {count} rows. Elapsed time: {elapsed_time:.2f} seconds.", end="", flush=True)
    print(f"\nBuilt index map with {len(tsv_dict)} entries.")
    ms_reader = msalign_reader.MsalignReader(input_msalign_file)
    ms_writer = msalign_writer.MsalignWriter(output_msalign_file)
    count = 0 
    output_count = 0
    proteoform_idx = header.index("TOPPIC_proteoform")
    database_seq_idx = header.index("TOPPIC_database_sequence")
    start_idx = header.index("TOPPIC_first_residue_position")
    fixed_mod_idx = header.index("TOPPIC_fixed_modifications")
    unexpected_mod_idx = header.index("TOPPIC_unexpected_modifications")
    collision_energy_idx = header.index("MZML_collision_energy") 
    e_value_idx = header.index("TOPPIC_e-value")
    instrument_idx = header.index("MZML_instrument")
    protein_accession_idx = header.index("TOPPIC_protein_accession")
    for spectrum in ms_reader.readmsalign_iter():
        key = (spectrum["meta"].get("DATASET_ID", ""),
               spectrum["meta"].get("MSALIGN_FILE_NAME", ""),
               spectrum["meta"].get("MS2_SCAN", ""))
        if key in tsv_dict:
            tsv_idx = tsv_dict[key]
            fields = tsv_array[tsv_idx].split("\t")
            field_length = len(fields)
            collision_energy = fields[collision_energy_idx]
            #print (f"\nAnnotating spectrum with key {key}. Field length: {field_length}.")  # Debug print 
            proteoform  = fields[proteoform_idx]
            database_seq = fields[database_seq_idx]
            first_residue_position = fields[start_idx]
            if first_residue_position != "":
                first_residue_position = int(first_residue_position)
            fixed_mod = fields[fixed_mod_idx]
            unexpected_mod = fields[unexpected_mod_idx]
            protein_accesion = fields[protein_accession_idx]
            e_value = fields[e_value_idx]
            if (input_option == "raw"):
                spectrum["meta_lines"].append(f"INSTRUMENT={fields[instrument_idx]}")
                spectrum["meta_lines"].append(f"COLLISION_ENERGY={collision_energy}")
                spectrum["meta_lines"].append(f"PROTEIN_ACCESSION={protein_accesion}")
                spectrum["meta_lines"].append(f"DATABASE_SEQUENCE={database_seq}")
                spectrum["meta_lines"].append(f"FIRST_RESIDUE_POSITION={first_residue_position}")
                spectrum["meta_lines"].append(f"PROTEOFORM={proteoform}")
                spectrum["meta_lines"].append(f"FIXED_MODIFICATIONS={fixed_mod}")
                spectrum["meta_lines"].append(f"UNEXPECTED_MODIFICATIONS={unexpected_mod}")
                spectrum["meta_lines"].append(f"E_VALUE={e_value}")
            if is_msalign:
                ms_writer.write(spectrum)
            else:
                ms_writer.write_mz_intensity(spectrum)   
            output_count += 1
        else:
            print(f"\nWarning: No matching entry found in TSV for spectrum with key {key}. Skipping annotation.")
        count += 1
        if count % 1000 == 0:
            print(f"\rProcessed {count} spectra. Filtered {output_count} spectra.", end="", flush=True)
    print(f"\nFinished processing {count} spectra. Wrote {output_count} spectra to {output_msalign_file}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Annotate MS2 spectra using provided tsv and msalign files.")
    parser.add_argument(
        "--tsv", required=True, type=str, help="Input tsv filename")
    parser.add_argument(
        "--msalign", required=True, type=str, help="Input msalign filename")
    parser.add_argument(
        "--out", type=str, default=None,
        help="Output annotated msalign filename (default: ms2_spectra_annot.msalign)")

    
    args = parser.parse_args()
    output_filename = args.out or "ms2_spectra_annot.msalign"
    merge_msalign_prsm(args.msalign, args.tsv, output_filename)    