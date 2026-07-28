#import pandas as pd
import numpy as np
import argparse
import msalign_reader
import msalign_writer

SEQ_LENGTH_CUTOFF = 200


def filter_msalign(input_msalign_file, output_msalign_file, evalue_cutoff=1e-4,
                   skip_seq_length_filter=False, long_seq_only=False):
    ms_reader = msalign_reader.MsalignReader(input_msalign_file)
    ms_writer = msalign_writer.MsalignWriter(output_msalign_file)
    count = 0 
    output_count = 0
    seq_length_filter_count = 0
    e_value_filter_count = 0
    seq_coverage_filter_count = 0
    precursor_charge_filter_count = 0
    for spectrum in ms_reader.readmsalign_iter():
        count += 1
        if count % 1000 == 0:
            print(f"\rProcessed {count} spectra. Sequence filter: {seq_length_filter_count}, E-value filter: {e_value_filter_count}, Precursor charge filter: {precursor_charge_filter_count}, Coverage filter: {seq_coverage_filter_count}. Output {output_count} spectra.", end="", flush=True)
        seq = spectrum["meta"].get("DATABASE_SEQUENCE")
        if not skip_seq_length_filter:
            if long_seq_only:
                keep_seq = len(seq) > SEQ_LENGTH_CUTOFF
            else:
                keep_seq = len(seq) <= SEQ_LENGTH_CUTOFF
            if not keep_seq:
                seq_length_filter_count += 1
                continue
        e_value = float(spectrum["meta"].get("E_VALUE", 10000))
        if (e_value > evalue_cutoff):
            e_value_filter_count += 1
            continue
        prec_charge_str = spectrum["meta"].get("PRECURSOR_CHARGE", "0")
        prec_charge_list = prec_charge_str.split(":")
        prec_charge = int(prec_charge_list[0]) if len(prec_charge_list) > 0 else 0
        if prec_charge > 30:
            precursor_charge_filter_count += 1
            continue
        # parse peaks 
        peak_lines = spectrum["peak_lines"]
        ion_coverage = np.zeros((2, len(seq)-1, 30), dtype=np.int8)  # 0: b-ion, 1: y-ion
        for i in range(len(peak_lines)):
            peak_lines[i] = peak_lines[i].strip()
            parts = peak_lines[i].split()
            if (len(parts) == 10):
                ion_str = parts[4][0].lower()[:1]
                if (ion_str == 'b' or ion_str == 'c'):
                    ion_type = 0
                elif (ion_str == 'y' or ion_str == 'z'):
                    ion_type = 1
                else:
                    print(f"Unknown ion type: {parts[4]}")
                    continue
                ion_pos = int(parts[5])
                ion_charge = int(parts[2])
                if ion_charge > prec_charge:
                    continue
                ion_coverage[ion_type, ion_pos-1, ion_charge-1] = 1
        # compute sequence coverage
        cover_count = ion_coverage.sum()
        cover_ratio = cover_count / (2*(len(seq)-1))

        if (cover_ratio < 0.2):
            #print(f"seq: {seq}, cover_count: {cover_count}, cover_ratio: {cover_ratio}")
            seq_coverage_filter_count += 1
            continue
        ms_writer.write(spectrum)
        output_count += 1
    print(f"\nFinished processing {count} spectra. Sequence filter: {seq_length_filter_count}, E-value filter: {e_value_filter_count}, Precursor charge filter: {precursor_charge_filter_count}, Coverage filter: {seq_coverage_filter_count}. Wrote {output_count} spectra to {output_msalign_file}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Annotate MS2 spectra using provided tsv and msalign files.")
    parser.add_argument(
        "--msalign", required=True, type=str, help="Input msalign filename")
    parser.add_argument(
        "--out", type=str, default=None,
        help="Output msalign filename (default: ms2_spectra_annot.msalign)")
    parser.add_argument(
        "--E_value_threshold", type=float, default=1e-4,
        help="E-value filter threshold (default: 1e-4)")
    # By default only sequences with a length <= 200 are kept. The two options below
    # change that, and cannot be combined.
    seq_length_group = parser.add_mutually_exclusive_group()
    seq_length_group.add_argument(
        "--skip_seq_length_filter", action='store_true',
        help=f"Do not filter on sequence length (default: keep length <= {SEQ_LENGTH_CUTOFF})")
    seq_length_group.add_argument(
        "--long_seq_only", action='store_true',
        help=f"Keep only sequences with a length > {SEQ_LENGTH_CUTOFF}")


    args = parser.parse_args()
    output_filename = args.out or "ms2_spectra_annot.msalign"
    filter_msalign(args.msalign, output_filename, evalue_cutoff = args.E_value_threshold,
                   skip_seq_length_filter = args.skip_seq_length_filter,
                   long_seq_only = args.long_seq_only)
