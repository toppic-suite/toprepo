#import pandas as pd
import numpy as np
import argparse
import msalign_reader
import msalign_writer


def filter_msalign(input_msalign_file, output_msalign_file): 
    ms_reader = msalign_reader.MsalignReader(input_msalign_file)
    ms_writer = msalign_writer.MsalignWriter(output_msalign_file)
    count = 0 
    output_count = 0
    for spectrum in ms_reader.readmsalign_iter():
        seq = spectrum["meta"].get("PROTEOFORM", "")
        if (seq.strip()) != "" and ("[" not in seq):
            ms_writer.write(spectrum)
            output_count += 1
        count += 1
        if count % 1000 == 0:
            print(f"\rProcessed {count} spectra. Filtered {output_count} spectra.", end="", flush=True)
    print(f"\nFinished processing {count} spectra. Wrote {output_count} spectra to {output_msalign_file}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Annotate MS2 spectra using provided tsv and msalign files.")
    parser.add_argument(
        "--msalign", required=True, type=str, help="Input msalign filename")
    parser.add_argument(
        "--out", type=str, default=None,
        help="Output msalign filename (default: ms2_spectra_annot.msalign)")

    args = parser.parse_args()
    output_filename = args.out or "ms2_spectra_annot.msalign"
    filter_msalign(args.msalign, output_filename)    
