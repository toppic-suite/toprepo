import argparse
import msalign_reader
import msalign_writer


def split_msalign_by_instrument(input_msalign_file):
    base_name = input_msalign_file.rsplit(".", 1)[0]
    ms_reader = msalign_reader.MsalignReader(input_msalign_file)

    # Dictionary to hold writers for each instrument type
    instrument_writers = {}
    instrument_counts = {}

    count = 0

    for spectrum in ms_reader.readmsalign_iter():
        instrument = spectrum["meta"].get("INSTRUMENT", "UNKNOWN")

        # Create a writer for this instrument type if it doesn't exist
        if instrument not in instrument_writers:
            instrument_str = instrument.replace(" ", "_")
            output_file = f"{base_name}_{instrument_str}.msalign"
            instrument_writers[instrument] = msalign_writer.MsalignWriter(output_file)
            instrument_counts[instrument] = 0
            print(f"Created output file: {output_file}")

        # Write spectrum to the appropriate instrument file
        instrument_writers[instrument].write(spectrum)
        instrument_counts[instrument] += 1

        count += 1
        if count % 1000 == 0:
            print(f"\rProcessed {count} spectra", end="", flush=True)

    # Close all writers
    for writer in instrument_writers.values():
        writer.close()

    print(f"\n\nFinished processing {count} spectra.")
    print(f"Split into {len(instrument_writers)} instrument type(s):")
    for instrument, inst_count in instrument_counts.items():
        print(f"  {instrument}: {inst_count} spectra")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Split msalign file by instrument type into separate files.")
    parser.add_argument(
        "--msalign", required=True, type=str, help="Input msalign filename")

    args = parser.parse_args()
    split_msalign_by_instrument(args.msalign)    
