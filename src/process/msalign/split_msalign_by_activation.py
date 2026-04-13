import argparse
import msalign_reader
import msalign_writer


def split_msalign_by_activation(input_msalign_file):
    base_name = input_msalign_file.rsplit(".", 1)[0]
    ms_reader = msalign_reader.MsalignReader(input_msalign_file)

    # Dictionary to hold writers for each activation type
    activation_writers = {}
    activation_counts = {}

    count = 0

    for spectrum in ms_reader.readmsalign_iter():
        activation = spectrum["meta"].get("ACTIVATION", "UNKNOWN")

        # Create a writer for this activation type if it doesn't exist
        if activation not in activation_writers:
            output_file = f"{base_name}_{activation}.msalign"
            activation_writers[activation] = msalign_writer.MsalignWriter(output_file)
            activation_counts[activation] = 0
            print(f"Created output file: {output_file}")

        # Write spectrum to the appropriate activation file
        activation_writers[activation].write(spectrum)
        activation_counts[activation] += 1

        count += 1
        if count % 1000 == 0:
            print(f"\rProcessed {count} spectra", end="", flush=True)

    # Close all writers
    for writer in activation_writers.values():
        writer.close()

    print(f"\n\nFinished processing {count} spectra.")
    print(f"Split into {len(activation_writers)} activation type(s):")
    for activation, inst_count in activation_counts.items():
        print(f"  {activation}: {inst_count} spectra")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Split msalign file by activation type into separate files.")
    parser.add_argument(
        "--msalign", required=True, type=str, help="Input msalign filename")

    args = parser.parse_args()
    split_msalign_by_activation(args.msalign)    
