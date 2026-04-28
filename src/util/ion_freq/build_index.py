import argparse
import pickle


def build_msalign_index(input_file, output_index):

    activation_lists = {
        "cid": [],
        "hcd": [],
        "etd": [],
        "ethcd": []
    }

    total_blocks = 0
    kept_blocks = 0

    with open(input_file, "r") as f:

        inside_block = False
        block_start = None
        activation = None
        proteoform = ""

        while True:
            pos = f.tell()
            line = f.readline()

            if not line:
                break

            line = line.rstrip()

            if line.startswith("BEGIN IONS"):
                inside_block = True
                block_start = pos
                activation = None
                proteoform = ""

            elif inside_block and line.startswith("ACTIVATION="):
                activation = line.split("=", 1)[1].strip().lower()

            elif inside_block and line.startswith("PROTEOFORM="):
                proteoform = line.split("=", 1)[1].strip()

            elif line.startswith("END IONS") and inside_block:

                inside_block = False
                total_blocks += 1

                # keep only spectra with activation + proteoform
                if activation in activation_lists and proteoform != "":
                    activation_lists[activation].append(block_start)
                    kept_blocks += 1

                if total_blocks % 10000 == 0:
                    print(f"Processed {total_blocks} spectra...")

    # save index
    with open(output_index, "wb") as f:
        pickle.dump(activation_lists, f)

    print("\nIndex built successfully")
    print(f"Total spectra scanned: {total_blocks}")
    print(f"Valid spectra indexed: {kept_blocks}")

    for k, v in activation_lists.items():
        print(f"{k.upper()} spectra: {len(v)}")

    print(f"\nIndex saved to: {output_index}")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Build index for msalign file")
    parser.add_argument("--input", required=True, help="Input msalign file")
    parser.add_argument("--output", required=True, help="Output index file (.pkl)")

    args = parser.parse_args()

    build_msalign_index(args.input, args.output)
    
    