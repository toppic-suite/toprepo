import argparse
import pickle
import random
import re

def is_clean_sequence(seq):
    # detect any modification patterns
    if re.search(r"\(.*?\)", seq):
        return False
    if re.search(r"\[.*?\]", seq):
        return False
    if re.search(r"\.", seq):
        return False
    return True


def extract_random_subset(input_file, index_file, output_file, activation, num):

    with open(index_file, "rb") as f:
        index = pickle.load(f)

    activation = activation.lower()

    if activation not in index:
        print("Activation not found in index.")
        return

    positions = index[activation]
    print(f"Total {activation.upper()} spectra available: {len(positions)}")

    valid_positions = []
    valid_count = 0
    with open(input_file, "r") as fin:

        for pos in positions:

            fin.seek(pos)

            proteoform = None
            clean_seq = None
            evalue = None

            while True:
                line = fin.readline()
                if not line:
                    break

                line = line.strip()

                if line.startswith("PROTEOFORM="):
                    proteoform = line.split("=",1)[1].strip()

                elif line.startswith("E_VALUE="):
                    val = line.split("=")[1].strip()
                    evalue = float(val) if val else None

                elif line.startswith("DATABASE_SEQUENCE="):
                    clean_seq = line.split("=")[1].strip().upper()

                if line.startswith("END IONS"):
                    break

            # filtering
            if not proteoform or not is_clean_sequence(proteoform):
                continue

            if clean_seq is None or len(clean_seq) > 200:
                continue

            if evalue is None or evalue > 1e-4:
                continue
            
            valid_count +=1
            valid_positions.append(pos)

    print(f"Total valid spectra: {len(valid_positions)}")
    
    random.seed(42)
    if len(valid_positions) >= num:
        selected_positions = sorted(random.sample(valid_positions, num))
    else:
        need = num - len(valid_positions) 
        valid_set = set(valid_positions) 
        invalid_pool = [p for p in positions if p not in valid_set]
        extra = random.sample(invalid_pool, need)
        selected_positions = sorted(valid_positions + extra)
        

    print(f"Randomly selecting {len(selected_positions)} spectra")

    count = 0

    with open(input_file, "r") as fin, open(output_file, "w") as fout:

        for pos in selected_positions:

            fin.seek(pos)

            while True:
                line = fin.readline()
                if not line:
                    break

                line = line.rstrip()

                if line.startswith("SEQUENCE_COVERAGE="):
                    continue

                if line and line[0].isdigit():
                    parts = line.split()
                    if len(parts) >= 4:
                        fout.write("\t".join(parts[:4]) + "\n")
                else:
                    fout.write(line + "\n")

                if line.startswith("END IONS"):
                    break

            count += 1

    print(f"Extracted {count} spectra to {output_file}")
    
     
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Extract random spectra from msalign using index")

    parser.add_argument("--input", required=True, help="Input msalign file")
    parser.add_argument("--index", required=True, help="Index file (.pkl)")
    parser.add_argument("--output", required=True, help="Output msalign file")
    parser.add_argument("--activation", required=True, help="CID / HCD / ETD / ETHCD")
    parser.add_argument("--num", type=int, default=10000, help="Number of spectra")

    args = parser.parse_args()

    extract_random_subset(args.input, args.index, args.output, args.activation, args.num)
    