from collections import Counter
import argparse
import re


def generalize_ion(ion):
    ion = ion.strip()  # remove hidden spaces/newlines    
    match = re.match(r"^(z_dot|[abcxyz])(\d+)([-+].*)?$", ion)
    if match:
        ion_type = match.group(1)
        suffix = match.group(3) if match.group(3) else ""
        return ion_type + suffix
    return ion


def process_annot_msalign(msalign_filename, out_filename):
    ion_counter = Counter()
    total_ion_num = 0
    with open(msalign_filename, "r") as f:
        inside_block = False
        for line in f:
            line = line.strip()
            if line.startswith("BEGIN IONS"):
                inside_block = True
                continue
            if line.startswith("END IONS"):
                inside_block = False
                continue
            if inside_block and line:
                parts = line.split()
                if len(parts) >= 5:  # check if annotation exists
                    ions = parts[4:]
                    for ion in ions:
                        gen = generalize_ion(ion)
                        if gen:
                            ion_counter[gen] += 1
            if inside_block and line.startswith("DATABASE_SEQUENCE"):
                part2 = line.split("=", 1)[1].strip()
                if part2:
                    total_ion_num +=  (len(part2) -1)
                    

    with open(out_filename, "w") as out_f:
        out_f.write("Ion\tCount\n")
        for ion, count in sorted(ion_counter.items(), key=lambda x: x[1], reverse=True):
            out_f.write(f"{ion}\t{count}\t{(count*100.0/total_ion_num):.2f}%\n")
        
        
        
if __name__ == "__main__":    
    parser = argparse.ArgumentParser(
    description="Count frequency for an annotated msalign file.")
    parser.add_argument(
        "--msalign", required=True, type=str, help="Input msalign filename")
    parser.add_argument(
        "--out", required=True, type=str, default="annot_frequency.tsv",
        help="Output annotated frequnecy tsv filename")
    args = parser.parse_args()
    process_annot_msalign(args.msalign, args.out)
    
