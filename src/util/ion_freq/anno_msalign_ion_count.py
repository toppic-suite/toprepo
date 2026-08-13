from collections import Counter
import csv
import re
import argparse


def generalize_ion(ion):
    ion = ion.strip() 
    match = re.match(r"^(z_dot|[abcxyz])(\d+)([-+].*)?$", ion)
    if match:
        ion_type = match.group(1)
        suffix = match.group(3) if match.group(3) else ""
        return ion_type + suffix
    return ion


def count_fragment_stats(msalign_file):

    total_exp_peaks = 0
    annotated_peaks = 0
    total_seq_len = 0
    ion_counts = Counter()

    with open(msalign_file) as f:

        for line in f:

            line = line.strip()
            
            if line.startswith("DATABASE_SEQUENCE="):
                clean_seq = line.split("=")[1].strip().upper()
                peptide_seq = len(clean_seq)-1    
                total_seq_len += peptide_seq
                continue

            if not line or not line[0].isdigit():
                continue

            cols = line.split()

            total_exp_peaks += 1

            if len(cols) > 4:

                annotated_peaks += 1

                ion = cols[4] 
                ion_type = generalize_ion(ion)
                ion_counts[ion_type] += 1

    
    print("Total experimental peaks:", total_exp_peaks)
    print("Annotated peaks:", annotated_peaks)
    fragmentation_rate = annotated_peaks / total_exp_peaks
    print(f"Annotated rate: {fragmentation_rate}")
    print("Total theoretical cleavge sites (sum(L-1)):", total_seq_len)

    print("\nIon type counts:")
    for ion, count in sorted(ion_counts.items(), key=lambda x: x[1], reverse=True):
        coverage = count / total_seq_len
        print(f"{ion}: {count}, coverage={coverage:.4f}")

    return ion_counts, total_seq_len, total_exp_peaks, annotated_peaks, fragmentation_rate


def write_ion_counts_tsv(ion_counter, total_seq_len, total_exp_peaks, annotated_peaks, fragmentation_rate, output_file):

    with open(output_file, "w", newline="") as f:
        
        f.write(f"Total experimental peaks: {total_exp_peaks}\n")
        f.write(f"Annotated peaks: {annotated_peaks}\n")
        f.write(f"Annotation rate: {fragmentation_rate:.4f}\n")
        f.write(f"Total theoretical cleavage sites (sum(L-1)): {total_seq_len}\n")

        writer = csv.writer(f, delimiter="\t")

        writer.writerow(["Theoretical ion type", "Count", "Coverage"])

        for ion, count in ion_counter.most_common():
            coverage = count / total_seq_len
            writer.writerow([ion, count, f"{coverage:.4f}"])                        


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Statical analysis for annotations")
    parser.add_argument("--input", required=True, help="Input annotated msalign file")
    parser.add_argument("--output", required=True, help="Output TSV file")

    args = parser.parse_args()
    ion_counter, total_seq_len, total_exp_peaks, annotated_peaks, fragmentation_rate = count_fragment_stats(args.input)
    write_ion_counts_tsv(ion_counter, total_seq_len, total_exp_peaks, annotated_peaks, fragmentation_rate, args.output)



