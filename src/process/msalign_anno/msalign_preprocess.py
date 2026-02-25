#!/usr/bin/env python3
import argparse
import sys
import os

def main():
    parser = argparse.ArgumentParser(
        description="Add DATASET_ID to each spectrum in an msalign file."
    )
    parser.add_argument("msalign_file", help="Input msalign file")
    parser.add_argument("dataset_id", help="Dataset ID string to add")
    parser.add_argument("output_file", help="Output msalign file")
    args = parser.parse_args()

    with open(args.msalign_file, "r") as fin, open(args.output_file, "w") as fout:
        for line in fin:
            if line.strip() == "BEGIN IONS":
                fout.write(line)
                fout.write(f"DATASET_ID={args.dataset_id}\n")
            elif line.find("=") != -1:
                if line.startswith("FILE_NAME="):
                    line = "MZML_" + line
                    fout.write(line)
                    fout.write(f"MSALIGN_FILE_NAME={os.path.basename(args.msalign_file)}\n")
                elif line.startswith("SPECTRUM_ID="):
                    parts = line.strip().split("=", 1)
                    if len(parts) == 2:
                        scan_info = parts[1]
                        fout.write(f"MS2_ID={scan_info}\n")
                elif line.startswith("SCANS="):
                    parts = line.strip().split("=", 1)
                    if len(parts) == 2:
                        scan_info = parts[1]
                        fout.write(f"MS2_SCAN={scan_info}\n")
                elif line.startswith("TITLE="):
                    continue
                elif line.startswith("RETENTION_TIME="):
                    line = "MS2_" + line
                    fout.write(line)
                elif line.startswith("MS_ONE_ID="):
                    parts = line.strip().split("=", 1)
                    if len(parts) == 2:
                        scan_info = parts[1]
                        fout.write(f"MS1_ID={scan_info}\n")
                elif line.startswith("MS_ONE_SCAN="):
                    parts = line.strip().split("=", 1)
                    if len(parts) == 2:
                        scan_info = parts[1]
                        fout.write(f"MS1_SCAN={scan_info}\n")
                elif line.startswith("PRECURSOR_MASS="):
                    parts = line.strip().split("=", 1)
                    if len(parts) == 2:
                        scan_info = parts[1]
                        fout.write(f"PRECURSOR_MONOISOTOPIC_MASS={scan_info}\n")
                elif line.startswith("PRECURSOR_MZ="):
                    parts = line.strip().split("=", 1)
                    if len(parts) == 2:
                        scan_info = parts[1]
                        fout.write(f"PRECURSOR_MONOISOTOPIC_MZ={scan_info}\n")
                else:
                    fout.write(line)
            else:
                fout.write(line) 
if __name__ == "__main__":
    main()
