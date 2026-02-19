#!/usr/bin/env python3
import argparse
import sys


def main():
    parser = argparse.ArgumentParser(
        description="Add DATASET_ID to each spectrum in an mgf file."
    )
    parser.add_argument("mgf_file", help="Input mgf file")
    parser.add_argument("dataset_id", help="Dataset ID string to add")
    parser.add_argument("output_file", help="Output mgf file")
    args = parser.parse_args()

    with open(args.mgf_file, "r") as fin, open(args.output_file, "w") as fout:
        for line in fin:
            if line.strip() == "BEGIN IONS":
                fout.write(line)
                fout.write(f"DATASET_ID={args.dataset_id}\n")
            else:
                fout.write(line) 
if __name__ == "__main__":
    main()
