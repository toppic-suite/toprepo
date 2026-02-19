import argparse
import csv
import os
import sys


def basename_data_file_names(rows):
    """Replace each 'Data file name' value with its basename."""
    for row in rows:
        if "Data file name" in row:
            row["Data file name"] = os.path.basename(row["Data file name"])
    return rows


def main():
    parser = argparse.ArgumentParser(description="Add a DATASET_ID column as the first column of a TSV file.")
    parser.add_argument("input", help="Input TSV file path")
    parser.add_argument("dataset_id", help="Value to use for the DATASET_ID column")
    parser.add_argument("-o", "--output", help="Output TSV file path (default: stdout)")
    args = parser.parse_args()

    with open(args.input, newline="") as infile:
        reader = csv.DictReader(infile, delimiter="\t")
        fieldnames = ["DATASET ID"] + reader.fieldnames

        outfile = open(args.output, "w", newline="") if args.output else sys.stdout
        try:
            writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            rows = [dict(row, **{"DATASET ID": args.dataset_id}) for row in reader]
            rows = basename_data_file_names(rows)
            for row in rows:
                writer.writerow(row)
        finally:
            if args.output:
                outfile.close()


if __name__ == "__main__":
    main()
