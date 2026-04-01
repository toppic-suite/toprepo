#!/usr/bin/env python3
"""
Filter rows in combined_info TSV using proteoform_info TSV as a filter.

Usage:
    python filter_combined_info.py <proteoform_info.tsv> <combined_info.tsv> <output.tsv>
"""

import sys
import csv


def build_filter_set(proteoform_file):
    """Read proteoform_info file and build a set of (DATASET_ID, Data_file_name, Scans) tuples."""
    filter_set = set()
    with open(proteoform_file, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            key = (row["DATASET ID"], row["Data file name"], row["Scan(s)"])
            filter_set.add(key)
    return filter_set


def filter_combined(combined_file, output_file, filter_set):
    """Filter combined_info file, keeping rows whose key is in filter_set."""
    kept = 0
    removed = 0
    with open(combined_file, newline='') as fin, open(output_file, 'w', newline='') as fout:
        reader = csv.DictReader(fin, delimiter='\t')
        writer = csv.DictWriter(fout, fieldnames=reader.fieldnames, delimiter='\t')
        writer.writeheader()
        for row in reader:
            key = (row["DATASET_id"], row["MSALIGN_file_name"], row["MZML_ms2_scan"])
            if key in filter_set:
                writer.writerow(row)
                kept += 1
            else:
                removed += 1
    return kept, removed


def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <proteoform_info.tsv> <combined_info.tsv> <output.tsv>")
        sys.exit(1)

    proteoform_file, combined_file, output_file = sys.argv[1], sys.argv[2], sys.argv[3]

    print(f"Reading filter keys from: {proteoform_file}")
    filter_set = build_filter_set(proteoform_file)
    print(f"  Loaded {len(filter_set)} unique (DATASET_ID, Data file name, Scans) keys")

    print(f"Filtering: {combined_file}")
    kept, removed = filter_combined(combined_file, output_file, filter_set)
    print(f"  Kept: {kept} rows, Removed: {removed} rows")
    print(f"Output written to: {output_file}")


if __name__ == "__main__":
    main()
