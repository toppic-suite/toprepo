#!/usr/bin/env python3
"""
Randomly select n scans without replacement from an msalign file.

Usage:
    python random_select_msalign.py <input_file> <output_file> <n>

Arguments:
    input_file  - Path to the input msalign file
    output_file - Path to the output msalign file
    n           - Number of scans to randomly select
"""

import argparse
import random
import sys


def parse_msalign_scans(filepath):
    """Parse an msalign file and return a list of scans.

    Each scan is stored as a string including BEGIN IONS and END IONS.
    """
    scans = []
    current_scan_lines = []
    in_scan = False

    with open(filepath, 'r') as f:
        for line in f:
            if line.strip() == 'BEGIN IONS':
                in_scan = True
                current_scan_lines = [line]
            elif line.strip() == 'END IONS':
                current_scan_lines.append(line)
                scans.append(''.join(current_scan_lines))
                current_scan_lines = []
                in_scan = False
            elif in_scan:
                current_scan_lines.append(line)

    return scans


def main():
    parser = argparse.ArgumentParser(
        description='Randomly select n scans without replacement from an msalign file.'
    )
    parser.add_argument('input_file', help='Path to the input msalign file')
    parser.add_argument('output_file', help='Path to the output msalign file')
    parser.add_argument('n', type=int, help='Number of scans to randomly select')

    args = parser.parse_args()

    # Parse input file
    scans = parse_msalign_scans(args.input_file)
    total_scans = len(scans)

    if args.n > total_scans:
        print(f"Error: Requested {args.n} scans, but input file only contains {total_scans} scans.",
              file=sys.stderr)
        sys.exit(1)

    if args.n <= 0:
        print("Error: n must be a positive integer.", file=sys.stderr)
        sys.exit(1)

    # Randomly select n scans without replacement
    selected_scans = random.sample(scans, args.n)

    # Write selected scans to output file
    with open(args.output_file, 'w') as f:
        for i, scan in enumerate(selected_scans):
            f.write(scan)
            # Add blank line between scans (but not after the last one)
            if i < len(selected_scans) - 1:
                f.write('\n')

    print(f"Successfully selected {args.n} scans from {total_scans} total scans.")
    print(f"Output written to: {args.output_file}")


if __name__ == '__main__':
    main()
