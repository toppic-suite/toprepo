#!/usr/bin/env python3
"""Produce nested random subsets of an msalign file with a fixed seed.

Draws one seeded random permutation of the scans, then writes each requested
size from the front of that permutation. Because every subset is a prefix of the
same shuffle, smaller subsets are guaranteed to be contained in larger ones
(e.g. 100K subset of 200K subset of all). Output format matches
random_select_msalign.py: scans separated by a blank line, none after the last.

Usage:
    python nested_select_msalign.py <input_file> <seed> <size>:<output_file> [<size>:<output_file> ...]
"""
import random
import sys


def parse_msalign_scans(filepath):
    scans = []
    current = []
    in_scan = False
    with open(filepath, "r") as f:
        for line in f:
            if line.strip() == "BEGIN IONS":
                in_scan = True
                current = [line]
            elif line.strip() == "END IONS":
                current.append(line)
                scans.append("".join(current))
                current = []
                in_scan = False
            elif in_scan:
                current.append(line)
    return scans


def write_scans(scans, output_file):
    with open(output_file, "w") as f:
        for i, scan in enumerate(scans):
            f.write(scan)
            if i < len(scans) - 1:
                f.write("\n")


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    input_file = sys.argv[1]
    seed = int(sys.argv[2])
    jobs = []
    for spec in sys.argv[3:]:
        size_str, out = spec.split(":", 1)
        jobs.append((int(size_str), out))

    scans = parse_msalign_scans(input_file)
    total = len(scans)
    max_size = max(size for size, _ in jobs)
    if max_size > total:
        sys.exit(f"Error: requested {max_size} scans, input only has {total}.")

    random.seed(seed)
    order = list(range(total))
    random.shuffle(order)  # one shuffle -> all subsets are nested prefixes

    for size, out in sorted(jobs):
        subset = [scans[i] for i in order[:size]]
        write_scans(subset, out)
        print(f"seed={seed}: selected {size} of {total} scans -> {out}")


if __name__ == "__main__":
    main()
