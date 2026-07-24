#!/usr/bin/env python3
"""
Assign a proteoform identification level using the 5-level classification system.

In top-down proteomics the ambiguity of a proteoform identification is described by a
5-level classification system. Each identification falls into one of eight groups
defined by three conditions, taken from these columns:

    hits  : TOPPIC_number_of_protein_hits                  (== 1 or > 1)
    mods  : TOPPIC_number_of_unexpected_modifications      (== 0 or > 0)
    match : PTM_matched_mod                                (non-empty or empty)

and is assigned a level as follows:

    hits    mods    match       level
    ---------------------------------
    1       0       (any)       1
    > 1     0       (any)       2D
    1       > 0     non-empty   2A
    1       > 0     empty       4
    > 1     > 0     non-empty   3
    > 1     > 0     empty       5

When mods == 0 there is no mass shift to match, so match is empty and the level depends
only on the number of protein hits.

The level is written to a new "Proteoform_classification" column and the resulting table
is saved as a TSV.

Usage:
    python3 proteoform_classification.py <input.tsv> <output.tsv>
"""

import argparse

import numpy as np
import pandas as pd

HITS_COL = "TOPPIC_number_of_protein_hits"
MODS_COL = "TOPPIC_number_of_unexpected_modifications"
MATCH_COL = "PTM_matched_mod"
CLASS_COL = "Proteoform_classification"

# Reporting order of the levels produced by the rules above.
LEVELS = ["1", "2A", "2D", "3", "4", "5"]


def classify_proteoforms(df):
    """Return a copy of df with the Proteoform_classification column added."""
    missing = [c for c in (HITS_COL, MODS_COL, MATCH_COL) if c not in df.columns]
    if missing:
        raise SystemExit(f"missing required column(s): {', '.join(missing)}")

    df = df.copy()

    # Coerce the counts so text or blank values cannot silently mis-classify a row; a
    # missing modification count means no unexpected modification.
    hits = pd.to_numeric(df[HITS_COL], errors="coerce")
    mods = pd.to_numeric(df[MODS_COL], errors="coerce").fillna(0)
    matched = df[MATCH_COL].notna() & (df[MATCH_COL].astype(str).str.strip() != "")

    one_hit, multi_hit = hits == 1, hits > 1
    no_mods, has_mods = mods == 0, mods > 0

    conditions = [
        one_hit & no_mods,                    # level 1
        multi_hit & no_mods,                  # level 2D
        one_hit & has_mods & matched,         # level 2A
        one_hit & has_mods & ~matched,        # level 4
        multi_hit & has_mods & matched,       # level 3
        multi_hit & has_mods & ~matched,      # level 5
    ]
    choices = ["1", "2D", "2A", "4", "3", "5"]

    # Rows whose hit count is missing or below 1 match no rule and stay empty.
    df[CLASS_COL] = np.select(conditions, choices, default="")
    return df


def main():
    parser = argparse.ArgumentParser(
        description="Classify proteoform identifications using the 5-level system"
    )
    parser.add_argument("input", help="Input TSV file with proteoform identifications")
    parser.add_argument("output", help="Output TSV file with the classification column")
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t", low_memory=False)
    df = classify_proteoforms(df)
    df.to_csv(args.output, sep="\t", index=False)

    counts = df[CLASS_COL].value_counts()
    print(f"Total proteoforms: {len(df)}")
    for level in LEVELS:
        print(f"  level {level:<2} {int(counts.get(level, 0)):>10}")
    unclassified = int(counts.get("", 0))
    if unclassified:
        print(f"  unclassified {unclassified:>7}")
    print(f"Output written to: {args.output}")


if __name__ == "__main__":
    main()
