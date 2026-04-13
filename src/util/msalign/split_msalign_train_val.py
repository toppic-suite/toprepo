#!/usr/bin/env python3
"""
Split an msalign file into training and validation sets based on
DATABASE_SEQUENCE membership in two proteoform list files.

Train list format:  SEQUENCE=<sequence>
Val list format:    DATABASE_SEQUENCE=<sequence>
msalign spectrum field matched: DATABASE_SEQUENCE=<sequence>
"""

import argparse
import os
import sys


def load_sequences(path):
    """Return a set of bare sequences from a proteoform list file.
    as well as plain bare sequences.
    """
    sequences = set()
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("DATABASE_SEQUENCE="):
                sequences.add(line[len("DATABASE_SEQUENCE="):])
    return sequences


def iter_spectra(path):
    """Yield (lines, db_sequence) tuples for each spectrum block."""
    block = []
    db_seq = None
    in_block = False

    with open(path) as fh:
        for line in fh:
            raw = line.rstrip("\n")

            if raw == "BEGIN IONS":
                block = [raw]
                db_seq = None
                in_block = True
                continue

            if in_block:
                block.append(raw)
                if raw.startswith("DATABASE_SEQUENCE="):
                    db_seq = raw[len("DATABASE_SEQUENCE="):]
                if raw == "END IONS":
                    yield block, db_seq
                    block = []
                    db_seq = None
                    in_block = False


def split_msalign(msalign_path, train_list_path, val_list_path,
                  train_out_path, val_out_path):
    print(f"Loading train sequences from: {train_list_path}")
    train_seqs = load_sequences(train_list_path)
    print(f"  {len(train_seqs):,} unique train sequences")

    print(f"Loading val sequences from:   {val_list_path}")
    val_seqs = load_sequences(val_list_path)
    print(f"  {len(val_seqs):,} unique val sequences")

    overlap = train_seqs & val_seqs
    if overlap:
        print(f"WARNING: {len(overlap):,} sequence(s) appear in both lists; "
              "they will be written to BOTH output files.")

    counts = {"train": 0, "val": 0, "skipped": 0}

    with open(train_out_path, "w") as train_fh, \
         open(val_out_path, "w") as val_fh:

        for block, db_seq in iter_spectra(msalign_path):
            in_train = db_seq in train_seqs
            in_val = db_seq in val_seqs

            if in_train:
                train_fh.write("\n".join(block) + "\n\n")
                counts["train"] += 1
            if in_val:
                val_fh.write("\n".join(block) + "\n\n")
                counts["val"] += 1
            if not in_train and not in_val:
                counts["skipped"] += 1

    total = counts["train"] + counts["val"] + counts["skipped"]
    print(f"\nDone. Processed {total} spectra total:")
    print(f"  -> train : {counts['train']:,}  ({train_out_path})")
    print(f"  -> val   : {counts['val']:,}  ({val_out_path})")
    print(f"  -> skipped (sequence not in either list): {counts['skipped']:,}")


def main():
    parser = argparse.ArgumentParser(
        description="Split an msalign file into train / val sets.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("msalign", help="Input msalign file")
    parser.add_argument("train_list", help="File listing train proteoform sequences")
    parser.add_argument("val_list", help="File listing validation proteoform sequences")
    args = parser.parse_args()

    for path in (args.msalign, args.train_list, args.val_list):
        if not os.path.exists(path):
            sys.exit(f"ERROR: file not found: {path}")

    basename = os.path.splitext(os.path.basename(args.msalign))[0]
    train_out = basename + "_train.msalign"
    val_out = basename + "_val.msalign"

    split_msalign(
        args.msalign,
        args.train_list,
        args.val_list,
        train_out,
        val_out,
    )


if __name__ == "__main__":
    main()
