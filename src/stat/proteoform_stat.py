#!/usr/bin/env python3
"""
Summarise a proteoform info / proteoform table TSV.

For each input file it counts, over the proteoform-identification rows:

  * total rows (proteoform IDs)
  * rows with a mass shift and rows without one
  * each of those split by number of protein hits (== 1 vs > 1)

A "mass shift" is an unexpected modification reported by TopPIC, i.e. a non-zero value
in the "number of unexpected modifications" column. The four leaf buckets

    with-shift & hits==1, with-shift & hits>1,
    without-shift & hits==1, without-shift & hits>1

partition the counted rows, so they always add up to the total (minus any rows whose
mass-shift or protein-hit field is blank/unparseable, which are reported separately).

Works on the aggregated proteoform table (TOPPIC_* namespaced columns, e.g. the files in
05_merged_info/04_proteoform_combined) and on a raw TopPIC proteoform info TSV
(#Protein hits / #unexpected modifications columns).

Usage:
    python3 proteoform_stat.py <file.tsv> [<file2.tsv> ...]
"""

import csv
import sys

# Column name candidates, tried in order; first present in the header wins.
_UNEXPECTED_COLS = ["TOPPIC_number_of_unexpected_modifications", "#unexpected modifications"]
_HITS_COLS       = ["TOPPIC_number_of_protein_hits", "#Protein hits"]
_ID_COLS         = ["TOPPIC_proteoform_id", "Proteoform ID"]

# Proteoform strings / database sequences can be long; lift the csv field cap.
csv.field_size_limit(min(sys.maxsize, 2**31 - 1))


def _resolve(fieldnames, candidates, path):
    for name in candidates:
        if name in fieldnames:
            return name
    raise SystemExit(
        f"{path}: none of the expected columns {candidates} found in header"
    )


def _to_int(value):
    """Return int(value), or None if it is blank / not an integer."""
    if value is None:
        return None
    value = value.strip()
    if value == "":
        return None
    try:
        return int(value)
    except ValueError:
        return None


def summarise(path):
    counts = {
        "total": 0,
        "with_shift": 0,
        "without_shift": 0,
        "with_shift_hits1": 0,
        "with_shift_hits_gt1": 0,
        "without_shift_hits1": 0,
        "without_shift_hits_gt1": 0,
        "id_blank": 0,
        "uncounted": 0,   # mass-shift or protein-hit field blank/unparseable
    }
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise SystemExit(f"{path}: empty file or missing header")
        unexpected_col = _resolve(reader.fieldnames, _UNEXPECTED_COLS, path)
        hits_col       = _resolve(reader.fieldnames, _HITS_COLS, path)
        id_col         = next((c for c in _ID_COLS if c in reader.fieldnames), None)

        for row in reader:
            counts["total"] += 1

            if id_col is not None:
                pid = row.get(id_col)
                if pid is None or pid.strip() == "":
                    counts["id_blank"] += 1

            shifts = _to_int(row.get(unexpected_col))
            hits = _to_int(row.get(hits_col))
            if shifts is None or hits is None:
                counts["uncounted"] += 1
                continue

            has_shift = shifts > 0
            one_hit = hits == 1
            if has_shift:
                counts["with_shift"] += 1
                counts["with_shift_hits1" if one_hit else "with_shift_hits_gt1"] += 1
            else:
                counts["without_shift"] += 1
                counts["without_shift_hits1" if one_hit else "without_shift_hits_gt1"] += 1

    counts["_unexpected_col"] = unexpected_col
    counts["_hits_col"] = hits_col
    return counts


def _pct(part, whole):
    return f"{100.0 * part / whole:5.1f}%" if whole else "  n/a "


def report(path, c):
    total = c["total"]
    print(f"File: {path}")
    print(f"  mass-shift column: {c['_unexpected_col']}")
    print(f"  protein-hit column: {c['_hits_col']}")
    print(f"  {'Total proteoform IDs':<34}{total:>12}")
    print(f"  {'With mass shift':<34}{c['with_shift']:>12}  ({_pct(c['with_shift'], total)})")
    print(f"    {'protein hits = 1':<32}{c['with_shift_hits1']:>12}")
    print(f"    {'protein hits > 1':<32}{c['with_shift_hits_gt1']:>12}")
    print(f"  {'Without mass shift':<34}{c['without_shift']:>12}  ({_pct(c['without_shift'], total)})")
    print(f"    {'protein hits = 1':<32}{c['without_shift_hits1']:>12}")
    print(f"    {'protein hits > 1':<32}{c['without_shift_hits_gt1']:>12}")
    if c["id_blank"]:
        print(f"  {'(rows with blank proteoform ID)':<34}{c['id_blank']:>12}")
    if c["uncounted"]:
        print(f"  {'(rows uncounted: blank/bad fields)':<34}{c['uncounted']:>12}")


def main(argv):
    if len(argv) < 2:
        print(f"Usage: {argv[0]} <file.tsv> [<file2.tsv> ...]")
        return 1
    for i, path in enumerate(argv[1:]):
        if i:
            print()
        report(path, summarise(path))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
