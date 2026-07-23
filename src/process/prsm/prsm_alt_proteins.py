import argparse
import csv
import os
import sys


def build_alternative_proteins(prsm_file):
    """Build a mapping from 'Prsm ID' to its 'Alternative proteins' string.

    In the (non-single) prsm TSV file a scan matched to more than one protein
    spans several consecutive rows sharing the same 'Prsm ID': the first row is
    the primary protein and each following row is an alternative protein with
    only the protein columns filled in. For every alternative protein the
    'Protein accession', 'First residue' and 'Last residue' are joined by ","
    and the proteins are joined by ";". A scan matched to a single protein has
    no alternatives and maps to an empty string.
    """
    groups = {}
    order = []
    with open(prsm_file, newline="") as infile:
        reader = csv.DictReader(infile, delimiter="\t")
        for row in reader:
            prsm_id = row["Prsm ID"]
            if prsm_id not in groups:
                groups[prsm_id] = []
                order.append(prsm_id)
            groups[prsm_id].append(row)

    alt_proteins = {}
    for prsm_id in order:
        rows = groups[prsm_id]
        # The first row is the primary protein; the rest are alternatives.
        alternatives = [
            ",".join([r["Protein accession"], r["First residue"], r["Last residue"]])
            for r in rows[1:]
        ]
        alt_proteins[prsm_id] = ";".join(alternatives)
    return alt_proteins


def basename_data_file_names(rows):
    """Replace each 'Data file name' value with its basename."""
    for row in rows:
        if "Data file name" in row:
            row["Data file name"] = os.path.basename(row["Data file name"])
    return rows


def main():
    parser = argparse.ArgumentParser(
        description="Add a DATASET ID column and an 'Alternative proteins' column "
        "(extracted from the second prsm TSV) to the single prsm TSV file."
    )
    parser.add_argument("single_input", help="Single prsm TSV file path (one row per scan)")
    parser.add_argument("prsm_input", help="Prsm TSV file path (alternative proteins on extra rows)")
    parser.add_argument("dataset_id", help="Value to use for the DATASET ID column")
    parser.add_argument("output", help="Output TSV file path")
    args = parser.parse_args()

    alt_proteins = build_alternative_proteins(args.prsm_input)

    with open(args.single_input, newline="") as infile:
        reader = csv.DictReader(infile, delimiter="\t")
        fieldnames = ["DATASET ID"] + reader.fieldnames + ["Alternative proteins"]

        with open(args.output, "w", newline="") as outfile:
            writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            rows = [
                dict(
                    row,
                    **{
                        "DATASET ID": args.dataset_id,
                        "Alternative proteins": alt_proteins.get(row["Prsm ID"], ""),
                    },
                )
                for row in reader
            ]
            rows = basename_data_file_names(rows)
            for row in rows:
                writer.writerow(row)


if __name__ == "__main__":
    main()
