"""
remove duplicate proteoform
"""

import pandas as pd
import sys
from collections import defaultdict


MASS_TOL = 2.2
BUCKET = MASS_TOL * 2  # bucket width so 3 adjacent buckets always cover ±MASS_TOL
OVERLAP_THRESH = 0.70


def calculate_overlap_percentage(r1_start, r1_end, r2_start, r2_end):
    overlap_start = max(r1_start, r2_start)
    overlap_end   = min(r1_end, r2_end)
    if overlap_start > overlap_end:
        return 0.0
    overlap_len = overlap_end - overlap_start + 1
    shorter_len = min(r1_end - r1_start + 1, r2_end - r2_start + 1)
    return overlap_len / shorter_len if shorter_len > 0 else 0.0


def deduplicate_group(df):
    # Pre-convert types once rather than per-comparison
    df = df.copy()
    df["TOPPIC_e-value"]                = pd.to_numeric(df["TOPPIC_e-value"],                errors="coerce")
    df["MSALIGN_precursor_mass"]        = pd.to_numeric(df["MSALIGN_precursor_mass"],        errors="coerce")
    df["TOPPIC_first_residue_position"] = pd.to_numeric(df["TOPPIC_first_residue_position"], errors="coerce")
    df["TOPPIC_last_residue_position"]  = pd.to_numeric(df["TOPPIC_last_residue_position"],  errors="coerce")

    files = []
    for _, sub in df.groupby("MSALIGN_file_name"):
        files.append(sub.sort_values("TOPPIC_e-value").to_dict("records"))

    # Survivors list with lazy-deletion via a boolean flag.
    # A mass-bucketed index maps bucket_id -> [survivor indices] for O(K) dup lookup.
    survivors: list[dict] = []
    deleted:   list[bool] = []
    bucket_idx: dict[int, list[int]] = defaultdict(list)

    def _bkt(mass: float) -> int:
        return int(mass / BUCKET)

    def _insert(row: dict) -> None:
        idx = len(survivors)
        survivors.append(row)
        deleted.append(False)
        b = _bkt(row["MSALIGN_precursor_mass"])
        bucket_idx[b - 1].append(idx)
        bucket_idx[b    ].append(idx)
        bucket_idx[b + 1].append(idx)

    def _find_dup(row: dict) -> int:
        """Return index of a duplicate survivor, or -1."""
        m         = row["MSALIGN_precursor_mass"]
        row_file  = row["MSALIGN_file_name"]
        row_start = row["TOPPIC_first_residue_position"]
        row_end   = row["TOPPIC_last_residue_position"]
        seen: set[int] = set()
        for idx in bucket_idx[_bkt(m)]:  # 3 buckets already merged at insert time
            if idx in seen or deleted[idx]:
                continue
            seen.add(idx)
            s = survivors[idx]
            if abs(s["MSALIGN_precursor_mass"] - m) > MASS_TOL:
                continue
            if s["MSALIGN_file_name"] == row_file:
                continue
            if calculate_overlap_percentage(
                s["TOPPIC_first_residue_position"], s["TOPPIC_last_residue_position"],
                row_start, row_end
            ) >= OVERLAP_THRESH:
                return idx
        return -1

    for file_rows in files:
        for row in file_rows:
            dup_idx = _find_dup(row)
            if dup_idx >= 0:
                if row["TOPPIC_e-value"] < survivors[dup_idx]["TOPPIC_e-value"]:
                    deleted[dup_idx] = True
                    _insert(row)
            else:
                _insert(row)

    return pd.DataFrame([survivors[i] for i in range(len(survivors)) if not deleted[i]])


def main(input_file, output_file):
    df = pd.read_csv(input_file, sep="\t", low_memory=False)
    print(f"{len(df)}")
    df["MSALIGN_precursor_mass"] = df["MSALIGN_precursor_monoisotopic_mass"].str.split(":").str[0]

    groups = [deduplicate_group(g) for _, g in df.groupby("TOPPIC_protein_accession")]
    deduped_chunk = pd.concat(groups, ignore_index=True)
    print(f"Producing {len(deduped_chunk)} deduplicated rows")

    deduped_chunk.drop(columns=["MSALIGN_precursor_mass"], inplace=True)
    num_ptms = deduped_chunk["TOPPIC_proteoform_id"].notna().sum()

    deduped_chunk.to_csv(output_file, sep="\t", index=False)

    print(f"Total spectra (rows): {len(df)}, total remove duplicates: {len(deduped_chunk)}")
    print(f"Rows with proteoform ID: {num_ptms}")
    print("Done.")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py <input_file> <output_file>")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2])
