"""
remove duplicate proteoform (fast, output-identical to remove_duplicate_proteoform.py)

This is a drop-in speedup of remove_duplicate_proteoform.py: it produces the exact
same output rows in the same order. It differs from remove_duplicate_proteoform_fast.py,
which uses a different (global-e-value, no-replace) algorithm and therefore can produce a
slightly different set of survivors.

Speedups, none of which change the result:
  * The O(n^2) survivor scan runs on plain Python primitives (floats/ints/strings)
    extracted once per group, instead of pandas Series access via DataFrame.iterrows().
  * Numeric conversions (precursor mass, residue positions, e-value) are done once up
    front rather than on every pairwise comparison.
  * A precursor-mass window (sorted list + binary search) limits each row's comparisons
    to survivors within +-2.2 Da. Since is_duplicate() requires |dmass| <= 2.2, survivors
    outside that window can never match, so pruning them cannot change the outcome. The
    window is widened by a tiny epsilon so it is always a superset of the true match set,
    and the exact |dmass| <= 2.2 test is still applied to every candidate.

The original's semantics are preserved exactly:
  * rows are processed grouped by MSALIGN_file_name (file-name order), and within each
    file sorted by e-value -- the same order the original builds and iterates;
  * for each row the FIRST survivor in insertion order that is a duplicate is selected;
  * if the incoming row has a strictly smaller e-value it REPLACES that survivor (keeping
    the slot, so later comparisons see the replacement); otherwise the incoming row is
    dropped; a row that matches no survivor is appended;
  * output row order is the survivor insertion order (with in-place replacements).
"""

import bisect

import pandas as pd
import sys


# ---- window epsilon --------------------------------------------------------
# Only used to widen the mass window so it is a guaranteed superset of the
# candidates the exact |dmass| <= 2.2 test would accept; it never relaxes that test.
_MASS_TOL = 2.2
_EPS = 1e-6


def calculate_overlap_percentage(r1_start, r1_end, r2_start, r2_end):
    try:
        r1_start, r1_end = int(r1_start), int(r1_end)
        r2_start, r2_end = int(r2_start), int(r2_end)

        overlap_start = max(r1_start, r2_start)
        overlap_end   = min(r1_end, r2_end)

        if overlap_start > overlap_end:
            return 0.0

        overlap_len = overlap_end - overlap_start + 1
        len1 = r1_end - r1_start + 1
        len2 = r2_end - r2_start + 1
        shorter_len = min(len1, len2)
        return overlap_len / shorter_len if shorter_len > 0 else 0
    except:
        return 0.0


def is_duplicate(a, b):
    """Return True if b is duplicate of a and should be removed.

    Kept for API parity with remove_duplicate_proteoform.py; deduplicate_group below
    inlines the same predicate over precomputed primitives for speed.
    """
    if a["MSALIGN_file_name"] == b["MSALIGN_file_name"]:
        return False

    if a["TOPPIC_protein_accession"] != b["TOPPIC_protein_accession"]:
        return False

    if abs(float(a['MSALIGN_precursor_mass']) -
           float(b['MSALIGN_precursor_mass'])) > 2.2:
        return False

    ov = calculate_overlap_percentage(
        a["TOPPIC_first_residue_position"], a["TOPPIC_last_residue_position"],
        b["TOPPIC_first_residue_position"], b["TOPPIC_last_residue_position"]
    )
    if ov < 0.70:
        return False

    return True


def _safe_int(x):
    """int(x) with the same fall-through-to-invalid behaviour as the original's
    try/except inside calculate_overlap_percentage."""
    try:
        return int(x)
    except:
        return None


def deduplicate_group(df):
    df = df.copy()
    df["TOPPIC_e-value"] = df["TOPPIC_e-value"].astype(float)
    df = df.reset_index(drop=True)

    # Exact processing order: iterate files in MSALIGN_file_name order (groupby sorts
    # keys), and within each file sort by e-value -- identical to the original, which
    # builds `files` the same way and iterates them in the same sequence.
    order_positions = []
    for _f, sub in df.groupby("MSALIGN_file_name"):
        order_positions.extend(sub.sort_values("TOPPIC_e-value").index.tolist())

    # Primitive columns, indexed by row position (0..n-1). Converting once here is what
    # replaces the per-comparison float()/int() and Series lookups in the original.
    file_col  = df["MSALIGN_file_name"].tolist()
    mass_col  = [float(x) for x in df["MSALIGN_precursor_mass"].tolist()]
    ev_col    = df["TOPPIC_e-value"].tolist()          # already float
    start_col = [_safe_int(x) for x in df["TOPPIC_first_residue_position"].tolist()]
    end_col   = [_safe_int(x) for x in df["TOPPIC_last_residue_position"].tolist()]
    resok_col = [start_col[i] is not None and end_col[i] is not None
                 for i in range(len(df))]

    # Survivors in insertion order (slot index == the original's enumerate index).
    s_pos, s_file, s_mass, s_ev, s_start, s_end, s_resok = [], [], [], [], [], [], []

    # Survivor masses kept sorted for the +-2.2 Da window, parallel to survivor slot ids.
    win_mass, win_slot = [], []

    tol = _MASS_TOL
    for p in order_positions:
        rf  = file_col[p]
        rm  = mass_col[p]
        rev = ev_col[p]
        rs  = start_col[p]
        re  = end_col[p]
        rok = resok_col[p]

        lo = bisect.bisect_left(win_mass, rm - tol - _EPS)
        hi = bisect.bisect_right(win_mass, rm + tol + _EPS)

        # First survivor in insertion order (smallest slot id) that is a duplicate.
        best = -1
        for w in range(lo, hi):
            j = win_slot[w]
            if j >= best and best != -1:
                # A smaller slot id already matched; a larger one can't beat it.
                continue
            if s_file[j] == rf:
                continue
            if abs(s_mass[j] - rm) > tol:
                continue
            if s_resok[j] and rok:
                a_s = s_start[j]; a_e = s_end[j]
                o_s = a_s if a_s > rs else rs
                o_e = a_e if a_e < re else re
                if o_s > o_e:
                    ov = 0.0
                else:
                    overlap_len = o_e - o_s + 1
                    len1 = a_e - a_s + 1
                    len2 = re - rs + 1
                    shorter = len1 if len1 < len2 else len2
                    ov = overlap_len / shorter if shorter > 0 else 0.0
            else:
                ov = 0.0
            if ov < 0.70:
                continue
            if best == -1 or j < best:
                best = j

        if best != -1:
            # Duplicate of survivor `best`: replace it iff the incoming row has a
            # strictly smaller e-value, otherwise drop the incoming row.
            if rev < s_ev[best]:
                old_mass = s_mass[best]
                s_pos[best]   = p
                s_file[best]  = rf
                s_mass[best]  = rm
                s_ev[best]    = rev
                s_start[best] = rs
                s_end[best]   = re
                s_resok[best] = rok
                # Move slot `best` to its new mass in the sorted window.
                k = bisect.bisect_left(win_mass, old_mass)
                while win_slot[k] != best or win_mass[k] != old_mass:
                    k += 1
                win_mass.pop(k); win_slot.pop(k)
                ins = bisect.bisect_left(win_mass, rm)
                win_mass.insert(ins, rm); win_slot.insert(ins, best)
            # else: drop
        else:
            # No duplicate: append as a new survivor.
            slot = len(s_pos)
            s_pos.append(p); s_file.append(rf); s_mass.append(rm); s_ev.append(rev)
            s_start.append(rs); s_end.append(re); s_resok.append(rok)
            ins = bisect.bisect_left(win_mass, rm)
            win_mass.insert(ins, rm); win_slot.insert(ins, slot)

    return df.iloc[s_pos]


def main(input_file, output_file):
    df = pd.read_csv(
        input_file,
        sep="\t",
        low_memory=False
    )
    print(f"{len(df)}")
    df['MSALIGN_precursor_mass'] = df['MSALIGN_precursor_monoisotopic_mass'].apply(lambda x: x.split(':')[0])

    total_rows = len(df)
    processed = 0
    groups = []
    for acc, group in df.groupby("TOPPIC_protein_accession"):
        groups.append(deduplicate_group(group))
        processed += len(group)
        if processed // 1000 > (processed - len(group)) // 1000:
            print(f"Processed {processed}/{total_rows} records ({processed/total_rows:.0%})", flush=True)

    deduped_chunk = pd.concat(groups, ignore_index=True)
    print(f"Producing {len(deduped_chunk)} deduplicated rows")

    deduped_chunk.drop(columns=["MSALIGN_precursor_mass"], axis=1, inplace=True)
    num_ptms = deduped_chunk['TOPPIC_proteoform_id'].notna().sum()

    deduped_chunk.to_csv(
        output_file,
        sep="\t",
        index=False)

    print(f"Total spectra (rows): {len(df)}, total remove duplicates: {len(deduped_chunk)}")
    print(f"Rows with proteoform ID: {num_ptms}")
    print("Done.")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    main(input_file, output_file)
