import pandas as pd
import numpy as np
import argparse


def extract_mass_shift(df):
    df = df.copy()
    # astype(str) so the .str accessor works even when the column is empty or all-NaN
    # (e.g. no mass-shift rows); "nan"/"None" simply fail the regex and stay NaN.
    mods = df["TOPPIC_unexpected_modifications"].astype(str)

    df["PTM_mass_shift"] = mods.str.extract(r'([+-]?\d+\.\d+)').astype(float)

    # Residue range, written either as [first-last] or a single position [pos]
    # (in which case first and last are the same residue).
    positions = mods.str.extract(r'\[(\d+)(?:-(\d+))?\]')
    df["PTM_first_residue_position"] = pd.to_numeric(
        positions[0], errors="coerce"
    ).astype("Int64")
    df["PTM_last_residue_position"] = pd.to_numeric(
        positions[1].fillna(positions[0]), errors="coerce"
    ).astype("Int64")
    return df


def load_mass_shift_dict(path):
    """Load the common-PTM mass shifts from a TSV file.

    The file has a modification-name column, a mass column, and the list of amino acids
    the modification can occur on. By default these are "Modification", a column whose
    name starts with "Mass", and a third column of residues. The residues may instead be
    appended to the mass cell (e.g. "14.0157 CHKNQRILDEST"), so if there is no separate
    residue column the tokens after the mass are used. Returns an ordered
    {name: (mass, aa_set)} dict, where aa_set is the set of modifiable amino acids
    (empty means no residue restriction).
    """
    df = pd.read_csv(path, sep="\t", dtype=str)
    cols = list(df.columns)
    name_col = "Modification" if "Modification" in cols else cols[0]
    mass_col = next(
        (c for c in cols if str(c).strip().lower().startswith("mass")),
        cols[1],
    )
    # Explicit residue column, else the residues are appended to the mass cell.
    aa_col = next(
        (c for c in cols
         if c not in (name_col, mass_col)
         and any(k in str(c).strip().lower() for k in ("aa", "residue", "amino"))),
        None,
    )
    if aa_col is None and len(cols) > 2:
        aa_col = next((c for c in cols if c not in (name_col, mass_col)), None)

    out = {}
    for _, row in df.iterrows():
        name, raw_mass = row[name_col], row[mass_col]
        if pd.isna(name) or pd.isna(raw_mass):
            continue
        tokens = str(raw_mass).split()          # "<mass> [<residues>]"
        mass = float(tokens[0])
        if aa_col is not None and pd.notna(row[aa_col]):
            aa = str(row[aa_col])
        else:
            aa = "".join(tokens[1:])
        aa_set = {c for c in aa.upper() if c.isalpha()}
        out[str(name)] = (mass, aa_set)
    return out


def match_mass_shift_by_precursor_tol(mass_shift, precursor_mass, subsequence,
                                      mass_shift_dict, ppm_tol=10):
    if pd.isna(mass_shift) or pd.isna(precursor_mass):
        return pd.Series([None, np.nan])

    tol_da = precursor_mass * ppm_tol * 1e-6
    subseq_residues = {c for c in str(subsequence).upper() if c.isalpha()}

    # Pick the closest PTM within tolerance, not merely the first one in dict order
    # that happens to fit (which could otherwise mis-assign when two PTMs are both
    # within tolerance, e.g. for very large precursors). A PTM only qualifies if the
    # modified subsequence contains at least one amino acid it can modify.
    best_name, best_err = None, None
    for name, (theo_mass, aa_set) in mass_shift_dict.items():
        da_error = abs(mass_shift - theo_mass)
        if da_error > tol_da:
            continue
        if aa_set and not (subseq_residues & aa_set):
            continue
        if best_err is None or da_error < best_err:
            best_name, best_err = name, da_error

    if best_name is None:
        return pd.Series([None, np.nan])

    return pd.Series([best_name, best_err])


def _modified_subsequence(sequence, first_pos, last_pos):
    """Residues of TOPPIC_database_sequence spanned by the mass shift, [first, last]
    (1-based, inclusive). Returns "" if the sequence or positions are missing."""
    if pd.isna(sequence) or pd.isna(first_pos) or pd.isna(last_pos):
        return ""
    return str(sequence)[int(first_pos) - 1:int(last_pos)]


def search_common_ptm(df, mass_shift_dict, ppm_tol=10):
    df = df.copy()

    # df.apply(axis=1) on an empty frame returns an empty frame with no columns, which
    # cannot be assigned to the named columns; handle that case explicitly.
    if df.empty:
        df["PTM_matched_mod"] = pd.Series(dtype=object)
        df["PTM_match_da_error"] = pd.Series(dtype=float)
        return df

    df[["PTM_matched_mod", "PTM_match_da_error"]] = df.apply(
        lambda row: match_mass_shift_by_precursor_tol(
            row["PTM_mass_shift"],
            row["MSALIGN_precursor_mass"],
            _modified_subsequence(
                row["TOPPIC_database_sequence"],
                row["PTM_first_residue_position"],
                row["PTM_last_residue_position"],
            ),
            mass_shift_dict,
            ppm_tol
        ),
        axis=1
    )

    # The PTM_* columns describe a match; blank the shift and its residue positions
    # where no PTM was matched.
    unmatched = df["PTM_matched_mod"].isna()
    df.loc[unmatched, "PTM_mass_shift"] = np.nan
    df.loc[unmatched, "PTM_first_residue_position"] = pd.NA
    df.loc[unmatched, "PTM_last_residue_position"] = pd.NA

    return df

def ptm_match_search(filename, mass_shift_filename, out_filename):
    mass_shift_dict = load_mass_shift_dict(mass_shift_filename)

    form_df = pd.read_csv(filename, sep='\t', low_memory=False)
    form_df['MSALIGN_precursor_mass'] = form_df['MSALIGN_precursor_monoisotopic_mass'].apply(lambda x: x.split(':')[0])
    form_df["MSALIGN_precursor_mass"] = pd.to_numeric(form_df["MSALIGN_precursor_mass"], errors="coerce")
    form_df['TOPPIC_sequence_length'] = form_df['TOPPIC_database_sequence'].str.len()

    # Match every row's mass shift (if any) to a common PTM. Rows without a mass shift
    # simply get empty PTM_* columns.
    form_df = extract_mass_shift(form_df)
    form_df = search_common_ptm(form_df, mass_shift_dict)

    # The residue positions are only needed for the residue constraint during matching.
    form_df = form_df.drop(
        columns=["PTM_first_residue_position", "PTM_last_residue_position"],
        errors="ignore",
    )
    form_df.to_csv(out_filename, sep="\t", index=False)
    print(f"Total rows: {len(form_df)}, "
          f"matched to a common PTM: {int(form_df['PTM_matched_mod'].notna().sum())}")
    print(f"Output written to: {out_filename}")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PTM analysis")
    parser.add_argument(
        "-f", "--form",
        required=True,
        help="Input form file(project level)"
    )

    parser.add_argument(
        "-m", "--mass",
        required=True,
        help="Input common mass-shift TSV file (columns: Modification, Mass)"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output TSV file for matched results"
    )

    args = parser.parse_args()

    ptm_match_search(
        filename=args.form,
        mass_shift_filename=args.mass,
        out_filename=args.output
    )



