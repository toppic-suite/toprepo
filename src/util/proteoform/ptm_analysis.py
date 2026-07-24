import pandas as pd
import numpy as np
import argparse


def form_type_gen(df):
    # Proteoforms that carry a mass shift (>= 1 unexpected modification). Coerce the
    # count to numeric and treat a missing/blank value as 0, so a NaN count is never
    # mistaken for a mass shift (which would otherwise slip through the != 0 test).
    n_unexpected = pd.to_numeric(
        df['TOPPIC_number_of_unexpected_modifications'], errors='coerce'
    ).fillna(0)
    return df[n_unexpected != 0]


def extract_mass_shift(df):
    df = df.copy()
    # astype(str) so the .str accessor works even when the column is empty or all-NaN
    # (e.g. no mass-shift rows); "nan"/"None" simply fail the regex and stay NaN.
    df["PTM_mass_shift"] = (
        df["TOPPIC_unexpected_modifications"]
        .astype(str)
        .str.extract(r'([+-]?\d+\.\d+)')
        .astype(float)
    )
    return df


def load_mass_shift_dict(path):
    """Load the common-PTM mass shifts from a TSV file.

    The file has a modification-name column and a mass column; by default these are
    "Modification" and a column whose name starts with "Mass" (as in common_mass.tsv),
    falling back to the first two columns. The mass cell may carry the applicable amino
    acids after the value (e.g. "14.0157 CHKNQRILDEST"), so the leading numeric token is
    taken as the mass. Returns an ordered {name: mass} dict.
    """
    df = pd.read_csv(path, sep="\t")
    cols = list(df.columns)
    name_col = "Modification" if "Modification" in cols else cols[0]
    mass_col = next(
        (c for c in cols if str(c).strip().lower().startswith("mass")),
        cols[1],
    )
    out = {}
    for n, m in zip(df[name_col], df[mass_col]):
        if pd.isna(n) or pd.isna(m):
            continue
        out[str(n)] = float(str(m).split()[0])  # value may be "<mass> <residues>"
    return out


def match_mass_shift_by_precursor_tol(mass_shift, precursor_mass, mass_shift_dict, ppm_tol=10):
    if pd.isna(mass_shift) or pd.isna(precursor_mass):
        return pd.Series([None, np.nan])

    tol_da = precursor_mass * ppm_tol * 1e-6

    # Pick the closest PTM within tolerance, not merely the first one in dict order
    # that happens to fit (which could otherwise mis-assign when two PTMs are both
    # within tolerance, e.g. for very large precursors).
    best_name, best_err = None, None
    for name, theo_mass in mass_shift_dict.items():
        da_error = abs(mass_shift - theo_mass)
        if da_error <= tol_da and (best_err is None or da_error < best_err):
            best_name, best_err = name, da_error

    if best_name is None:
        return pd.Series([None, np.nan])

    return pd.Series([best_name, best_err])


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
            mass_shift_dict,
            ppm_tol
        ),
        axis=1
    )

    # The PTM_* columns describe a match; blank the shift where no PTM was matched.
    df.loc[df["PTM_matched_mod"].isna(), "PTM_mass_shift"] = np.nan

    return df

def ptm_match_search(filename, mass_shift_filename, out_filename):
    mass_shift_dict = load_mass_shift_dict(mass_shift_filename)

    form_df = pd.read_csv(filename, sep='\t', low_memory=False)
    form_df['MSALIGN_precursor_mass'] = form_df['MSALIGN_precursor_monoisotopic_mass'].apply(lambda x: x.split(':')[0])
    form_df["MSALIGN_precursor_mass"] = pd.to_numeric(form_df["MSALIGN_precursor_mass"], errors="coerce")
    form_df['TOPPIC_sequence_length'] = form_df['TOPPIC_database_sequence'].str.len()

    # Keep the proteoforms that carry a mass shift and match each shift to a common PTM.
    form_df2 = form_type_gen(form_df)
    form_df2 = extract_mass_shift(form_df2)
    form_df2 = search_common_ptm(form_df2, mass_shift_dict)

    form_df2.to_csv(out_filename, sep="\t", index=False)
    print(f"Mass-shift proteoforms: {len(form_df2)}, "
          f"matched to a common PTM: {int(form_df2['PTM_matched_mod'].notna().sum())}")
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



