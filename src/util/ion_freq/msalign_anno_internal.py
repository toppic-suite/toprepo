"""
Annotate internal ions in an msalign file.

For a spectrum whose proteoform identification has a sequence of length n, an internal
ion is any subsequence that neither starts at the first amino acid nor ends at the last
one, i.e. residues i..j with 2 <= i <= j <= n-1. There are (n-1)(n-2)/2 such
subsequences and the mass of each internal ion is the sum of the residue masses of the
amino acids it contains.

Each deconvoluted mass in the spectrum is matched against that list within a ppm
tolerance (10 ppm by default). Masses that already carry an annotation are left exactly
as they are and are never given an internal-ion annotation, so this script only fills in
the masses that no other ion type explained.

An annotated internal ion is written with the same six trailing columns used by the
terminal-ion annotators, so an annotated peak line stays at 10 columns:

    mass  intensity  charge  ?  ion  aa_num  pos  shift  delta_da  delta_ppm

with ion = "I<first>-<last>" (1-based, inclusive), aa_num = the number of residues in the
subsequence, pos = its first residue, and shift = 0. Note that the leading "I" is not a
terminal ion type, so code that switches on the first character of the ion label (e.g.
filter_msalign_ml.py) will not recognize these lines.

When more than one internal ion falls within the tolerance the closest by mass is
reported. Internal ions are heavily degenerate -- any two subsequences with the same
residue composition have exactly the same mass -- so the reported subsequence is often
only one of several equally good candidates, and which one is reported is arbitrary. The
number of candidate masses also grows quadratically with sequence length, so a match to
an internal ion is much weaker evidence than a match to a terminal ion.

Internal ion masses are built from unmodified residue masses. Fixed and unexpected
modifications recorded in the spectrum metadata are NOT applied, so internal ions
spanning a modified residue will not match.

--shuffle turns the run into a decoy control: each sequence is randomly shuffled before
its internal ions are built, so the internal ion masses no longer correspond to any real
subsequence of the proteoform while the spectra, the tolerance and the set of masses
offered for matching all stay the same. Comparing the number of matches with and without
--shuffle estimates how many matches are coincidence. Decoy annotations are labelled
"DI<first>-<last>" instead of "I<first>-<last>" so that a decoy output file can never be
mistaken for a real annotation. Use --seed to reproduce a decoy run.

Note that an internal ion mass depends only on the composition of the subsequence, and
shuffling preserves the composition of the whole sequence, so a shuffled sequence still
produces many of the same masses. The decoy match count is therefore an over-estimate of
the false match rate rather than a strict null.

Usage:
    python3 msalign_anno_internal.py --msalign <in.msalign> --out <out.msalign> [--error_tol 10]
"""

import argparse
import random
import time

import numpy as np

from process.msalign import msalign_reader
from process.msalign import msalign_writer

AMINO_ACID_MASSES = {
    "A":71.03711,
    "R":156.10111,
    "N":114.04293,
    "D":115.02694,
    "C":103.00918,
    "E":129.04259,
    "Q":128.05858,
    "G":57.02146,
    "H":137.05891,
    "I":113.08406,
    "L":113.08406,
    "K":128.09496,
    "M":131.04049,
    "F":147.06841,
    "P":97.05276,
    "S":87.03203,
    "T":101.04768,
    "W":186.07931,
    "Y":163.06333,
    "V":99.06841,
    "U":150.95363,
    "O":132.08988
}

# A peak line holds this many columns before any annotation is appended.
BASE_COL_NUM = 4


def build_internal_mass_table(seq):
    """Return (masses, starts, ends) for every internal ion of seq, sorted by mass.

    Residues are numbered from 1. An internal ion covers residues start..end with
    2 <= start <= end <= len(seq)-1, giving (n-1)(n-2)/2 ions. Unknown amino acids
    contribute mass 0 and are reported by the caller through count_unknown_residues().
    """
    n = len(seq)
    if n < 3:
        empty = np.empty(0, dtype=float)
        return empty, np.empty(0, dtype=int), np.empty(0, dtype=int)

    residue_masses = np.array([AMINO_ACID_MASSES.get(aa, 0.0) for aa in seq], dtype=float)
    # prefix[k] is the summed mass of the first k residues, so residues i..j sum to
    # prefix[j] - prefix[i-1].
    prefix = np.concatenate(([0.0], np.cumsum(residue_masses)))

    mass_blocks = []
    start_blocks = []
    end_blocks = []
    for start in range(2, n):
        ends = np.arange(start, n)
        mass_blocks.append(prefix[ends] - prefix[start - 1])
        start_blocks.append(np.full(ends.shape, start))
        end_blocks.append(ends)

    masses = np.concatenate(mass_blocks)
    starts = np.concatenate(start_blocks)
    ends = np.concatenate(end_blocks)

    order = np.argsort(masses)
    return masses[order], starts[order], ends[order]


def count_unknown_residues(seq):
    return sum(1 for aa in seq if aa not in AMINO_ACID_MASSES)


def find_closest_internal_ion(exp_mass, masses, ppm_tol):
    """Return the index of the internal ion closest to exp_mass, or -1 if none is within
    ppm_tol. masses must be sorted, so only the two neighbours of the insertion point can
    be the closest."""
    idx = np.searchsorted(masses, exp_mass)

    best_idx = -1
    best_diff = None
    for cand in (idx - 1, idx):
        if cand < 0 or cand >= len(masses):
            continue
        theo_mass = masses[cand]
        if theo_mass <= 0:
            continue
        diff = abs(exp_mass - theo_mass)
        if best_diff is None or diff < best_diff:
            best_diff = diff
            best_idx = cand

    if best_idx < 0:
        return -1
    theo_mass = masses[best_idx]
    if abs((exp_mass - theo_mass) / theo_mass) * 1e6 > ppm_tol:
        return -1
    return best_idx


def annot_one_spectrum(spectrum, ppm_tol=10.0, rng=None):
    """Annotate the unannotated masses of one spectrum. Returns
    (spectrum, mass_num, annotated_num, kept_num, internal_ion_num).

    Passing an rng turns this into a decoy run: the sequence is shuffled before its
    internal ions are built and the annotations are labelled "DI..." instead of "I...".
    """
    peak_lines = spectrum.get("peak_lines", [])
    seq = spectrum["meta"].get("DATABASE_SEQUENCE", None)

    mass_num = 0
    annotated_num = 0
    kept_num = 0

    # Without an identification there is no sequence to build internal ions from, so the
    # peaks are passed through untouched but still counted.
    if seq in (None, ""):
        for line in peak_lines:
            parts = line.split()
            if len(parts) >= BASE_COL_NUM:
                mass_num += 1
                if len(parts) > BASE_COL_NUM:
                    kept_num += 1
        return spectrum, mass_num, annotated_num, kept_num, 0

    # Decoy control: shuffling keeps the length and composition of the sequence, so the
    # number of internal ions and the overall mass range are unchanged, but the masses
    # no longer correspond to real subsequences of the proteoform.
    ion_prefix = "I"
    if rng is not None:
        seq = "".join(rng.sample(seq, len(seq)))
        ion_prefix = "DI"

    masses, starts, ends = build_internal_mass_table(seq)
    internal_ion_num = len(masses)

    out_lines = []
    for line in peak_lines:
        stripped = line.strip()
        parts = stripped.split()

        if len(parts) < BASE_COL_NUM:
            out_lines.append(stripped)
            continue

        mass_num += 1

        # An existing annotation is kept as it is and never replaced.
        if len(parts) > BASE_COL_NUM:
            kept_num += 1
            out_lines.append(stripped)
            continue

        if len(masses) == 0:
            out_lines.append(stripped)
            continue

        exp_mass = float(parts[0])
        match_idx = find_closest_internal_ion(exp_mass, masses, ppm_tol)
        if match_idx < 0:
            out_lines.append(stripped)
            continue

        theo_mass = masses[match_idx]
        start = int(starts[match_idx])
        end = int(ends[match_idx])
        delta_da = exp_mass - theo_mass
        delta_ppm = (delta_da / theo_mass) * 1e6

        annotation = [
            f"{ion_prefix}{start}-{end}",
            str(end - start + 1),
            str(start),
            "0",
            f"{delta_da:.4f}",
            f"{delta_ppm:.2f}"
        ]
        out_lines.append("\t".join(parts + annotation))
        annotated_num += 1

    spectrum["peak_lines"] = out_lines
    return spectrum, mass_num, annotated_num, kept_num, internal_ion_num


def annot_msalign(input_msalign, output_file, ppm_tol=10.0, shuffle=False, seed=0):
    ms_reader = msalign_reader.MsalignReader(input_msalign)
    ms_writer = msalign_writer.MsalignWriter(output_file)

    # A single generator seeded once, so the whole decoy run is reproducible from --seed.
    rng = random.Random(seed) if shuffle else None
    if shuffle:
        print(f"DECOY RUN: sequences are shuffled before annotation (seed {seed}). "
              f"Annotations are labelled DI<first>-<last> and are not real identifications.")

    spectrum_num = 0
    no_seq_num = 0
    unknown_residue_num = 0
    total_mass_num = 0
    total_annotated_num = 0
    total_kept_num = 0
    total_internal_ion_num = 0

    start_time = time.time()
    for spectrum in ms_reader.readmsalign_iter():
        seq = spectrum["meta"].get("DATABASE_SEQUENCE", None)
        if seq in (None, ""):
            no_seq_num += 1
        else:
            unknown_residue_num += count_unknown_residues(seq)

        spectrum, mass_num, annotated_num, kept_num, internal_ion_num = annot_one_spectrum(
            spectrum, ppm_tol=ppm_tol, rng=rng)
        ms_writer.write(spectrum)

        spectrum_num += 1
        total_mass_num += mass_num
        total_annotated_num += annotated_num
        total_kept_num += kept_num
        total_internal_ion_num += internal_ion_num

        if spectrum_num % 1000 == 0:
            elapsed_time = time.time() - start_time
            print(f"\rAnnotated {spectrum_num} spectra. Elapsed time: {elapsed_time:.2f} seconds.",
                  end="", flush=True)
    ms_writer.close()

    elapsed_time = time.time() - start_time
    print(f"\rAnnotated {spectrum_num} spectra. Elapsed time: {elapsed_time:.2f} seconds.")
    print("\n========== Internal Ion Annotation Summary ==========")
    if shuffle:
        print(f"Mode                                   : DECOY (shuffled, seed {seed})")
    else:
        print(f"Mode                                   : target")
    print(f"Spectra                                : {spectrum_num}")
    print(f"Spectra without an identified sequence : {no_seq_num}")
    print(f"Internal ions generated for matching   : {total_internal_ion_num}")
    print(f"Total masses                           : {total_mass_num}")
    print(f"Masses with an existing annotation     : {total_kept_num}")
    print(f"Masses annotated as internal ions      : {total_annotated_num}")
    if total_mass_num > 0:
        print(f"Internal ions as a share of all masses : "
              f"{total_annotated_num / total_mass_num * 100:.2f}%")
    if unknown_residue_num > 0:
        print(f"Warning: {unknown_residue_num} unrecognized amino acids were given mass 0.")
    print(f"Output written to                      : {output_file}")
    print("====================================================")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Annotate internal ions of all spectra in an msalign file.")
    parser.add_argument(
        "--msalign", required=True, type=str, help="Input msalign filename")
    parser.add_argument(
        "--out", required=True, type=str, help="Output annotated msalign filename")
    parser.add_argument(
        "--error_tol", required=False, type=float, default=10.0,
        help="Error tolerance in ppm (default: 10.0)")
    parser.add_argument(
        "--shuffle", required=False, action='store_true',
        help="Decoy control: shuffle each proteoform sequence before building its "
             "internal ions, and label the annotations DI<first>-<last>")
    parser.add_argument(
        "--seed", required=False, type=int, default=0,
        help="Random seed used by --shuffle (default: 0)")

    args = parser.parse_args()
    annot_msalign(args.msalign, args.out, ppm_tol=args.error_tol,
                  shuffle=args.shuffle, seed=args.seed)
