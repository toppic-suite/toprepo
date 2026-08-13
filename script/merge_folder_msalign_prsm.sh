#!/bin/bash

DIR="$(dirname "$(realpath "$0")")"

shopt -s nullglob
for f in "${DIR}"/*ms2_preprocess.msalign; do
    basename="$(basename "$f")"
    base="${basename%_ms2_preprocess.msalign}"
    combined_info=${base}_combined_info.tsv
    out="${base}_ms2_prsm.msalign"
    python3 merge_msalign_prsm.py --tsv "$combined_info" --msalign "$f" --out "$out"
    echo "Processed: $(basename "$f") -> $out"
done

echo "Done."
