#!/bin/bash

DIR="$(dirname "$(realpath "$0")")"

for f in "$DIR"/*.tsv; do
    basename="$(basename "$f")"
    base="${basename%_ms2_toppic_prsm_single_replaced.tsv}"
    project_id="${basename%%_*}"
    out="${base}_toppic_info.tsv"
    python3 prsm_preprocess.py $f $project_id --output $out
    echo "Processed: $(basename "$f") $project_id -> $(basename "$out")"
done

echo "Done."
