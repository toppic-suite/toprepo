#!/bin/bash

DIR="$(dirname "$(realpath "$0")")"

SINGLE_DIR="${1:-$DIR/03_single_prsm_remove_parameters}"
PRSM_DIR="${2:-$DIR/03_prsm_remove_parameters}"
OUT_DIR="${3:-$DIR/05_prsm_info}"

mkdir -p "$OUT_DIR"

for f in "$SINGLE_DIR"/*_ms2_toppic_prsm_single_replaced.tsv; do
    basename="$(basename "$f")"
    base="${basename%_ms2_toppic_prsm_single_replaced.tsv}"
    project_id="${basename%%_*}"
    prsm="$PRSM_DIR/${base}_ms2_toppic_prsm_replaced.tsv"
    out="$OUT_DIR/${base}_toppic_info.tsv"
    python3 prsm_alt_proteins.py "$f" "$prsm" "$project_id" "$out"
    echo "Processed: $basename + $(basename "$prsm") $project_id -> $(basename "$out")"
done

echo "Done."
