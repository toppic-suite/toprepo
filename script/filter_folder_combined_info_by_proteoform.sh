#!/bin/bash

DIR="$(dirname "$(realpath "$0")")"

for f in "$DIR"/*toppic_proteoform_info.tsv; do
    basename="$(basename "$f")"
    base="${basename%_toppic_proteoform_info.tsv}"
    info="${base}_combined_info.tsv"
    out="${base}_proteoform_combined_info.tsv"
    python3 filter_combined_info_by_proteoform.py $f $info $out
    echo "Processed: $(basename "$f") -> $(basename "$out")"
done

echo "Done."
