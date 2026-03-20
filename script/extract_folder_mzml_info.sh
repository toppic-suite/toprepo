#!/bin/bash

DIR="$(dirname "$(realpath "$0")")"

for f in "$DIR"/*.mzML; do
    basename="$(basename "$f")"
    base="${basename%.mzML}"
    project_id="${basename%%_*}"
    out="${base}_mzml_info.tsv"
    python3 extract_mzml_info.py $project_id $f $out
    echo "Processed: $(basename "$f") $project_id -> $(basename "$out")"
done

echo "Done."
