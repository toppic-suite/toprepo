#!/bin/bash

DIR="$(dirname "$(realpath "$0")")"

for f in "$DIR"/*.feature; do
    basename="$(basename "$f")"
    base="${basename%_ms2.feature}"
    project_id="${basename%%_*}"
    out="${base}_feature_info.tsv"
    python3 extract_feature_info.py $project_id $f $out
    echo "Processed: $(basename "$f") $project_id -> $(basename "$out")"
done

echo "Done."
