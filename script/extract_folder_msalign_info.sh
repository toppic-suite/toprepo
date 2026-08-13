#!/bin/bash

DIR="$(dirname "$(realpath "$0")")"

for f in "$DIR"/*ms2.msalign; do
    basename="$(basename "$f")"
    base="${basename%_ms2.msalign}"
    project_id="${basename%%_*}"
    out="${base}_msalign_info.tsv"
    python3 extract_msalign_info.py $project_id $f $out
    echo "Processed: $(basename "$f") $project_id -> $(basename "$out")"
done

echo "Done."
