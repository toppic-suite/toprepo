#!/bin/bash

DIR="$(dirname "$(realpath "$0")")"

for f in "$DIR"/*ms2.msalign; do
    basename="$(basename "$f")"
    base="${basename%_ms2.msalign}"
    project_id="${basename%%_*}"
    out="${base}_ms2_preprocess.msalign"
    python3 msalign_preprocess.py $f $project_id $out
    echo "Processed: $(basename "$f") $project_id -> $(basename "$out")"
done

echo "Done."
