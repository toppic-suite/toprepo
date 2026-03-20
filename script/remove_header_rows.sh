#!/bin/bash

DIR="$(dirname "$(realpath "$0")")"

for f in "$DIR"/*.tsv; do
    base="${f%.tsv}"
    out="${base}_replaced.tsv"
    python3 remove_params.py $f $out
    echo "Processed: $(basename "$f") -> $(basename "$out")"
done

echo "Done."
