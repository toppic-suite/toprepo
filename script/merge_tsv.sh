#!/usr/bin/env bash

if [ $# -ne 1 ]; then
    echo "Usage: $0 <output_file>"
    exit 1
fi

output="$1"
first=1

for f in *.tsv; do
    [ -e "$f" ] || { echo "No TSV files found."; exit 1; }
    if [ $first -eq 1 ]; then
        cat "$f"
        first=0
    else
        tail -n +2 "$f"
    fi
done > "$output"

echo "Merged into $output"
