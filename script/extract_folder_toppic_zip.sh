#!/bin/bash
#
# Extract every TopPIC output zip found under a folder (e.g. 00_toprepo_toppic_output).
#
# Each <dataset_id>/<run>_ms2_toppic_output.zip is unpacked into
# <output_dir>/<dataset_id>/<run>_ms2/, so the dataset layout of the input folder
# is preserved. Runs that already have a non-empty output directory are skipped
# unless -f is given.
#
# Usage: extract_folder_toppic_zip.sh [-o OUTPUT_DIR] [-p PATTERN] [-j JOBS] [-f] [-n] [INPUT_DIR]
#
#   INPUT_DIR   folder holding <dataset_id>/*.zip (default: current directory)
#   -o          output folder (default: INPUT_DIR)
#   -p          only extract archive members matching this pattern, e.g. '*.tsv'
#   -j          number of parallel unzip jobs (default: 4)
#   -f          re-extract runs that already have output
#   -n          dry run: list what would be extracted
#
# Extracting all 4,615 zips of TopRepo 1.2.1 needs roughly 48 GB; use -p '*.tsv'
# (~7 GB) when the XML and raw_prsm members are not needed.

set -euo pipefail

# Worker mode: invoked once per zip by xargs below.
if [ "${1:-}" = "--extract-one" ]; then
    zip="$2"
    dest="$3"
    pattern="$4"
    mkdir -p "$dest"
    if [ -n "$pattern" ]; then
        unzip -q -o "$zip" "$pattern" -d "$dest"
    else
        unzip -q -o "$zip" -d "$dest"
    fi
    echo "Extracted: $(basename "$zip")"
    exit 0
fi

output_dir=""
pattern=""
jobs=4
force=0
dry_run=0

while getopts ":o:p:j:fnh" opt; do
    case "$opt" in
        o) output_dir="$OPTARG" ;;
        p) pattern="$OPTARG" ;;
        j) jobs="$OPTARG" ;;
        f) force=1 ;;
        n) dry_run=1 ;;
        h) sed -n '2,25p' "$0"; exit 0 ;;
        \?) echo "Unknown option: -$OPTARG" >&2; exit 1 ;;
        :) echo "Option -$OPTARG requires an argument." >&2; exit 1 ;;
    esac
done
shift $((OPTIND - 1))

input_dir="${1:-$PWD}"
[ -d "$input_dir" ] || { echo "No such directory: $input_dir" >&2; exit 1; }
input_dir="$(realpath "$input_dir")"
output_dir="${output_dir:-$input_dir}"
mkdir -p "$output_dir"
output_dir="$(realpath "$output_dir")"

self="$(realpath "$0")"
work_list="$(mktemp)"
trap 'rm -f "$work_list"' EXIT

total=0
skipped=0
while IFS= read -r -d '' zip; do
    total=$((total + 1))
    rel="${zip#"$input_dir"/}"                        # e.g. PXD017265/run_ms2_toppic_output.zip
    dest="$output_dir/${rel%_toppic_output.zip}"      # e.g. <out>/PXD017265/run_ms2
    if [ "$force" -eq 0 ] && [ -d "$dest" ] && [ -n "$(ls -A "$dest" 2>/dev/null)" ]; then
        skipped=$((skipped + 1))
        continue
    fi
    printf '%s\0%s\0' "$zip" "$dest" >> "$work_list"
done < <(find "$input_dir" -type f -name '*.zip' -print0 | sort -z)

queued=$(( total - skipped ))
echo "Found $total zip file(s) under $input_dir"
echo "Skipping $skipped already extracted, extracting $queued into $output_dir"
[ -n "$pattern" ] && echo "Member pattern: $pattern"

if [ "$dry_run" -eq 1 ]; then
    xargs -0 -n2 printf '  %s -> %s\n' < "$work_list"
    echo "Dry run: nothing extracted."
    exit 0
fi

if [ "$queued" -gt 0 ]; then
    export TOPPIC_ZIP_SELF="$self" TOPPIC_ZIP_PATTERN="$pattern"
    xargs -0 -n2 -P "$jobs" bash -c \
        'exec "$TOPPIC_ZIP_SELF" --extract-one "$0" "$1" "$TOPPIC_ZIP_PATTERN"' < "$work_list"
fi

echo "Done."
