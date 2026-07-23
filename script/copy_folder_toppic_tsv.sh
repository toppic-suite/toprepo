#!/bin/bash
#
# Collect one kind of TopPIC result file out of an extracted TopPIC output tree into a new
# folder that mirrors its dataset structure. Used to build 01_prsm and 02_proteoform from
# 00_toprepo_toppic_output.
#
# Input tree (as produced by extract_folder_toppic_zip.sh):
#     00_toprepo_toppic_output/<dataset_id>/<run>_ms2/<run>_ms2_toppic_prsm.tsv
# Output tree (one flat folder of TSVs per dataset, the run name is already in the filename):
#     01_prsm/<dataset_id>/<run>_ms2_toppic_prsm.tsv
#
# Usage: copy_folder_toppic_tsv.sh [-s SUFFIX] [-l] [-n] INPUT_DIR OUTPUT_DIR
#
#   -s   filename suffix to collect (default: prsm.tsv)
#   -l   hard link instead of copying — instant and free, but the two trees then share
#        storage, so editing one file changes both
#   -n   dry run: report how many files would be copied per dataset
#
# Note that the default suffix prsm.tsv matches <run>_ms2_toppic_prsm.tsv but NOT
# <run>_ms2_toppic_prsm_single.tsv; pass -s prsm_single.tsv for those. Likewise
# proteoform.tsv versus proteoform_single.tsv.

set -uo pipefail

suffix="prsm.tsv"
link=0
dry_run=0

while getopts ":s:lnh" opt; do
    case "$opt" in
        s) suffix="$OPTARG" ;;
        l) link=1 ;;
        n) dry_run=1 ;;
        h) sed -n '2,24p' "$0"; exit 0 ;;
        \?) echo "Unknown option: -$OPTARG" >&2; exit 1 ;;
        :) echo "Option -$OPTARG requires an argument." >&2; exit 1 ;;
    esac
done
shift $((OPTIND - 1))

input_dir="${1:-}"
output_dir="${2:-}"
if [ -z "$input_dir" ] || [ -z "$output_dir" ]; then
    echo "Usage: $(basename "$0") [-s SUFFIX] [-l] [-n] INPUT_DIR OUTPUT_DIR" >&2
    exit 1
fi
[ -d "$input_dir" ] || { echo "No such directory: $input_dir" >&2; exit 1; }
input_dir="$(realpath "$input_dir")"

copied=0
missing=0
for dataset_path in "$input_dir"/*/; do
    [ -d "$dataset_path" ] || continue
    dataset="$(basename "$dataset_path")"
    files="$(find "$dataset_path" -type f -name "*$suffix" | sort)"
    n=$(printf '%s\n' "$files" | grep -c .)
    if [ "$n" -eq 0 ]; then
        echo "  $dataset: no *$suffix found"
        missing=$((missing + 1))
        continue
    fi
    echo "  $dataset: $n file(s)"
    copied=$((copied + n))
    [ "$dry_run" -eq 1 ] && continue
    mkdir -p "$output_dir/$dataset"
    while IFS= read -r f; do
        if [ "$link" -eq 1 ]; then
            ln -f "$f" "$output_dir/$dataset/$(basename "$f")"
        else
            cp -p "$f" "$output_dir/$dataset/"
        fi
    done <<< "$files"
done

if [ "$dry_run" -eq 1 ]; then
    echo "Dry run: $copied file(s) would be collected into $output_dir ($missing dataset(s) with no match)."
else
    echo "Collected $copied file(s) matching *$suffix into $output_dir ($missing dataset(s) with no match)."
fi
