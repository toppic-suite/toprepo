#!/bin/bash
#
# Run remove_header_rows.sh over a folder of TopPIC TSVs, or over every dataset subfolder of
# one, stripping the "********** Parameters **********" preamble from each file.
#
# remove_header_rows.sh globs its own directory and calls "python3 remove_params.py" from the
# current working directory, so it only ever processes one flat folder and only ever writes
# <name>_replaced.tsv beside its input. This wrapper drives it across a folder tree and can
# redirect the results into a separate output folder.
#
# Usage: remove_folder_header_rows.sh [-o OUTPUT_DIR] [-j JOBS] [-f] [-n] TARGET_DIR
#
#   TARGET_DIR  folder of *.tsv, or a folder of <dataset_id>/*.tsv subfolders
#   -o          write results into OUTPUT_DIR instead of beside the inputs; TARGET_DIR is left
#               untouched and its subfolder structure is reproduced under OUTPUT_DIR
#   -j          number of folders processed in parallel (default: 8)
#   -f          process folders that already have output
#   -n          dry run: list the folders that would be processed
#
# Output filenames are whatever remove_params.py wrote — <name>.tsv becomes <name>_replaced.tsv
# — in both modes. This wrapper never renames them.
#
# Without -o, a folder that already holds *_replaced.tsv is skipped: remove_header_rows.sh globs
# *.tsv, so a second pass would turn that output into *_replaced_replaced.tsv. Inputs already
# named *_replaced.tsv are never fed to the script in -o mode either.
#
# In -o mode the inputs are hard linked into the output folder as scratch, the script is run
# there, and the links are removed once it finishes — so no input data is copied or modified.

set -uo pipefail

DIR="$(dirname "$(realpath "$0")")"
SH="$DIR/remove_header_rows.sh"
PY="$DIR/../src/util/tsv/remove_params.py"

output_dir=""
jobs=8
force=0
dry_run=0

while getopts ":o:j:fnh" opt; do
    case "$opt" in
        o) output_dir="$OPTARG" ;;
        j) jobs="$OPTARG" ;;
        f) force=1 ;;
        n) dry_run=1 ;;
        h) sed -n '2,28p' "$0"; exit 0 ;;
        \?) echo "Unknown option: -$OPTARG" >&2; exit 1 ;;
        :) echo "Option -$OPTARG requires an argument." >&2; exit 1 ;;
    esac
done
shift $((OPTIND - 1))

target="${1:-}"
[ -n "$target" ] || { echo "Usage: $(basename "$0") [-o OUTPUT_DIR] [-j JOBS] [-f] [-n] TARGET_DIR" >&2; exit 1; }
[ -d "$target" ] || { echo "No such directory: $target" >&2; exit 1; }
[ -f "$SH" ] || { echo "Missing $SH" >&2; exit 1; }
[ -f "$PY" ] || { echo "Missing $PY" >&2; exit 1; }

target="$(realpath "$target")"
SH="$(realpath "$SH")"
PY="$(realpath "$PY")"
if [ -n "$output_dir" ]; then
    [ "$dry_run" -eq 0 ] && mkdir -p "$output_dir"
    [ -d "$output_dir" ] && output_dir="$(realpath "$output_dir")"
fi
export SH PY force output_dir target

# Strip the preamble from every *.tsv in $1, leaving the results in $1 (in-place mode).
process_in_place() {
    dir="$1"
    name="$2"
    if [ "$force" -eq 0 ] && ls "$dir"/*_replaced.tsv >/dev/null 2>&1; then
        echo "SKIP (already processed): $name"
        return 0
    fi
    cp "$SH" "$PY" "$dir"/ || { echo "FAIL (copy): $name"; return 1; }
    ( cd "$dir" && bash "$(basename "$SH")" >/dev/null 2>&1 ) || echo "FAIL (run): $name"
    rm -f "$dir/$(basename "$SH")" "$dir/$(basename "$PY")"
    echo "OK $name ($(ls "$dir"/*_replaced.tsv 2>/dev/null | wc -l) files)"
}

# Strip the preamble from every *.tsv in $1, leaving the results in $3 and $1 untouched.
process_to_output() {
    dir="$1"
    name="$2"
    dest="$3"
    if [ "$force" -eq 0 ] && [ -d "$dest" ] && [ -n "$(ls -A "$dest" 2>/dev/null)" ]; then
        echo "SKIP (output exists): $name"
        return 0
    fi
    mkdir -p "$dest" || { echo "FAIL (mkdir): $name"; return 1; }

    # Stage the inputs inside the output folder so the script's own-directory glob sees them.
    # Hard links keep this free; fall back to copying across filesystems.
    staged=0
    for f in "$dir"/*.tsv; do
        [ -e "$f" ] || continue
        case "$(basename "$f")" in *_replaced.tsv) continue ;; esac
        ln "$f" "$dest/$(basename "$f")" 2>/dev/null || cp -p "$f" "$dest/$(basename "$f")"
        staged=$((staged + 1))
    done
    if [ "$staged" -eq 0 ]; then
        echo "SKIP (no TSV): $name"
        rmdir "$dest" 2>/dev/null
        return 0
    fi

    cp "$SH" "$PY" "$dest"/ || { echo "FAIL (copy): $name"; return 1; }
    ( cd "$dest" && bash "$(basename "$SH")" >/dev/null 2>&1 ) || echo "FAIL (run): $name"
    rm -f "$dest/$(basename "$SH")" "$dest/$(basename "$PY")"

    # Drop the staged inputs, keeping only what the script produced.
    for f in "$dest"/*.tsv; do
        case "$(basename "$f")" in *_replaced.tsv) continue ;; esac
        rm -f "$f"
    done

    echo "OK $name ($(ls "$dest"/*_replaced.tsv 2>/dev/null | wc -l) files)"
}

process_one() {
    dir="$1"
    name="$(basename "$dir")"
    if ! ls "$dir"/*.tsv >/dev/null 2>&1; then
        echo "SKIP (no TSV): $name"
        return 0
    fi
    if [ -z "$output_dir" ]; then
        process_in_place "$dir" "$name"
    elif [ "$dir" = "$target" ]; then
        process_to_output "$dir" "$name" "$output_dir"
    else
        process_to_output "$dir" "$name" "$output_dir/$name"
    fi
}
export -f process_in_place process_to_output process_one

# A folder of TSVs with no subfolders is processed as a single unit.
if [ -z "$(find "$target" -mindepth 1 -maxdepth 1 -type d -print -quit)" ]; then
    folders="$target"
else
    folders="$(find "$target" -mindepth 1 -maxdepth 1 -type d | sort)"
fi

count=$(printf '%s\n' "$folders" | grep -c .)
echo "Processing $count folder(s) under $target with $jobs parallel job(s)"
[ -n "$output_dir" ] && echo "Writing results to $output_dir"

if [ "$dry_run" -eq 1 ]; then
    printf '  %s\n' $folders
    echo "Dry run: nothing written."
    exit 0
fi

printf '%s\n' "$folders" | xargs -d '\n' -n1 -P "$jobs" bash -c 'process_one "$0"'

echo "Done."
