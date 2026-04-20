#!/bin/bash

DIR="$(dirname "$(realpath "$0")")"

for msalign_fname in "$DIR"/*msalign_info.tsv; do
  msalign_basename="$(basename "$msalign_fname")"
  msalign_base="${msalign_basename%_msalign_info.tsv}"
  echo $msalign_base
  feature_fname="${msalign_base}_feature_info.tsv"
  toppic_fname="${msalign_base}_toppic_info.tsv"
  out="${msalign_base}_combined_info.tsv"
  echo "python3 merge_mzml_msalign_toppic_info.py $msalign_basename $feature_fname $toppic_fname toprepo_file_info_v1.2.0.tsv $out"
  python3 merge_mzml_msalign_toppic_info.py $msalign_basename $feature_fname $toppic_fname $out
  echo "Processed: $msalign_basename -> $(basename "$out")"
done

echo "Done."
