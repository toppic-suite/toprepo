#!/bin/bash

DIR="$(dirname "$(realpath "$0")")"

for mzml_fname in "$DIR"/*mzml_info.tsv; do
  mzml_basename="$(basename "$mzml_fname")"
  mzml_base="${mzml_basename%_mzml_info.tsv}"
  #echo "$mzml_base"
  for msalign_fname in "$DIR"/"$mzml_base"*_msalign_info.tsv; do
    msalign_basename="$(basename "$msalign_fname")"
    msalign_base="${msalign_basename%_msalign_info.tsv}"
    echo $msalign_base
    feature_fname="${msalign_base}_feature_info.tsv"
    toppic_fname="${msalign_base}_toppic_info.tsv"
    out="${msalign_base}_combined_info.tsv"
    echo "python3 merge_mzml_msalign_toppic_info.py $mzml_basename $msalign_basename $feature_fname $toppic_fname $out"
    python3 merge_mzml_msalign_toppic_info.py $mzml_basename $msalign_basename $feature_fname $toppic_fname $out
    echo "Processed: $msalign_basename -> $(basename "$out")"
  done
done

echo "Done."
