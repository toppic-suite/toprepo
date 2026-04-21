"""
remove duplicate proteoform
"""

import pandas as pd
import sys
import remove_duplicate_proteoform_fast


def main(input_file, output_file):
    df = pd.read_csv(
        input_file,
        sep="\t",
        low_memory=False
    )
    print(f"Total input rows: {len(df)}")
    df['MSALIGN_precursor_mass'] = df['MSALIGN_precursor_monoisotopic_mass'].apply(lambda x: x.split(':')[0])
    total_rows = len(df)
    rows_processed = 0
    last_printed = 0
    groups = []
    for (proj, acc), group in df.groupby(
        ["PROJECT_id","TOPPIC_protein_accession"]
    ):
        groups.append(remove_duplicate_proteoform_fast.deduplicate_group(group))
        rows_processed += len(group)
        if rows_processed - last_printed >= 1000:
            print(f"Progress: {rows_processed}/{total_rows} rows ({rows_processed/total_rows*100:.1f}%)")
            last_printed = rows_processed
    
    deduped_chunk = pd.concat(groups, ignore_index=True)
    print(f"Producing {len(deduped_chunk)} deduplicated rows")
    
    deduped_chunk.drop(columns=["MSALIGN_precursor_mass"], axis=1, inplace=True)
    num_ptms = deduped_chunk['TOPPIC_proteoform_id'].notna().sum()
    
    deduped_chunk.to_csv(
        output_file,
        sep="\t",
        index=False)
    
    
    print(f"Total spectra (rows): {len(df)}, total remove duplicates: {len(deduped_chunk)}")
    print(f"Rows with proteoform ID: {num_ptms}")
    print("Done.")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    main(input_file, output_file)
