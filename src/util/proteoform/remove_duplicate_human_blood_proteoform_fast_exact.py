"""
remove duplicate proteoform (human blood subset)

Fast, output-identical version of remove_duplicate_human_blood_proteoform.py. Only the
human-blood datasets are deduplicated (grouped by TOPPIC_protein_accession); all other
rows are passed through unchanged. The per-group dedup is delegated to
remove_duplicate_proteoform_fast_exact.deduplicate_group, which reproduces the original
deduplicate_group byte-for-byte while avoiding pandas per-row access and pruning
comparisons with a precursor-mass window. main() is unchanged from the original, so the
output is identical.
"""

import pandas as pd
import sys
import remove_duplicate_proteoform_fast_exact


HUMAN_BLOOD_DATASETS = list(set([
    'PXD026123','PXD026124','PXD026126','PXD026127','PXD026128','PXD026129',
    'PXD026130','PXD026131','PXD026132','PXD026133','PXD026134','PXD026135',
    'PXD026136','PXD026137','PXD026138','PXD026139','PXD026140','PXD026141',
    'PXD026142','PXD026143','PXD026144','PXD026145','PXD026146','PXD026147',
    'PXD026148','PXD026149','PXD026150','PXD026151','PXD026152','PXD026153',
    'PXD026154','PXD026155','PXD026156','PXD026157','PXD026158','PXD026159',
    'PXD026160','PXD026161','PXD026162','PXD026163','PXD026164','PXD026165',
    'PXD026166','PXD026167','PXD026168','PXD026169','PXD026170','PXD026171',
    'PXD026172','PXD026173','PXD026174','PXD026175','PXD026176','PXD026177',
    'PXD026178'
]))


def process_human_blood(df):
    human_blood_df = df[df['DATASET_id'].isin(HUMAN_BLOOD_DATASETS)].copy()
    rest_df = df[~df['DATASET_id'].isin(HUMAN_BLOOD_DATASETS)].copy()

    return human_blood_df, rest_df


def main(input_file, output_file):
    df = pd.read_csv(
        input_file,
        sep="\t",
        low_memory=False
    )
    print(f"{len(df)}")
    human_df, rest_df = process_human_blood(df)
    print(f"Human blood rows: {len(human_df)}")
    print(f"Rest rows: {len(rest_df)}")

    # --- Deduplicate human blood ---
    human_df['MSALIGN_precursor_mass'] = human_df['MSALIGN_precursor_monoisotopic_mass'].apply(lambda x: x.split(':')[0])
    groups = []
    for acc, group in human_df.groupby(
        ["TOPPIC_protein_accession"]
    ):
        groups.append(remove_duplicate_proteoform_fast_exact.deduplicate_group(group))

    dedup_human_df = pd.concat(groups, ignore_index=True)
    print(f"Producing {len(dedup_human_df)} deduplicated rows")

    dedup_human_df.drop(columns=["MSALIGN_precursor_mass"], axis=1, inplace=True)
    final_df = pd.concat([rest_df, dedup_human_df], ignore_index=True)

    num_ptms = final_df['TOPPIC_proteoform_id'].notna().sum()
    final_df.to_csv(
        output_file,
        sep="\t",
        index=False)

    print(f"Total spectra (rows): {len(final_df)}")
    print(f"Rows with proteoform ID: {num_ptms}")
    print("Done.")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    main(input_file, output_file)
