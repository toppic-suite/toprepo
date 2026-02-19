import os
import argparse
import mgf_anno_util
from multiprocessing import Pool, cpu_count
import time


def annotation_processing(theo_file, msalign_filename, mgf_filename, out_filename, num_workers=None):
    # start_time = time.time()    
    # Set workers
    num_workers = num_workers or max(cpu_count() - 1, 1)
    ppm_tol = 20

    # get ms2 data
    print(f"Annotation for file: {os.path.basename(mgf_filename)}")  
    ms2_df = mgf_anno_util.load_msalign_data(msalign_filename)
    #print(f"Loaded msalign data in {time.time() - time_start:.2f} seconds")
    time_start = time.time()
    mgf_df = mgf_anno_util.load_mgf_data(mgf_filename)
    print(f"Loaded mgf data in {time.time() - time_start:.2f} seconds")
    form_df = mgf_anno_util.combined_msalign_mgf(ms2_df, mgf_df)
    #print(f"Combined msalign and mgf data in {time.time() - time_start:.2f} seconds")   
    # Multiprocessing
    print(f"Processing {len(form_df)} spectra with {num_workers} workers...")

    tasks = [
        (row._asdict(), theo_file, ppm_tol)
        for row in form_df.itertuples(index=False)
    ]

    annotated_block_count = 0
    start_time = time.time()    
    with Pool(num_workers) as pool, open(out_filename, "w", encoding="utf-8", buffering=1024*1024*1024) as out:
        for result in pool.imap(mgf_anno_util.process_one_spectrum, tasks, chunksize=50):
            meta = result["meta"]
            peaks = result["peaks"]
            meta_lines = result["meta_lines"]
            out.write("BEGIN IONS\n")
            for k, v in meta.items():
                out.write(f"{k}={v}\n")
            for meta_line in meta_lines:
                out.write(meta_line + "\n")
            for line in peaks:
                out.write(line + "\n")
            out.write("END IONS\n\n")
            
            annotated_block_count += 1            
            if annotated_block_count % 100 == 0:
                print(f"\rAnnotated {annotated_block_count} spectra...", end='', flush=True)
                    
    end_time = time.time()
    elapsed = end_time - start_time
    mins = elapsed / 60
                    
    print("\n========== Annotation Summary ==========")
    print(f"Spectra written to file : {annotated_block_count}")
    print(f"Time elapsed            : {elapsed:.2f} seconds ({mins:.2f} min)")
    print("========================================\n")  
            
                                  

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Annotate MGF files with theoretical fragment patterns"
    )
    parser.add_argument(
        "--theo_file", "-t",
        required=True,
        help="Path to the theoretical envelope file (theo_patt.txt)"
    )
    parser.add_argument(
        "--mgf_file", "-m",
        required=True,
        help="Path to a single MGF file (*_ms2.mgf)"
    )
    parser.add_argument(
        "--msalign_file", "-s",
        required=True,
        help="Path to a single annotated msalign file (*_ms2_annot.msalign)"
    )
    parser.add_argument(
        "--out", "-o",
        required=True,
        help="Output path for annotated MGF file"
    )
    parser.add_argument(
        "--num_workers", "-n",
        type=int,
        default=None,
        help="Number of CPU workers (default: max CPUs - 1)"
    )

    args = parser.parse_args()

    annotation_processing(
        args.theo_file,
        args.msalign_file,
        args.mgf_file,
        args.out,
        args.num_workers
    )

