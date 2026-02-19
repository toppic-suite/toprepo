import os
import argparse
import mgf_anno_util
from multiprocessing import Pool, cpu_count
import time


def annotation_batch_processing(theo_file, msalign_dir, mgf_dir, out_dir, num_workers=None):
    """
    Parameters:
        theo_file [str]: theoretical envelop file "theo_patt.txt".
        msalign_dir [str]: annotated msalign file folder.
        mgf_dir [str]: mgf folder stores mgf files to be annotated.
        out_dir[str]: output directory
        num_workers [int]: number of cpu threads will be used, default = max(cpu)-1 
    """
    # start_time = time.time()    
    # Set workers
    num_workers = num_workers or max(cpu_count() - 1, 1)
    count = 0
    os.makedirs(out_dir, exist_ok=True)

    msalign_files = [
        f for f in os.listdir(msalign_dir) if f.endswith("_ms2_annot.msalign")
    ]

    mgf_files = [
        f for f in os.listdir(mgf_dir) if f.endswith("_ms2.mgf")
    ]
    # find matched file name
    mgf_file_dict = {os.path.splitext(f)[0]: f for f in mgf_files}
    msalign_file_dict = {os.path.splitext(f)[0].replace('_annot', ''): f for f in msalign_files}
    # Match files
    common_keys = set(mgf_file_dict.keys()) & set(msalign_file_dict.keys())    
    ppm_tol = 20
    for key in common_keys:
        mgf_path = os.path.join(mgf_dir, mgf_file_dict[key])
        msalign_path = os.path.join(msalign_dir, msalign_file_dict[key])
        output_file = mgf_file_dict[key] + "_annot.mgf"
        output_path = os.path.join(out_dir, output_file)
        filenames = [msalign_path, mgf_path]
        if all(os.path.isfile(f) for f in filenames):
            # get ms2 data
            print(f"Annotation for file: {os.path.basename(mgf_path)}")  
            ms2_df = mgf_anno_util.load_msalign_data(msalign_path)
            #print(f"Loaded msalign data in {time.time() - time_start:.2f} seconds")
            time_start = time.time()
            mgf_df = mgf_anno_util.load_mgf_data(mgf_path)
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
            with Pool(num_workers) as pool, open(output_path, "w", encoding="utf-8", buffering=1024*1024*1024) as out:
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
            count += 1
                    
            print("\n========== Annotation Summary ==========")
            print(f"Processed file number   : {count}")
            print(f"Spectra written to file : {annotated_block_count}")
            print(f"Time elapsed            : {elapsed:.2f} seconds ({mins:.2f} min)")
            print("========================================\n")  
            
        else:
            print(f"One of mgf or msalign files missing for {key}")
                                  

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
        "--mgf_dir", "-m",
        required=True,
        help="Directory containing MGF files (*_ms2.mgf)"
    )
    parser.add_argument(
        "--msalign_dir", "-s",
        required=True,
        help="Directory containing annotated msalign files (*_ms2_annot.msalign)"
    )
    parser.add_argument(
        "--out_dir", "-o",
        required=True,
        help="Output directory for annotated MGF files"
    )
    parser.add_argument(
        "--num_workers", "-n",
        type=int,
        default=None,
        help="Number of CPU workers (default: max CPUs - 1)"
    )

    args = parser.parse_args()
    annotation_batch_processing(
        args.theo_file,
        args.msalign_dir,
        args.mgf_dir,
        args.out_dir,
        args.num_workers
    )

