import json
import pandas as pd
import mgf_anno 


def load_mgf_data(mgf_file):
    """
    Read a mgf file and return a Python DataFrame containing all information in a mgf file
    """
    rows = []
    dataset_id = None
    mzml_filename = None
    scan = None
    title = None
    rtinseconds = None
    pepmass_mz = None
    charge = None

    mz_all = []
    intensity_all = []
   
    with open(mgf_file, "r", encoding="utf-8", errors="ignore", buffering=1024*1024*1024) as fh:
        for line in fh:
            line = line.strip()

            if len(line)> 0 and line[0] <= "9" and line[0] >= "0":
                # Peak line
                parts = line.split()
                if len(parts) == 2:
                    # mandatory fields
                    mz_all.append(float(parts[0]))
                    intensity_all.append(float(parts[1]))

            elif line == "BEGIN IONS":
                # reset for new spectrum
                dataset_id = None
                mzml_filename = None
                scan = None
                title = None
                rtinseconds = None
                pepmass_mz = None
                charge = None

                mz_all = []
                intensity_all = []
            elif line.find("=") != -1:
                if line.startswith("DATASET_ID="):
                    dataset_id = line.split("=", 1)[1]
                elif line.startswith("MZML_FILE_NAME="):
                    mzml_filename = line.split("=", 1)[1]

                elif line.startswith("SCAN="):
                    scan = int(line.split("=", 1)[1])
                
                elif line.startswith("TITLE="):
                    title = line.split("=", 1)[1]
                
                elif line.startswith("RTINSECONDS="):
                    rtinseconds = float(line.split("=", 1)[1])
                
                elif line.startswith("PEPMASS_MZ="):
                    pepmass_mz = float(line.split("=", 1)[1])
            
                elif line.startswith("CHARGE="):
                    charge = line.split("=", 1)[1]

            elif line == "END IONS":
                if scan is not None:
                    rows.append({
                        "dataset_id": dataset_id,
                        "mzml_file_name": mzml_filename,
                        "scan": scan,
                        "title": title,
                        "rtinseconds": rtinseconds,
                        "pepmass_mz": pepmass_mz,
                        "charge": charge,
                        "mz_array": mz_all,
                        "intensity_array": intensity_all
                    })

                # clean up explicitly
                dataset_id = None
                mzml_filename = None
                scan = None
                title = None
                rtinseconds = None
                pepmass_mz = None
                charge = None
                
    print(f"Loaded {len(rows)} spectra from {mgf_file}")
            
    mgf_df = pd.DataFrame(rows)
    return mgf_df



def load_msalign_data(msalign_file):
    """
    Read a msalign file and return a Python DataFrame containing all information in a msalign file
    """
    
    rows = []
    dataset_id = None
    mzml_filename = None
    scan = None

    mass_all = []
    intensity_all = []
    charge_all = []
    confidence_all = []
    ms2_deconv_label = []

    with open(msalign_file, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            line = line.strip()

            if line == "BEGIN IONS":
                # reset for new spectrum
                dataset_id = None
                mzml_filename = None
                scan = None

                mass_all = []
                intensity_all = []
                charge_all = []
                confidence_all = []
                ms2_deconv_label = []
                meta_line_list = []

            elif line.startswith("DATASET_ID="):
                dataset_id = line.split("=", 1)[1]

            elif line.startswith("MZML_FILE_NAME="):
                mzml_filename = line.split("=", 1)[1]

            elif line.startswith("MS2_SCAN="):
                scan = int(line.split("=", 1)[1])
            elif line.find("=") != -1 and not line.startswith("MS2_RETENTION_TIME="):
                meta_line_list.append(line)
            elif line == "END IONS":
                if scan is not None:
                    #print("label for scan {}: {}".format(scan, ms2_deconv_label))
                    rows.append({
                        "dataset_id": dataset_id,
                        "mzml_file_name": mzml_filename,
                        "scan": scan,
                        "meta_lines": meta_line_list,
                        "mass_all": mass_all,
                        "intensity_all": intensity_all,
                        "charge_all": charge_all,
                        "confidence_all": confidence_all,
                        "ms2_deconv_label": ms2_deconv_label
                    })

                # clean up explicitly
                dataset_id = None
                mzml_filename = None
                scan = None
            elif line == "":
                continue
            elif line.find("=") == -1:
                # Peak line
                parts = line.split()
                # mandatory fields
                mass_all.append(float(parts[0]))
                intensity_all.append(float(parts[1]))
                charge_all.append(int(parts[2]))
                confidence_all.append(float(parts[3]))
                if len(parts) > 4:
                    ms2_deconv_label.append(" ".join(parts[4:]))
                else:
                    ms2_deconv_label.append("")
                   
    msalign_df = pd.DataFrame(rows)
    
    return msalign_df


def combined_msalign_mgf(ms2_df, mgf_df):
    """
    merge MS2 data from msalign and centriod data from mgf into a Python DataFrame 
    Return a DataFrame with an additional column 'ms2_deconv_label' containing deconvoluted fragment annotations.
    """
    #print(ms2_df.columns.tolist())
    #print(mgf_df.columns.tolist())
    form_df = ms2_df.merge(mgf_df,
                          on=['dataset_id','mzml_file_name', 'scan'], 
                          how='left')
    form_df['ms2_deconv_label'] = form_df['ms2_deconv_label'].apply(json.dumps)
    return form_df
    

def process_one_spectrum(args):
    """
    Worker function: process ONE spectrum (one row of form_df)
    """
    # row, theo_file, ppm_tol = args
    # Build a single-row DataFrame
    row_dict, theo_file, ppm_tol = args
    form_df_one = pd.DataFrame([row_dict])

    # Run annotation
    labels = mgf_anno.get_ms2_centroid_label(theo_file, form_df_one, ppm_tol)
    labels = mgf_anno.filter_ms2_centroid_labels(labels)

    # Extract result
    ms2_centroid_label = labels[0]

    # Build MGF output
    meta = {
        "DATASET_ID": row_dict["dataset_id"],
        "MZML_FILE_NAME": row_dict["mzml_file_name"],
        "SCAN": row_dict["scan"]
    }
    if row_dict["title"] is not None:
        meta["TITLE"] = row_dict["title"]
    if row_dict["rtinseconds"] is not None:
        meta["RTINSECONDS"] = row_dict["rtinseconds"]
    if row_dict["pepmass_mz"] is not None:
        meta["PEPMASS_MZ"] = row_dict["pepmass_mz"]
    if row_dict["charge"] is not None:
        meta["CHARGE"] = row_dict["charge"]
    
    msalign_meta_lines = row_dict["meta_lines"]

    peaks_out = []
    for mz, inten, annot in zip(row_dict["mz_array"], row_dict["intensity_array"], ms2_centroid_label):
        if annot == "" or annot is None:
            peaks_out.append(f"{mz} {inten}")
        else:
            (
                exp_int,
                theo_mass,
                ms2_id,
                ch,
                theo_mz,
                inte_flag,
                idx,
                theo_intensity,
                inte_per,
                label
            ) = annot
            
            peaks_out.append(
                f"{mz:.5f} {inten:.2f} {ms2_id+1:d} {theo_mass:.5f} "
                f"{ch:d} {theo_mz:.5f} {inte_flag} {idx + 1:d} "
                f"{theo_intensity:.2f} {inte_per:.3f} {label}"
            )


    return {"meta": meta, "peaks": peaks_out, "meta_lines": msalign_meta_lines}



