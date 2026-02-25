import json
import bisect
from collections import defaultdict
from time import time
# from xxlimited import new

# read theoretical envelope file "theo_patt.txt" 
def read_theo_patterns(file_path):
    """
    Reads the theo_patt.txt file and returns a list of envelopes.
    Each envelope is a dict with 'mass' and 'peaks': [(mass, intensity), ...]
    """
    envelopes = []
    with open(file_path, 'r') as f:
        lines = f.readlines()

    mono_mass_list = []
    i = 0
    while i < len(lines):
        if lines[i].startswith("formula:"):
            header = lines[i].strip()
            parts = header.split()
            mono_mass = float(parts[-1])
            i += 1
            peaks = []
            while i < len(lines) and lines[i].strip():
                m, inten = map(float, lines[i].strip().split())
                peaks.append((m, inten))
                i += 1
            mono_mass_list.append(mono_mass)               
            envelopes.append({'mono_mass': mono_mass, 'peaks': peaks})
        i += 1
    return envelopes, mono_mass_list


_THEO_CACHE = {}

def get_theo_envelopes_cached(theo_file):
    """
    read all theoretical envelopes out in the file "theo_patt.txt" and store it as a dict
    """
    if theo_file not in _THEO_CACHE:
        _THEO_CACHE[theo_file] = read_theo_patterns(theo_file)
    return _THEO_CACHE[theo_file]

# format conversion
def safe_json_load(value):
    """convert JSON string to Python object if it's a string; otherwise return as is."""
    if isinstance(value, str):
        return json.loads(value)
    return value

# find closest theoretical envolope according to the experimental fragment mass 
def find_closest_env(envelopes, mono_mass_list, target_mono_mass):
    """
    Finds the envelope whose monoisotopic mass is closest to target_mono_mass.
    Return the closest envolope index and corresponding 'mass' and 'peaks' in all envelopes
    """
    idx = bisect.bisect_left(mono_mass_list, target_mono_mass) # return nearest index
    closest_env = envelopes[idx]
    return closest_env


def get_annotated_mz_intensity(mono_mass, charge, envelopes, mono_mass_list):
    """
    Given a monoisotopic monoisotopic mass and charge, return a list of (m/z, intensity) values for the envelope.
    """
    env = find_closest_env(envelopes, mono_mass_list, mono_mass)
    mono_mass_err = env['mono_mass'] - mono_mass # get the mass difference between closest theoretical envelope and fragment mass.
    mz_err = mono_mass_err / charge
    PROTON_MASS = 1.007276 
    # Now recalculate m/z values from env
    annotated_peaks = []
    total_intensity = sum(intensity for mass, intensity in env['peaks']) # calculate the total theoretical intensity
    
    for i, (mass, intensity) in enumerate(env['peaks']):
        mz = (mass / charge) + PROTON_MASS  # convert mass to m/z value  
        mz -= mz_err # offset the m/z error 
        intensity_percent = intensity / total_intensity  # calculate the ratio of each intensity to total intensity
        annotated_peaks.append((round(mz, 5), intensity, round(intensity_percent, 3)))  # store annotated results (mz with 5 decimal places and intensity ratio with 3 decimal places) 
    return annotated_peaks


def get_ms2_centroid_label(theo_file, form_df, ppm_tol):
    """
    Annotate centroided MS/MS spectra using theoretical isotopic peaks.
    
    For each centroided MS/MS spectrum in 'form_df', theoretical peaks from 'theo_file' are matched to experimental peaks within a given mass tolerance 'ppm_tol'. 
   
    Parameters
    ----------
    theo_file : str 
        file name: 'theo_patt.txt'
    form_df : pandas.DataFrame
        DataFrame containing centroided MS/MS spectra and deconvoluted masses.
    ppm_tol : float
        Mass error tolerance in ppm, ppm=20 by default.

    Returns
    -------
        The input DataFrame with an additional column containing centroid spectrum annotations in nested lists.
    """
    envelopes, mono_mass_list = get_theo_envelopes_cached(theo_file)
    
    ms2_centroid_label_all = []
    PROTON_MASS = 1.007276 
    #time_1 = 0
    #time_2 = 0
    #time_3 = 0
    for ss in range(len(form_df)):
        # if ss % 100 == 0:
            # print(ss)
        time_1_start = time()
        ms2_mass_value = form_df['mass_all'].iloc[ss]                     # column 'mass_all' stores MS2 deconvoluted masses in a list 
        ms2_inte_value = form_df['intensity_all'].iloc[ss]                # column 'intensity_all' stores MS2 deconvoluted intensities in a list
        ms2_ch_value = form_df['charge_all'].iloc[ss]                     # column 'charge_all' stores MS2 deconvoluted charges in a list
        ms2_mz_centroid_value = form_df['mz_array'].iloc[ss]              # column 'mz_array' stores centriod m/z values in a list
        ms2_inte_centroid_value = form_df['intensity_array'].iloc[ss]     # column 'intensity_array' stores centriod intensity values in a list
        ms2_deconv_label_value = form_df['ms2_deconv_label'].iloc[ss]     # column 'ms2_deconv_label' stores annotations of deconvoluted fragment peaks in a list, ex. 'b-H2O 8 0 0.0001 0.1150'
        
        # format conversion to ensure values are loaded as Python lists
        ms2_mass_list = safe_json_load(ms2_mass_value)
        ms2_inte_list = safe_json_load(ms2_inte_value)
        ms2_ch_list = safe_json_load(ms2_ch_value)
        ms2_mz_centroid_list = safe_json_load(ms2_mz_centroid_value)
        ms2_inte_centroid_list = safe_json_load(ms2_inte_centroid_value)
        ms2_deconv_label_list = safe_json_load(ms2_deconv_label_value)               

        #time_1 += time() - time_1_start

        #time_2_start = time()
        all_theo_mz = []
        all_theo_data = []  
    
        for index, (mass, intensity, charge) in enumerate(zip(ms2_mass_list, ms2_inte_list, ms2_ch_list)):
            annotated_peaks = get_annotated_mz_intensity(mass, charge, envelopes, mono_mass_list) # Given a (mass, charge), get annotated theoretical peaks.  
            ann_idx_list = range(len(annotated_peaks))
            ann_idx = 0
            
            for mz, inte, inte_per in annotated_peaks:
                theo_mass = charge * (mz - PROTON_MASS) # convert m/z to mass 
                if inte == 100:  # if theoretical isotopic intensity is "100", then set a flag to "yes" meaning the maximum intensity, otherwise, "no".  
                    flag = 'yes'
                else:
                    flag = 'no'
                
                theo_intensity = intensity * inte_per  # compute theoretical peak intensity by experimental intensity × theoretical isotopic relative intensity  
                
                all_theo_mz.append(mz)
                all_theo_data.append((mz, round(theo_mass, 5), index, charge, flag, ann_idx_list[ann_idx], theo_intensity, inte_per)) 
                #  (theoretical m/z, theoretical mass, MS2 peak index, charge state, isotopic envelope ID, monoisotopic flag, isotopic peak index, 
                #  theoretical intensity, isotopic relative intensity)
                ann_idx += 1
        #time_2 += time() - time_2_start
        
        #time_3_start = time()
        # Sort by mz for bisect search (faster)
        sorted_indices = sorted(range(len(all_theo_mz)), key=lambda x: all_theo_mz[x]) # sort all theoretical m/z values
        sorted_mz = [all_theo_mz[i] for i in sorted_indices]
        sorted_data = [all_theo_data[i] for i in sorted_indices]  
            
        # For each centroid m/z, find if any match exists
        ms2_centroid_annot = []

        # start to annotate for centroid peaks within error tolerance (ppm=20)
        for exp_mz, exp_inte in zip(ms2_mz_centroid_list, ms2_inte_centroid_list):
            matched = ""
            tol_da = exp_mz * ppm_tol * 1e-6
            
            idx = bisect.bisect_left(sorted_mz, exp_mz - tol_da) # return nearest index
            while idx < len(sorted_mz) and sorted_mz[idx] <= exp_mz + tol_da:
                if abs(sorted_mz[idx] - exp_mz) <= tol_da: 
                    theo_mz, matched_mass, ms2_id, matched_charge, inte_flag, matched_idx, theo_intensity, inte_per = sorted_data[idx]
                    if ms2_deconv_label_list[ms2_id] != "": 
                        # check if the deconvoluted fragment mass has a valid annotation (not '?'), then include it to centroid annotation
                        matched = (exp_inte, matched_mass, ms2_id, matched_charge, theo_mz, inte_flag, matched_idx, theo_intensity, inte_per, ms2_deconv_label_list[ms2_id])
                    else:
                        matched = (exp_inte, matched_mass, ms2_id, matched_charge, theo_mz, inte_flag, matched_idx, theo_intensity, inte_per, "")
                    break
                idx += 1
            ms2_centroid_annot.append(matched)
    
        ms2_centroid_label_all.append(ms2_centroid_annot)
        #time_3 += time() - time_3_start
        #print(f"Processed spectrum {ss+1}/{len(form_df)}: annotation time = {time_1:.2f} and {time_2:.2f} seconds, search time = {time_3:.2f} seconds" )
    return ms2_centroid_label_all


def filter_ms2_centroid_labels(ms2_centroid_label_all):
    """
    Filter matched centroid annotations based on theoretical intensity thresholds.
    Parameters:
    ms2_centroid_label_all: containing annotated lists
    """
    filtered_ms2_centroid_label_all = []

    for annot_list in ms2_centroid_label_all:
        # Step 1: Group all experimental intensities by ms2_id
        id_to_exp_intensities = defaultdict(list)

        for entry in annot_list:
            if isinstance(entry, tuple) and len(entry) > 7:
                exp_intensity = entry[0]
                ms2_id = entry[2]
                id_to_exp_intensities[ms2_id].append(exp_intensity)

        # Step 2: Compute min intensities
        id_to_min_exp_intensity = {}
        for ms2_id, intensities in id_to_exp_intensities.items():
            id_to_min_exp_intensity[ms2_id] = min(intensities)


        # Step 3: Filter out matches with theo_intensity < min_exp_intensity
        filtered_list = []
        for entry in annot_list:
            if isinstance(entry, tuple) and len(entry) > 7:
                ms2_id = entry[2]
                theo_intensity = entry[7]
                min_exp_intensity = id_to_min_exp_intensity.get(ms2_id, 0)
                if theo_intensity >= min_exp_intensity: # if greater, then valid entry
                    filtered_list.append(entry)
                else:
                    filtered_list.append("") # otherwise, set it to be empty
            else:
                # Keep unmatched (empty string) entries
                filtered_list.append(entry)

        filtered_ms2_centroid_label_all.append(filtered_list)
        
    return filtered_ms2_centroid_label_all



