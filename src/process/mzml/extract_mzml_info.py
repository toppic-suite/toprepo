"""
this version extracts all avaiable information from mzml files

"""
from pyteomics import mzml
import os
import sys
import re
from pathlib import Path
import pandas as pd


def get_instrument_name_safe(reader):
    try:
        reader.reset()
        elem = next(reader.iterfind('referenceableParamGroup[@id="CommonInstrumentParams"]'))
    except StopIteration:
        return None  # no such section found

    for key in elem.keys():
        key_lower = key.lower()
        # skip meta or serial number fields
        if key_lower in ('id', 'instrument serial number', 'accession', 'cvparam'):
            continue
        # skip empty or trivial keys
        if not key.strip():
            continue
        return key

    return None



def mzML_ms_extract_with_ms1(dataset_id, mzml_filename):
    result = []  
    with mzml.MzML(mzml_filename) as reader:
        instrument_name = get_instrument_name_safe(reader)
        print(f"instrument: {instrument_name}")
        last_ms1 = None  # last MS1 spectrum
        for spectrum in reader:
            # mz_array = spectrum['m/z array']
            # intensity_array = spectrum['intensity array']
            scan_id = int(spectrum['id'][spectrum['id'].find('scan='):][5:])     
            scan_lower = float(str(spectrum['scanList']['scan'][0]['scanWindowList']['scanWindow'][0]['scan window lower limit']))
            scan_upper = float(str(spectrum['scanList']['scan'][0]['scanWindowList']['scanWindow'][0]['scan window upper limit']))
            ret_time = float(spectrum['scanList']['scan'][0]['scan start time']) * 60
            total_ion_current = spectrum.get('total ion current')  
            mass_solving_power = (float(spectrum['scanList']['scan'][0]['mass resolving power']) 
                             if 'mass resolving power' in spectrum['scanList']['scan'][0] else None)                 
            ion_injection_time = float(spectrum['scanList']['scan'][0]['ion injection time'])
            lower_obsevered_mz = spectrum['lowest observed m/z']
            highest_obsevered_mz = spectrum['highest observed m/z'] 
            mzml_fullname = os.path.basename(mzml_filename)
            if dataset_id in mzml_fullname:
                mzml_filename_extract = mzml_fullname.replace(f"{dataset_id}_", "", 1)
            else:
                mzml_filename_extract = mzml_fullname
            
            if spectrum['ms level'] == 1:
                last_ms1 = {
                    'ms1_scan_id': scan_id,
                    'ms1_scan_begin': scan_lower,
                    'ms1_scan_end': scan_upper,
                    'ms1_retention_time': ret_time,
                    'ms1_injection_time': ion_injection_time,
                    'ms1_resolution': mass_solving_power,
                    'ms1_total_ion_current':  total_ion_current,
                    'ms1_lower_obsevered_mz': lower_obsevered_mz,
                    'ms1_highest_obsevered_mz': highest_obsevered_mz,
                    # 'ms1_mz_array': mz_array,
                    # 'ms1_intensity_array': intensity_array
                }
            elif spectrum['ms level'] == 2:
                if last_ms1 is not None:
                    if 'spectrum title' in spectrum:
                        title = spectrum['spectrum title'].split(',')[0]
                        scan_str = spectrum['spectrum title'].split(',')[1]
                        scan_num = int(scan_str[scan_str.find('scan='):-1][5:])
                    else:
                        id_str = spectrum.get('id', '')
                        if 'scan=' in id_str:
                            try:
                                scan_num = int(id_str.split('scan=')[-1])
                            except:
                                scan_num = None
                        else:
                            scan_num = None
                        title = Path(mzml_filename).stem  
                    selected_ion_mz = float(spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z'])
                    selected_charge = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0].get('charge state', None)
                    peak_intensity = (float(spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['peak intensity']) 
                                if 'peak intensity' in spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0] else None)                 
                    collision_energy = float(spectrum['precursorList']['precursor'][0]['activation'].get('collision energy'))
                    iso_win_target_mz = float(spectrum['precursorList']['precursor'][0]['isolationWindow']['isolation window target m/z'])
                    iso_win_lower_offset = float(spectrum['precursorList']['precursor'][0]['isolationWindow']['isolation window lower offset'])
                    iso_win_upper_offset = float(spectrum['precursorList']['precursor'][0]['isolationWindow']['isolation window upper offset'])
                    # activation = spectrum['precursorList']['precursor'][0]['activation'].get('collision energy')

                    result.append({
                        'dataset_id': dataset_id,
                        'instrument': instrument_name,
                        'file_name': mzml_filename_extract,
                        'title': title,
                        'ms2_scan_id': scan_num,
                        'ms2_scan_begin': scan_lower,
                        'ms2_scan_end': scan_upper,
                        'ms2_retention_time': ret_time,
                        'pepmass_mz': selected_ion_mz,
                        'selected_ion_charge': selected_charge,
                        'peak_intensity': peak_intensity,
                        'collision_energy': collision_energy,
                        'ms2_total_ion_current': total_ion_current,
                        'ms2_lower_obsevered_mz': lower_obsevered_mz,
                        'ms2_highest_obsevered_mz': highest_obsevered_mz,
                        'ms2_injection_time': ion_injection_time,
                        'ms2_resolution': mass_solving_power,
                        'isolation_window_mz': iso_win_target_mz,
                        'isolation_window_lower_offset': iso_win_lower_offset,
                        'isolation_window_upper_offset': iso_win_upper_offset,
                        # 'ms2_mz_array': mz_array,
                        # 'ms2_intensity_array': intensity_array,
                        'ms1_scan_id': last_ms1['ms1_scan_id'],
                        'ms1_scan_begin': last_ms1['ms1_scan_begin'],
                        'ms1_scan_end': last_ms1['ms1_scan_end'],
                        'ms1_retention_time': last_ms1['ms1_retention_time'],
                        'ms1_total_ion_current': last_ms1['ms1_total_ion_current'],
                        'ms1_injection_time': last_ms1['ms1_injection_time'],
                        'ms1_resolution': last_ms1['ms1_resolution'],
                        'ms1_lower_obsevered_mz': last_ms1['ms1_lower_obsevered_mz'],
                        'ms1_highest_obsevered_mz': last_ms1['ms1_highest_obsevered_mz']
                        # 'ms1_mz_array': last_ms1['ms1_mz_array'],
                        # 'ms1_intensity_array': last_ms1['ms1_intensity_array']
                    })
    return result

 

def process_mzml_folder(dataset_id, mzml_filename: str, output_filename: str):

    print(f"Extracting metadata from: {mzml_filename} ({dataset_id})")

    result = mzML_ms_extract_with_ms1(dataset_id, mzml_filename)
        
    if not result:
        print("No spectra found across all mzML files.")
        return

    df = pd.DataFrame.from_records(result)
    df.to_csv(output_filename, sep='\t', index=False)
    print(f"\nMetadata extracted for {len(df)} MS2 spectra from {mzml_filename}.")
    print(f"Saved to: {output_filename}")



if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <dataset_id> <input_mzml_filename> <output_tsv_filename>")
        sys.exit(1)

    dataset_id = sys.argv[1]
    input_mzml_filename = sys.argv[2]
    output_tsv_filename = sys.argv[3]
    process_mzml_folder(dataset_id, input_mzml_filename, output_tsv_filename)

        
    
