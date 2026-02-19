"""
this version includes ms1 peak list
"""
from pyteomics import mzml, mgf
import os
import sys
# import re
from pathlib import Path


def mzML_ms_extract_with_ms1(mzml_filename):
    result = []  
    with mzml.MzML(mzml_filename) as reader:
        last_ms1 = None  # last MS1 spectrum
        for spectrum in reader:
            mz_array = spectrum['m/z array']
            intensity_array = spectrum['intensity array']
            scan_id = int(spectrum['id'][spectrum['id'].find('scan='):][5:])     
            scan_lower = float(str(spectrum['scanList']['scan'][0]['scanWindowList']['scanWindow'][0]['scan window lower limit']))
            scan_upper = float(str(spectrum['scanList']['scan'][0]['scanWindowList']['scanWindow'][0]['scan window upper limit']))
            ret_time = float(spectrum['scanList']['scan'][0]['scan start time']) * 60
            if spectrum['ms level'] == 1:
                last_ms1 = {
                    'ms1_scan_id': scan_id,
                    'ms1_scan_begin': scan_lower,
                    'ms1_scan_end': scan_upper,
                    'ms1_retention_time': ret_time,
                    'ms1_mz_array': mz_array,
                    'ms1_intensity_array': intensity_array
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
                    # ms2_scan_id = int(spectrum['spectrum title'].split(',')[1][spectrum['spectrum title'].split(',')[1].find('scan='):-1][5:])
                    # collison_energy = float(spectrum['precursorList']['precursor'][0]['activation']['collision energy'])
                    pepmass_mz = float(spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z'])
                    charge = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0].get('charge state', None)
                    peak_intensity = (float(spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['peak intensity']) 
                                if 'peak intensity' in spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0] else None)                 
                    collision_energy = spectrum['precursorList']['precursor'][0]['activation'].get('collision energy')
                    total_ion_current = spectrum.get('total ion current')  
                    result.append({
                        'title': title,
                        'ms2_scan_id': scan_num,
                        'ms2_scan_begin': scan_lower,
                        'ms2_scan_end': scan_upper,
                        'ms2_retention_time': ret_time,
                        'pepmass_mz': pepmass_mz,
                        'charge': charge,
                        'intensity': peak_intensity,
                        'collision_energy': collision_energy,
                        'total_ion_current': total_ion_current,
                        'ms2_mz_array': mz_array,
                        'ms2_intensity_array': intensity_array,
                        'ms1_scan_id': last_ms1['ms1_scan_id'],
                        'ms1_scan_begin': last_ms1['ms1_scan_begin'],
                        'ms1_scan_end': last_ms1['ms1_scan_end'],
                        'ms1_retention_time': last_ms1['ms1_retention_time'],
                        'ms1_mz_array': last_ms1['ms1_mz_array'],
                        'ms1_intensity_array': last_ms1['ms1_intensity_array']
                    })
    return result

    
def mgf_write(mzml_spectrum, mgf_filename):
    # write it to a mgf file
    with open(mgf_filename, 'w') as out:
        for spec in mzml_spectrum:
            mgf.write([spec], out)
            


def extend_mgf_write_with_ms1(spectrum2, mzml_filename, mgf_filename):
    count_written = 0
    with open(mgf_filename, 'w') as out:
        for spec in spectrum2:
            out.write("BEGIN IONS\n")
            out.write(f"MZML_FILE_NAME={mzml_filename}\n")
            out.write(f"SCAN={spec['ms2_scan_id']}\n")
            if "title" in spec:
                out.write(f"TITLE={spec['title']}\n")
            #out.write(f"MS_ONE_SCAN={spec['ms1_scan_id']}\n")
            #out.write(f"MS_ONE_RETENTION_TIME={spec['ms1_retention_time']:.5f}\n")              
            #out.write(f"MS_ONE_SCAN_WINDOW_LOWER_LIMIT={spec['ms1_scan_begin']}\n")
            #out.write(f"MS_ONE_SCAN_WINDOW_UPPER_LIMIT={spec['ms1_scan_end']}\n")         
            out.write(f"RTINSECONDS={spec['ms2_retention_time']:.5f}\n")   
            if "pepmass_mz" in spec:
                pepmass_mz = spec["pepmass_mz"]
                out.write(f"PEPMASS_MZ={pepmass_mz}\n")
            if "charge" in spec and spec["charge"] is not None:
                out.write(f"CHARGE={spec['charge']}+\n")
            '''
            if spec.get('collision_energy') is not None:
                out.write(f"COLLISION_ENERGY={spec['collision_energy']}\n")
            if spec.get('total_ion_current') is not None:
                out.write(f"TOTAL_ION_CURRENT={spec['total_ion_current']}\n")   
            
            out.write(f"MS_TWO_SCAN={spec['ms2_scan_id']}\n")
            out.write(f"MS_TWO_SCAN_WINDOW_LOWER_LIMIT={spec['ms2_scan_begin']}\n")
            out.write(f"MS_TWO_SCAN_WINDOW_UPPER_LIMIT={spec['ms2_scan_end']}\n")
            '''
            
            # write peaks (m/z and intensity) for ms1
            '''
            out.write("MS_ONE_PEAKS_BEGIN\n")
            for mz, inten in zip(spec["ms1_mz_array"], spec["ms1_intensity_array"]):
                out.write(f"{mz:.5f} {inten:.2f}\n")
            out.write("MS_ONE_PEAKS_END\n")
            '''
            
            #out.write("MS_TWO_PEAKS_BEGIN\n")
            # write peaks (m/z and intensity) for ms2
            for mz, inten in zip(spec["ms2_mz_array"], spec["ms2_intensity_array"]):
                out.write(f"{mz:.5f} {inten:.2f}\n")
            #out.write("MS_TWO_PEAKS_END\n")
    
            out.write("END IONS\n\n")
            count_written += 1 


if __name__ == "__main__":
    if len(sys.argv) <3:
        print("Usage: python script.py <input_mzml_filename> <output_mgf_filename> <mgf_format>")
        sys.exit()
    else:
        mzml_filename = sys.argv[1] # input the mzML files folder 
        mgf_filename = sys.argv[2]
        if len(sys.argv)==4:
            convert_method = sys.argv[3]
        else:
            convert_method = "extended"
        print(f"converting: {mzml_filename}")
        mzml_spectrum = mzML_ms_extract_with_ms1(mzml_filename)  
        if convert_method == "standard":
            mgf_write(mzml_spectrum, mgf_filename)
        else:
            extend_mgf_write_with_ms1(mzml_spectrum, mzml_filename, mgf_filename)
