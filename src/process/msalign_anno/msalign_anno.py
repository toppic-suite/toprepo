import numpy as np
import time
import argparse
from process.msalign import msalign_reader
from process.msalign import msalign_writer

H_MASS = 1.00782503223
ISOTOPIC_MASS = 1.00235
H2O_MASS = 18.010564683704
PROTON_MASS = 1.007276 
ACETYL_MASS = 42.0106
NH3_MASS = 17.026549
FIXED_PTM_MASS = {
    "Carbamidomethylation": 57.021464
}

AMINO_ACID_MASSES = {
    "A":71.03711,
    "R":156.10111,
    "N":114.04293,
    "D":115.02694,
    "C":103.00918,
    "E":129.04259,
    "Q":128.05858,
    "G":57.02146,
    "H":137.05891,
    "I":113.08406,
    "L":113.08406,
    "K":128.09496,
    "M":131.04049,
    "F":147.06841,
    "P":97.05276,
    "S":87.03203,
    "T":101.04768,
    "W":186.07931,
    "Y":163.06333,
    "V":99.06841, 
    "U":150.95363
}

# We need to double check the shifts
ION_TYPES = {
    'b': 0,
    'a': -26.9871 - H_MASS, 
    'c': 17.0265,
    'x': 44.9977 - H_MASS,
    'y': 18.01056468362, 
    'z': 1.9918 - H_MASS,
    'z_dot': 1.9919
}


# if we need to increase i, return true, otherwise, return false
def increase_I(i, j, deviation, exp_masses, theo_masses): 
    # we assume that each exp peak is matched to at most one theoretical
    # peak, so we do not check i and j+1
    if deviation <= 0:
        return True
    # severl exp peak can be matched to the same theoretical peak
    if i >= len (exp_masses) - 1:
        return False

    next_pos = exp_masses[i + 1]
    if j >= len(theo_masses) - 1:
        return True
    else:
        # check which theoretical mass is closer 
        return abs(next_pos - theo_masses[j]) < abs(next_pos - theo_masses[j + 1])
  
def comp_exp_mass_errors(exp_masses, theo_masses):
    min_distances = [float('inf')] * len(exp_masses)
    theo_indexes = [int(-1)] * len(exp_masses)
    i = 0
    j = 0
    while (i < len(exp_masses) and j < len(theo_masses)):
        d = exp_masses[i] - theo_masses[j];
        #print(f"Comparing exp_masses[{i}]={exp_masses[i]:.4f} to theo_masses[{j}]={theo_masses[j]:.4f}, d={d:.4f}"  )
        if (abs(d) <= abs(min_distances[i])):
            min_distances[i] = d
            theo_indexes[i] = j
        if increase_I(i, j, d, exp_masses, theo_masses): 
            i = i + 1
        else: 
            j = j + 1
    return theo_indexes

# ---------- MATCH OBSERVED TO THEORETICAL ----------
def match_observed_to_theo(exp_mass_table, theo_mass_table, seq, ppm_tol=20.0):
    exp_mass = np.array([e['mass'] for e in exp_mass_table])
    theo_mass = np.array([t['mass'] for t in theo_mass_table])
    theo_indexes = comp_exp_mass_errors(exp_mass, theo_mass)
    anno_peak_lines = []
    bond_num = len(seq) - 1
    coverage = [False] * bond_num
    for i in range(len(exp_mass)):
        if theo_indexes[i] >= 0:
            theo_m = theo_mass[theo_indexes[i]]
            delta_da = exp_mass[i] - theo_m
            delta_ppm = (delta_da / theo_m) * 1e6
            if abs(delta_ppm) <= ppm_tol or abs(delta_da) <= 0.01:
                pos = theo_mass_table[theo_indexes[i]]['pos']
                anno_line = exp_mass_table[i]['line'].strip() + "\t" \
                             f"{theo_mass_table[theo_indexes[i]]['ion']}" + "\t" \
                             f"{theo_mass_table[theo_indexes[i]]['aa_num']}" + "\t" \
                             f"{theo_mass_table[theo_indexes[i]]['pos']}" + "\t" \
                             f"{theo_mass_table[theo_indexes[i]]['shift']}" + "\t" \
                             f"{delta_da:.4f}" + "\t" \
                             f"{delta_ppm:.2f}"
                anno_peak_lines.append(anno_line)
                coverage[pos - 1] = True
                continue    
        anno_peak_lines.append(exp_mass_table[i]['line'])
    # compute sequence coverage
    covered_bonds = sum(1 for c in coverage if c)
    return anno_peak_lines, covered_bonds

def build_mass_table(clean_seq, selected_ions, n_term_acetyl=False, fixed_mod_list=[], unexpected_mod_list=[]):
    n = len(clean_seq)
    # Prepare list of fixed PTMs per position
    fixed_mod_mass_list = [0.0] * n
    if fixed_mod_list is not None:
        for position, mod_name in fixed_mod_list:
            idx = position - 1  # convert to 0-based index
            if 0 <= idx < n:
                if (mod_name in FIXED_PTM_MASS):
                    fixed_mod_mass_list[idx] += FIXED_PTM_MASS[mod_name]
                else:
                    print(f"Warning: modification '{mod_name}' not recognized. Skipping.")
            else:
                print(f"Warning: position {position} out of range for sequence length {n}")
    # Prepare list of unexpected masses per position
    unexpected_mass_list = [0.0] * n
    if unexpected_mod_list is not None:
        for position, mass in unexpected_mod_list:
            idx = position - 1  # convert to 0-based index
            if 0 <= idx < n:
                unexpected_mass_list[idx] += mass
            else:
                print(f"Warning: position {position} out of range for sequence length {n}")
    # Calculate residue masses with modifications
    residue_mass_list = [0.0] * n
    for i in range(n):
        if clean_seq[i] in AMINO_ACID_MASSES:
            residue_mass_list[i] = AMINO_ACID_MASSES[clean_seq[i]] + fixed_mod_mass_list[i] + unexpected_mass_list[i]
        else:
            print(f"Warning: amino acid '{clean_seq[i]}' not recognized. Using mass 0.")
            residue_mass_list[i] = fixed_mod_mass_list[i] + unexpected_mass_list[i]
    # N-terminal acetylation
    if n_term_acetyl:
        residue_mass_list[0] += ACETYL_MASS 

    # Compute prefix and suffix masses 
    prefix_mass_list = [0.0] * n
    prefix_sum = 0.0
    for i in range(n):
        prefix_sum += residue_mass_list[i]
        prefix_mass_list[i] = prefix_sum
    suffix_mass_list = [0.0] * n
    suffix_sum = 0.0
    for i in range(n-1, -1, -1):
        suffix_sum += residue_mass_list[i]
        suffix_mass_list[i] = suffix_sum

    mass_table = []
    for ion in selected_ions:
        shift = selected_ions[ion]
        ion_first_char = ion[0]  # Get the base ion type (e.g., 'b' from 'b-H2O')
        # Build mass series: N-terminal 
        if ion_first_char in ['a', 'b', 'c']:
            for i in range(n - 1):
                aa_num = i + 1
                neutral_mass = prefix_mass_list[i] + shift
                #print(f"Prefix ion: {ion}{aa_num}, mass: {prefix_mass:.4f}")
                mass_table.append({'ion': f"{ion}{aa_num}", 'pos': aa_num, 'aa_num': aa_num, 'shift': 0, 'mass': neutral_mass})
                # ±1.00235 Da variants
                mass_table.append({'ion': f"{ion}{aa_num}+1", 'pos': aa_num, 'aa_num': aa_num, 'shift': 1, 'mass': neutral_mass + ISOTOPIC_MASS})
                mass_table.append({'ion': f"{ion}{aa_num}-1", 'pos': aa_num, 'aa_num': aa_num, 'shift': -1, 'mass': neutral_mass - ISOTOPIC_MASS})

        else:  # x, y, z, zo are from C-terminal
            for i in range(n - 1, 0, -1):
                aa_num = n - i
                neutral_mass = suffix_mass_list[i] + shift
                mass_table.append({'ion': f"{ion}{aa_num}", 'pos': n - aa_num, 'aa_num': aa_num, 'shift': 0, 'mass': neutral_mass})
                # ±1.00235 Da variants
                mass_table.append({'ion': f"{ion}{aa_num}+1", 'pos': n - aa_num, 'aa_num': aa_num, 'shift': 1, 'mass': neutral_mass + ISOTOPIC_MASS})
                mass_table.append({'ion': f"{ion}{aa_num}-1", 'pos': n - aa_num, 'aa_num': aa_num, 'shift': -1, 'mass': neutral_mass - ISOTOPIC_MASS})
    
    # Ensure masses are in increasing order
    mass_table.sort(key=lambda x: x['mass'])
    return mass_table

def parse_proteoform(spectrum_meta): 
    proteoform_str = spectrum_meta.get("PROTEOFORM")
    n_term_acetyl = proteoform_str.startswith('[Acetyl]-')

    #------get fixed modification's position and modification-----------
    fixed_ptms = spectrum_meta.get("FIXED_PTMS")
    fixed_mod_list = []
    if isinstance(fixed_ptms, str) and fixed_ptms != "":
        fixed_ptms = fixed_ptms.strip()
        print(f"FIXED_PTMS: {fixed_ptms}")
        for ptm in fixed_ptms.split(';'):
            ptm_name = ptm.split(':')[0]
            #print(f"Parsed unexpected mass: {mass}")
            pos_str = ptm.split(':')[1].replace('[','').replace(']','')
            #print(f"Parsed unexpected mass position: {pos_str}")
            first_part = pos_str.split('-')[0]
            first_part_int = int(first_part)
            #print(f"First part of position: {first_part}, integer value: {first_part_int}")
            fixed_mod_list.append((first_part_int, ptm_name)) 
    
    #------get unexpected mass shift's mass value and position-------------- 
    mod_list = spectrum_meta.get("UNEXPECTED_MODIFICATIONS")
    # unexpected mass 
    unexpected_mod_list = []
    #print(f"UNEXPECTED_MODIFICATIONS: {mod}")
    # if mod:
    if isinstance(mod_list, str) and mod_list != "":
        mod_list = mod_list.strip()
        for mod in mod_list.split(';'):
            mass = float(mod.split(':')[0])
            #print(f"Parsed unexpected mass: {mass}")
            pos_str = mod.split(':')[1].replace('[','').replace(']','')
            #print(f"Parsed unexpected mass position: {pos_str}")
            first_part = pos_str.split('-')[0]
            first_part_int = int(first_part)
            #print(f"First part of position: {first_part}, integer value: {first_part_int}")
            unexpected_mod_list.append((first_part_int, mass)) 
    else:
        unexpected_mod_list = None
    return n_term_acetyl, fixed_mod_list, unexpected_mod_list
        

def annot_one_spectrum(spectrum, activation_ions=None, ppm_tol=20.0):
    seq = spectrum["meta"].get("DATABASE_SEQUENCE", None)
    if seq in (None, ""):
        spectrum["meta_lines"].append("SEQUENCE_COVERAGE=")
        return spectrum
    activation = spectrum["meta"].get("ACTIVATION", "").lower()
    peak_lines = spectrum.get("peak_lines", [])
    exp_mass_table = []
    for line in peak_lines:
        stripped = line.strip()
        parts = stripped.split()
        if len(parts) >= 4:
            exp_mass_table.append({
                'mass': float(parts[0]),
                'line': line
            })
    exp_mass_table.sort(key=lambda x: x['mass'])
    selected_ions = activation_ions.get(activation, None)
    if selected_ions is None:
        print(f"Warning: No ion types selected for activation method '{activation}'. No annotation will be performed.")
        print(f"Available activation methods: {list(activation_ions.keys())}")
        return spectrum
    n_term_acetyl, fixed_mod_list, unexpected_mod_list = parse_proteoform(spectrum["meta"])
    theo_mass_table = build_mass_table(seq, selected_ions = selected_ions, 
                                       n_term_acetyl=n_term_acetyl, fixed_mod_list=fixed_mod_list, 
                                       unexpected_mod_list=unexpected_mod_list) 
    annot_peak_lines, covered_bonds = match_observed_to_theo(exp_mass_table, theo_mass_table, seq=seq, ppm_tol=ppm_tol)
    spectrum["meta_lines"].append(f"SEQUENCE_COVERAGE={covered_bonds}")
    spectrum["peak_lines"] = annot_peak_lines
    return spectrum

# ---------- WRITE ANNOTATED MSALIGN WITH MULTIPROCESSING ----------
def annot_msalign(input_msalign, output_file, activation_ions, ppm_tol=20.0):
    ms_reader = msalign_reader.MsalignReader(input_msalign)
    ms_writer = msalign_writer.MsalignWriter(output_file)
    count = 0
    start_time = time.time()
    for spectrum in ms_reader.readmsalign_iter():
        spectrum = annot_one_spectrum(spectrum, activation_ions=activation_ions, ppm_tol=ppm_tol)
        ms_writer.write(spectrum)
        count += 1
        if count % 1000 == 0:
            elapsed_time = time.time() - start_time
            print(f"\rAnnotated {count} spectra. Elapsed time: {elapsed_time:.2f} seconds.", end="", flush=True)

def get_ion_list(ion_mode, activation, include_ion_loss):
    # Selecte ion types based on activation method and ion mode
    if ion_mode == "basic":
        if activation in ['hcd', 'cid', 'ethcd']:
            ion_types = ['b', 'y']
        elif activation in ['etd', 'ecd']:
            ion_types = ['c', 'z_dot']
        else:
            print(f"Unknown activation method '{activation}', defaulting to 'b' and 'y' ions.")
    elif ion_mode == "all":
        ion_types = list(ION_TYPES.keys())
    else:
        raise ValueError(f"Invalid ion_mode: {ion_mode}. Choose 'basic' or 'all'.")
    
    if (include_ion_loss):
        ion_losses = {"-H2O": - H2O_MASS,
                      "-NH3": - NH3_MASS}
    else:
        ion_losses = {}
    
    selected_ions = {}
    for ion in ion_types:
        if ion not in ION_TYPES:
            print(f"Warning: ion type '{ion}' not recognized. Skipping.")
            continue
        ion_shift = ION_TYPES[ion]
        selected_ions[ion] = ion_shift
        for loss in ion_losses:
            ion = f"{ion}{loss}"
            selected_ions[ion] = ion_shift + ion_losses[loss]
    return selected_ions 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Annotate MS2 spectra using provided tsv and msalign files.")
    parser.add_argument(
        "--msalign", required=True, type=str, help="Input msalign filename")
    parser.add_argument(
        "--out", required=True, type=str, default="annot.msalign",
        help="Output annotated msalign filename (default: ms2_spectra_annot.msalign)")
    parser.add_argument(
        "--ion_type", required=False, type=str, choices = ['basic', 'all'], help="Ion type (basic/all)", default='basic')
    parser.add_argument(
        "--neutral_loss", required=False, action='store_true', help="Include ion neutral losses (e.g., -H2O, -NH3)")

    args = parser.parse_args()
    output_filename = args.out or "ms2_spectra_annot.msalign"

    ion_mode = args.ion_type
    include_ion_losses = args.neutral_loss

    # Prepare ion types for each activation method
    activation_ions = {}
    for activation in ['hcd', 'cid', 'etd', 'ecd', 'ethcd']:
        selected_ions = get_ion_list(ion_mode, activation, include_ion_losses) 
        activation_ions[activation] = selected_ions
        #print(f"Annotating spectra with activation method: {activation}")

    annot_msalign(args.msalign, output_filename, activation_ions)    