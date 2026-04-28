import numpy as np
import time
import argparse
import re
import pandas as pd
from process.msalign import msalign_reader
from process.msalign import msalign_writer


ISOTOPIC_MASS = 1.00235
H2O_MASS =     18.0105647
ACETYL_MASS =  42.010565
NH3_MASS =     17.0265491
FIXED_PTM_MASS = {
    "Carbamidomethylation": 57.021464
}

AMINO_ACID_MASSES = {
    "A":71.0371137846,
    "R":156.10111102302,
    "N":114.04292744016,
    "D":115.0269430228,
    "C":103.0091849595,
    "E":129.04259308732,
    "Q":128.05857750468,
    "G":57.02146372008,
    "H":137.05891185752,
    "I":113.08406397816,
    "L":113.08406397816,
    "K":128.09496301462,
    "M":131.04048508854,
    "F":147.06841391364,
    "P":97.05276384912,
    "S":87.0320284037,
    "T":101.04767846822,
    "W":186.0793129501,
    "Y":163.06332853274,
    "V":99.06841391364,
    "U":150.9536363846,
    "O":132.08987763372
}


# B: 0
# a: -CO       - 27.9949146
# c: +NH3      + 17.0265491
# y: +H2O      + 18.0105647
# x: +O + CO   + 43.9898292
# z: H2O - NH3 +  0.9840156
# z_dot: Z + H +  1.9918406
ION_TYPES = {
    'b':     0,
    'a':   -27.9949146, 
    'c':    17.0265491,
    'x':    43.9898292,
    'y':    18.0105647, 
    'z':     1.9940156,
    'z_dot': 1.9918406
}

LOSS_MASS = {
    "-H2O": -H2O_MASS,
    "-NH3": -NH3_MASS
}


def generalize_ion(ion):
    ion = ion.strip()      
    match = re.match(r"^(z_dot|[abcxyz])(\d+)([-+].*)?$", ion)
    if match:
        ion_type = match.group(1)
        suffix = match.group(3) if match.group(3) else ""
        return ion_type + suffix
    return ion

def load_annotation_tables(tsv_path):
    anno_df = pd.read_csv(tsv_path, sep='\t')
    act_types = anno_df['activation'].unique()
    label_dict = {}
    for act in act_types:
        key = act.lower()
        df = anno_df[anno_df['activation'] == act]
        label_dict[key] = {
            label: (count, coverage)
            for label, count, coverage in zip(df["label"], df["count"], df["coverage"])
        }

    return label_dict

def ion_mode_selection(ion_mode, losses):
    ion_mode = ion_mode.lower()

    if ion_mode == "basic":
        if losses:
            selected_losses = ["-H2O"]   
        else:
            selected_losses = []

    elif ion_mode == "all":
        selected_losses = ["-H2O", "-NH3"]

    else:
        raise ValueError("please select ion mode = 'basic' or 'all'.")

    return ion_mode, selected_losses


def match_all_annotations(exp_mass_table, theo_mass_table, seq, activation,
                          label_freq_dict, cov_rate, ppm_tol=20.0, da_tol=0.01):

    anno_peak_lines = []
    bond_num = len(seq) - 1
    coverage = [False] * bond_num

    i = 0
    j = 0

    while i < len(exp_mass_table) and j < len(theo_mass_table):

        exp = exp_mass_table[i]
        theo = theo_mass_table[j]

        exp_m = exp['mass']
        theo_m = theo['mass']
        if theo_m <= 0:
            j += 1
            continue

        delta_da = exp_m - theo_m
        delta_ppm = (delta_da / theo_m) * 1e6

        if (
            (ppm_tol is not None and abs(delta_ppm) <= ppm_tol) or
            (da_tol is not None and abs(delta_da) <= da_tol)
        ):

            matches = []

            k = j

            # scan neighboring theoretical ions
            while k < len(theo_mass_table):

                theo2 = theo_mass_table[k]

                delta_da = exp_m - theo2['mass']
                delta_ppm = (delta_da / theo2['mass']) * 1e6

                if (
                    (ppm_tol is not None and abs(delta_ppm) <= ppm_tol) or
                    (da_tol is not None and abs(delta_da) <= da_tol)
                ):

                    ion_label = theo2['ion'] # e.g. y12-H2O-1
                    ion_gen_label = generalize_ion(ion_label)  # y-H2O-1
                    
                    # Skip if no table
                    if activation.lower() not in label_freq_dict:
                        k += 1
                        continue
                    
                    freq_dict = label_freq_dict[activation.lower()]
                    # only keep if label exists
                    if ion_gen_label not in freq_dict:
                        k += 1
                        continue
                    
                    matches.append({
                        'ion': ion_label,
                        'aa_num': theo2['aa_num'],
                        'pos': theo2['pos'],
                        'shift': theo2['shift'],
                        'delta_da': delta_da,
                        'delta_ppm': delta_ppm,
                        'freq': freq_dict[ion_gen_label][0],
                        'rate': float(freq_dict[ion_gen_label][1].rstrip('%'))
                    })

                    k += 1
                else:
                    break

            base_cols = exp['cols']

            if matches:
                # pick highest frequency annotation
                best = max(matches, key=lambda x: x['freq'])
                if best['rate'] >= cov_rate:
                    coverage[best['pos'] - 1] = True
                    flat = [
                        best['ion'],
                        str(best['aa_num']),
                        str(best['pos']),
                        str(best['shift']),
                        f"{best['delta_da']:.4f}",
                        f"{best['delta_ppm']:.2f}"
                    ]
    
                    line_out = "\t".join(base_cols + flat)
                else:
                    line_out = "\t".join(base_cols)
            else:            
                line_out = "\t".join(base_cols)

            anno_peak_lines.append(line_out)
            i += 1

        elif exp_m < theo_m:        
            anno_peak_lines.append("\t".join(exp['cols']))
            i += 1
        else:
            j += 1

    # remaining experimental peaks
    while i < len(exp_mass_table):
        exp = exp_mass_table[i]
        anno_peak_lines.append("\t".join(exp['cols']))
        i += 1
        
    covered_bonds = sum(coverage)

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
        # ion_first_char = ion[0]  # Get the base ion type (e.g., 'b' from 'b-H2O')
        # -------- split base ion and loss --------
        if "-" in ion or "+" in ion:
            if "-" in ion:
                base, loss = ion.split("-", 1)
                loss = "-" + loss
            else:
                base, loss = ion.split("+", 1)
                loss = "+" + loss
        else:
            base = ion
            loss = ""
        
        ion_first_char = base[0]
        
        # Build mass series: N-terminal 
        if ion_first_char in ['a', 'b', 'c']:
            for i in range(n - 1):
                aa_num = i + 1
                neutral_mass = prefix_mass_list[i] + shift
                #print(f"Prefix ion: {ion}{aa_num}, mass: {prefix_mass:.4f}")
                mass_table.append({'ion': f"{base}{aa_num}{loss}", 'pos': aa_num, 'aa_num': aa_num, 'shift': 0, 'mass': neutral_mass})
                # ±1.00235 Da variants
                mass_table.append({'ion': f"{base}{aa_num}{loss}+1", 'pos': aa_num, 'aa_num': aa_num, 'shift': 1, 'mass': neutral_mass + ISOTOPIC_MASS})
                mass_table.append({'ion': f"{base}{aa_num}{loss}-1", 'pos': aa_num, 'aa_num': aa_num, 'shift': -1, 'mass': neutral_mass - ISOTOPIC_MASS})

        else:  # x, y, z, zo are from C-terminal
            for i in range(n - 1, 0, -1):
                aa_num = n - i
                neutral_mass = suffix_mass_list[i] + shift
                mass_table.append({'ion': f"{base}{aa_num}{loss}", 'pos': n - aa_num, 'aa_num': aa_num, 'shift': 0, 'mass': neutral_mass})
                # ±1.00235 Da variants
                mass_table.append({'ion': f"{base}{aa_num}{loss}+1", 'pos': n - aa_num, 'aa_num': aa_num, 'shift': 1, 'mass': neutral_mass + ISOTOPIC_MASS})
                mass_table.append({'ion': f"{base}{aa_num}{loss}-1", 'pos': n - aa_num, 'aa_num': aa_num, 'shift': -1, 'mass': neutral_mass - ISOTOPIC_MASS})
    
    # Ensure masses are in increasing order
    mass_table.sort(key=lambda x: x['mass'])
    return mass_table

def parse_proteoform(spectrum_meta): 
    proteoform_str = spectrum_meta.get("PROTEOFORM")
    n_term_acetyl = proteoform_str.startswith('[Acetyl]-')

    #------get fixed modification's position and modification-----------
    fixed_ptms = spectrum_meta.get("FIXED_MODIFICATIONS")
    fixed_mod_list = []
    if isinstance(fixed_ptms, str) and fixed_ptms != "":
        fixed_ptms = fixed_ptms.strip()
        # print(f"FIXED_MODIFICATIONS: {fixed_ptms}")
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
        

def annot_one_spectrum(spectrum, label_freq_dict, cov_rate, activation_ions=None, ppm_tol=20.0, da_tol=0.01):
    seq = spectrum["meta"].get("DATABASE_SEQUENCE", None)
    if seq in (None, ""):
        spectrum["meta_lines"].append("SEQUENCE_COVERAGE=")
        return spectrum
    activation = spectrum["meta"].get("ACTIVATION", "").lower()
    peak_lines = spectrum.get("peak_lines", [])
    exp_mass_table = []
    for line in peak_lines:
        stripped = line.strip()
        # stripped = line.rstrip("\n") 
        parts = stripped.split()
        if len(parts) >= 4:
            exp_mass_table.append({
                'mass': float(parts[0]),
                'cols': parts[:4]
                # 'line': line
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
    annot_peak_lines, covered_bonds = match_all_annotations(exp_mass_table, theo_mass_table, seq, activation, label_freq_dict, cov_rate, ppm_tol=ppm_tol, da_tol=da_tol)      
    spectrum["meta_lines"].append(f"SEQUENCE_COVERAGE={covered_bonds}")
    spectrum["peak_lines"] = annot_peak_lines
    return spectrum

# ---------- WRITE ANNOTATED MSALIGN WITH MULTIPROCESSING ----------
def annot_msalign(input_msalign, label_freq_dict, cov_rate, output_file, activation_ions, ppm_tol=20.0, da_tol=0.01):
    ms_reader = msalign_reader.MsalignReader(input_msalign)
    ms_writer = msalign_writer.MsalignWriter(output_file)
    
    count = 0
    start_time = time.time()
    for spectrum in ms_reader.readmsalign_iter():
        spectrum = annot_one_spectrum(spectrum, label_freq_dict, cov_rate, activation_ions=activation_ions, ppm_tol=ppm_tol, da_tol=da_tol)
        ms_writer.write(spectrum)
        count += 1
        if count % 1000 == 0:
            elapsed_time = time.time() - start_time
            print(f"\rAnnotated {count} spectra. Elapsed time: {elapsed_time:.2f} seconds.", end="", flush=True)

def get_ion_list(ion_mode, activation, ion_losses):
    # Selecte ion types based on activation method and ion mode
    if ion_mode=='basic':
        if activation in ['hcd', 'cid']:
            ion_types= ['b','y']
        else:
            ion_types = ['c','z_dot']
    else:
        ion_types = ['a','b', 'c', 'x', 'y','z','z_dot']
    
    selected_ions = {}
    for ion in ion_types:
        if ion not in ION_TYPES:
            print(f"Warning: ion type '{ion}' not recognized. Skipping.")
            continue
        ion_shift = ION_TYPES[ion]
        selected_ions[ion] = ion_shift
        for loss in ion_losses:
            if loss not in LOSS_MASS:
                print(f"Warning: loss type '{loss}' not recognized. Skipping.")
                continue
            
            ion_name = f"{ion}{loss}"
            selected_ions[ion_name] = ion_shift + LOSS_MASS[loss]
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
        "--table", required=True, type=str, default="combine_anno_table.tsv",
        help="A TSV filename containing ion type and frequency")
    parser.add_argument(
        "--ion_type", required=False, type=str, choices = ['basic', 'all'], help="Ion type (basic/all)", default='all')
    parser.add_argument(
        "--neutral_loss", required=False, action='store_true', help="Include ion neutral losses (e.g., -H2O, -NH3)")
    parser.add_argument(
        "--error_tol", required=False, type=float, help="Error tolerance in ppm", default=20.0)
    parser.add_argument("--da_tol", required=False, type=float, help="Error tolerance in Da", default=0.01)
    parser.add_argument("--rate", required=False, type=float, help="Coverage rate (between 0-0.2)", default=0.0)

    args = parser.parse_args()
    output_filename = args.out or "ms2_spectra_annot.msalign"

    ion_mode, selected_losses = ion_mode_selection(args.ion_type, args.neutral_loss)
    ppm_tol = args.error_tol
    da_tol = args.da_tol
    cov_rate = float(args.rate)*100
    
    label_freq_dict = load_annotation_tables(args.table)
    # Prepare ion types for each activation method
    activation_ions = {}
    for activation in ['hcd', 'cid', 'etd', 'ecd', 'ethcd']:
        selected_ions = get_ion_list(ion_mode, activation, selected_losses)
        activation_ions[activation] = selected_ions
        # if activation == 'hcd':
            # print(f"Annotating spectra with activation method: {activation_ions[activation]}")
    annot_msalign(args.msalign, label_freq_dict, cov_rate, output_filename, activation_ions, ppm_tol, da_tol)    
