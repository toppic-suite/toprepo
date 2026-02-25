# TopRepo

TopRepo is a top-down spectral repository containing more than 18 million MS/MS spectra from 12 species. For each MS raw file in TopRepo, msconvert is used to convert the raw file to a centroided mzML file, and then TopFD is employed to deconvolute spectra in the centroided mzML file to one or several msalign files and proteoform feature files. The msalign files are searched against its corresponding proteome sequence database for spectral identification using TopPIC. The identification results are stored in TSV files. 

Using mzML files, msalign files, feature files, and spectral identification (TSV) files generated from the data analysis pipeline, Python scripts in this repository are used to generate TSV files with comprehensive spectral information, annotated msalign files, and annotated mgf files. 

## 1. Generate TSV files with comprehensive spectral information

We use an example mzML file spectra.mzML with dataset id PXD029703, and its corresponding msalign file spectra_ms2.msalign, feature file spectra_ms2.feature, and spectral identification file spectra_ms2_toppic_prsm_single.tsv to explain the method.

**1.1 Extract spectral information from mzML file**

This step extracts spectral information from the mzML file and saves it into a TSV file.

* Usage:
python3 extract_mzml_info.py <dataset_id> <input_mzml_filename> <output_tsv_filename>

* Input:
   * dataset_id: MS dataset ID
   * input_mzml_filename: Input mzML filename
   * output_tsv_filename: Output TSV filename for storing spectral information extracted from mzML

Run the command: 
```
python3 toprepo/src/process/mzml/extract_mzml_info.py PXD029703 spectra.mzML spectra_mzml_info.tsv
```

**1.2 Extract spectral information from msalign file**

This step extracts MS2 spectral information from the msalign file and saves it into a TSV file.

* Usage python3 extract_msalign_info.py <dataset_id> <input_msalign_filename> <output_tsv_filename>

* Input:
   * dataset_id: MS dataset ID
   * input_msalign_filename: Input msalign filename
   * output_tsv_filename: Output TSV filename for storing extracted spectral information

Run the command:
```
python3 toprepo/src/process/msalign/extract_msalign_info.py PXD029703 spectra_ms2.msalign spectra_msalign_info.tsv
```

**1.3 Extract feature information from feature file**

This step extracts MS2 feature information (e.g., feature intensity, feature score, feature apex time) from the feature file and saves it into a TSV file.

* Usage python3 extract_feature_info.py <dataset_id> <input_feature_filename> <output_tsv_filename>

* Input:
   * dataset_id: MS dataset ID
   * input_feature_filename: Input feature filename
   * output_tsv_filename: Output TSV filename for storing extracted feature information

Run the command:
```
python3 toprepo/src/process/feature/extract_feature_info.py PXD029703 spectra_ms2.feature spectra_feature_info.tsv
```

**1.4 Preprocess TSV file containing PrSM identifications reported by TopPIC**

This step preprocesses the TSV file containing PrSM identifications and add a "dataset id" column to the file.

* Usage python3 prsm_preprocess.py <prsm_tsv_filename> <dataset_id> --output <output_tsv_filename>

* Input:
   * prsm_tsv_filename: Path to the TSV file containing spectral identifications reported by TopPIC
   * dataset_id: MS dataset id 
   * output_tsv_filename: Path to the TSV output filename 

Run the command:
```
python3 toprepo/src/process/prsm/prsm_preprocess.py spectra_ms2_toppic_prsm_single.tsv PXD029703 --output spectra_toppic_info.tsv
```

**1.5 Merge msalign spectral and feature information**

This step merges msalign spectral and feature information into a single TSV file.

* Usage python3 merge_msalign_feature_info.py <msalign_info_filename> <feature_info_filename> <output_tsv_filename>

* Input:
   * msalign_info_filename: Path to the TSV file containing MS2 spectral information extracted from msalign files.
   * feature_info_filename: Path to the TSV file containing MS2 feature information extracted from feature files.
   * output_tsv_filename: Path to the output TSV file where the merged spectral information will be written.

Run the command:
```
python3 toprepo/src/process/tsv/merge_msalign_feature_info.py spectra_msalign_info.tsv spectra_feature_info.tsv spectra_msalign_feature_info.tsv
```

**1.6 Merge spectral information extracted from mzML, msalign, and feature files**

This step merges the spectral information extracted from mzML, msalign, and feature files into a single TSV file.

* Usage python3 merge_mzml_msalign_info.py <mzml_info_filename> <msalign_feature_info_filename> <output_tsv_filename>

* Input:
   * mzml_info_filename: Path to the TSV file containing spectral information extracted from mzML files
   * msalign_feature_info_filename: Path to the TSV file containing MS2 spectral and feature information extracted from msalign and feature files 
   * output_tsv_filename: Path to the output TSV file for storing merged spectral information

Run the command:
```
python3 toprepo/src/process/tsv/merge_mzml_msalign_info.py spectra_mzml_info.tsv spectra_msalign_feature_info.tsv spectra_mzml_msalign_feature_info.tsv
```


**1.7 Merge all spectral information**

This step merges the spectral information obtained from mzML, msalign, and feature files with the spectral identification results generated by TopPIC.

* Usage python3 merge_mzml_msalign_toppic_info.py <toppic_info_filename> <mzml_msalign_feature_info_filename> <output_tsv_filename>

* Input:
  * toppic_info_filename: A preprocessed TSV file containing spectral identification results reported by TopPIC 
  * mzml_msalign_feature_info_filename: A TSV file containing merged spectral information extracted from mzML, msalign, and feature files 
  * output_tsv_filename: A TSV output file for storing merged spectral information extracted from mzML, msalign, feature, and TopPIC identification files 

Run the command: 
```
python3 toprepo/src/process/tsv/merge_mzml_msalign_toppic_info.py spectra_toppic_info.tsv spectra_mzml_msalign_feature_info.tsv  spectra_mzml_msalign_feature_toppic_info.tsv
```

## 2. Generate annotated msalign files 

We use an msalign file spectra_ms2.msalign with dataset id PXD029703 and its spectral information file spectra_mzml_msalign_feature_toppic_info.tsv to explain the method.

**2.1 Preprocess msalign file**
This step adds dataset id information to each scan in the msalign file and updates the spectral information in the msalign file.

```
python3 toprepo/src/process/msalign_anno/msalign_preprocess.py spectra_ms2.msalign PXD029703 spectra_preprocess_ms2.msalign
```

**2.2 Add PrSM identification information to msalign file**
This step adds spectral identification information to the msalign file. 
```
python3 toprepo/src/process/msalign_anno/merge_msalign_prsm.py --tsv spectra_mzml_msalign_feature_toppic_info.tsv --msalign spectra_preprocess_ms2.msalign --out spectra_prsm_ms2.msalign
```

**2.3 Annotate msalign file** 
This step adds annotations to the msalign file.  
```
python3 toprepo/src/process/msalign_anno/msalign_anno.py --msalign spectra_prsm_ms2.msalign --out spectra_anno_ms2.msalign
```


## 3. Generate annotated mgf files.

We use an example mzML file spectra.mzML and its corresponding annotated msalign file spectra_anno_ms2.msalign to explain the method.

**3.1 Convert matched mzML files to MGF file**  
The Python script convert_mzml_mgf.py converts matched mzML files into MGF format. Two modes are supported: 1) Standard MGF; 2) Extended MGF (including additional fields such as collision energy, scan limits, total ion current, etc).

* Usage python3 convert_mzml_to_mgf.py <input_mzml_filename> <output_mgf_filename>

Run the command:

```
python3 toprepo/src/process/mzml/convert_mzml_to_mgf.py spectra.mzML spectra_ms2.mgf 
```

**3.2 Add dataset id to MGF file**  

```
python3 toprepo/src/process/mgf/mgf_add_dataset_id.py spectra_ms2.mgf PXD029703 spectra_dataset_id_ms2.mgf
```

**3.3 Annotation**
The Python script mgf_anno_multi_thread.py annotates MGF files by matching experimental centroid peaks to theoretical fragment peaks and theoretical ions.

* Usage: python3 mgf_anno_file.py --theo_file <theo_file> --mgf_file <input_mgf_file> --msalign_file <input_msalign_file> --out <output_annotated mgf_file> 

Input:
* A theoretical_envelope_file: A txt file containing theoretical isotope envelope distributions. "toprepo/resources/theo_patt.txt".
* An MGF file: ```spectra_ms2.mgf``` 
* An annotated msalign file: ```spectra_anno_ms2.msalign```

Output:
* An annotated mgf file

Run the command:
```
python3 toprepo/src/process/mgf/mgf_anno_file.py --theo_file toprepo/resources/theo_patt.txt --mgf_file spectra_dataset_id_ms2.mgf --msalign_file spectra_anno_ms2.msalign --out spectra_anno_ms2.mgf
```
