# irmage-tools
Project to share useful codes and snippets or jupyter notebooks for neuroimaging. 

Please note that thoses codes are proposed as example and should be modified according to need.


## Table of contents
1. [data_management](#data_management)
    - [philips_compute_slice_timing.py](#philips_compute_slice_timing.py)
    - [philips_compute_slice_timing_dataset_bids.py](#philips_compute_slice_timing_dataset_bids.py)
    - [sort_dicom_directory.py](#sort_dicom_directory.py)
    - [bids_data_with_dcm2bids.py](#bids_data_with_dcm2bids.py)
    - [edit_dicom_tag.py](#edit_dicom_tag.py)
    - [mri_json_to_xls.py](#mri_json_to_xls.py)
2. [physiological_data](#physiological_data)
    - [HRV_philips_scanphyslog](#HRV_philips_scanphyslog)


## data_management 

Useful codes for data managment (read DICOM, sort DICOM, get info...)

### philips_compute_slice_timing.py <a name="philips_compute_slice_timing.py"></a>
Compute SliceTiming metadata for json BIDS (Philips system)

Example: 

```
python philips_compute_slice_timing.py -tr 2.5 -slices 42 -order ascending
```

### philips_compute_slice_timing_dataset_bids.py
Compute SliceTiming metadata for json BIDS (Philips system) and add information to a BIDS dataset.

Example: add SliceTiming for all the resting state (task rest) functional MRI data of the BIDS datatset study: 

```
python philips_compute_slice_timing_dataset_bids.py -bids_path '/home/username/study' -task rest -tr 2.5 -slices 42 -order ascending
```

### sort_dicom_directory.py
Sort DICOM files in a directory using Patient's Name, Study date and Series Description tag

Example: 

```
python sort_dicom_directory.py -d /home/username/dicom_directory_to_sort -o /home/username/dicom_directory_sorted
```

### bids_data_with_dcm2bids.py
Run dcm2bids on several subjects to organised data in BIDS.  

The sourcedata folder should be organised as follow (where "subject1" / "subject2" are the desired name in the BIDS architecture and "V1" / "V2" are the sessions names): 
```
├──DICOM_sourcedata/
|   ├──subject1_V1/
|   |   ├──DICOM/
|   ├──subject1_V2/
|   |   ├──DICOM/
|   ├──subject2_V1/
|   |   ├──DICOM/
|   ├──subject2_V2/
|   |   ├──DICOM/
```

Example: 

```
    python bids_data_with_dcm2bids.py -s /DICOM_sourcedata/ -o /output_bids_direcory/ -c config_file_dcm2bids.jsob
```

### edit_dicom_tag.py
Edit a specific DICOM tag (change the value). 
All the DICOM of the directory will be copied in the output directory and the tag value will be edit in the DICOM of the output directory.

The tag (-t) should be an str corresponding the element's keyword in the DICOM (as "PatientName", see [pydicom](https://pydicom.github.io/pydicom/dev/reference/generated/pydicom.tag.Tag.html))

Example:
```
    python edit_dicom_tag.py -d /input_DICOM_directory/ -o /output_DICOM_direcory/ -t ImagesInAcquisition -v 1752
```


### mri_json_to_xls.py {#mri_json_to_xls.py}
Obtain an excel file with all the value of the MRI json. 
The input folder could be a BIDS folder (--bids) or a folder with NIfTI/json (--path). 

It is possible to specify the output directory (--ouput_directory) and the output file name (--SCANPHYSLOGxxx)

Example:
```
    mri_json_to_xlsx.py -b path/to/bids/folder -f study.xlsx
```

## Physiological_data

### HRV_philips_scanphyslog

#### ppu4fmri.m
Read SCANPHYSLOG file (from Philips system) containing PPU data (PPU signal reflects the heart rate at fingertips), computing the tachogram (HRV) and extracting and resampling the LF-HRV to be used as regressors in fMRI (task or rs) analysis.

You need spm in your Matlab path to run this program.

Launch the program in matlab, select the SCANPHYSLOG file (SCANPHYSLOGxxx.log) and enter acquisition parameters. 4 files are obtained:
- SCANPHYSLOGxxx_QC.png: plots for peaks of PPU , tachogram, LF and HR component
- SCANPHYSLOGxxx_frequency_spectrum.png: frequency spectrum
- SCANPHYSLOGxxx_results.tsv: tachogram analysis results (heart rate, RMSSD, SDNN, total power ..)
- SCANPHYSLOGxxx_reg_LF_HRV.tsv: LF-regressors for fMRI 

Examples of SCANPHYSLOG files can be found [here](physiological_data/HRV_philips_scanphyslog/examples)


