# irmage-tools
Project to share useful codes and snippets or jupyter notebooks for neuroimaging. 

Please note that thoses codes are proposed as example and should be modified according to need.


Table of contents

[data_management ](#data_management)

<a name="data_management"></a>
## data_management 

Useful codes for data managment (read DICOM, sort DICOM, get info...)

### philips_compute_slice_timing.py 
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

The sourcedata folder should be organised as follow: 
```
├──DICOM_sourcedata/
|   ├──subject1/
|   |   ├──01/
|   |   |   ├──DICOM/
|   |   ├──02/
|   |   |   ├──DICOM/
|   ├──subject2/
|   |   ├──01/
|   |   |   ├──DICOM/
|   |   ├──02/
|   |   |   ├──DICOM/
```

Example: 

```
    python bids_data_with_dcm2bids.py -s /DICOM_sourcedata/ -o /output_bids_direcory/ -c cobfig_file_dcm2bids.jsob
```


