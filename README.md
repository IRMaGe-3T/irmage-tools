# irmage-tools
Project to share useful codes and snippets or jupyter notebooks for neuroimaging. 

Please note that thoses codes are proposed as example and should be modified according to need.

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

Example: Add SliceTiming for all the resting state (task rest) functional MRI data of the BIDS datatset study: 

```
python philips_compute_slice_timing_dataset_bids.py -bids_path '/home/username/study' -task rest -tr 2.5 -slices 42 -order ascending
```


