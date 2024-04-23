''' 
Compute SliceTiming metadata for json BIDS (Philips system) 
and add information to a BIDS dataset

Ressources: 
- https://en.wikibooks.org/wiki/SPM/Slice_Timing 
- https://neurostars.org/t/warning-1-code-13-slice-timing-not-defined/5273/3 
- https://neurostars.org/t/heudiconv-no-extraction-of-slice-timing-data-based-on-philips-dicoms/2201/12 
- https://mblab.si/blog/empirical-slice-timing-order-of-philips-achieva-multiband-bold-images/

'''

import argparse
import json
import math
import os

from bids import BIDSLayout


def compute_bids_slice_timing(tr_sec,
                              numb_slices,
                              slice_scan_order,
                              multiband_factor=1,
                              delay_between_vol_sec=0,
                              package=1):
    '''
        Compute SliceTiming metadata for json BIDS 
        https://bids-specification.readthedocs.io/en/stable/glossary.html#slicetiming-metadata

        SliceTiming is the time at which each slice was acquired within each 
        volume (frame) of the acquisition.

        tr_sec (float): repetion time in seconds
        numb_slices (int): number of slice 
        slice_order (either  ascending, descending, default, interleaved): slice order
        multiband_factor (int): multiband_factor (1 = no multiband) 
        delay_between_vol_sec (int) : pause between final slice of volume and 
                                start of next in seconds (sparse)
        package (int): nimber of package
    
    '''

    if multiband_factor != 1:
        print('Multibande option is not implemented yet')
        return
    if package !=1:
        print('Not implemented for several packages')
        return
    # acquistion time of one slice
    ta = (tr_sec - delay_between_vol_sec) / numb_slices

    # bids slice timing for ascending acquistion
    slice_timing = [ta * i for i in range(numb_slices)]
    order = list(range(1, numb_slices + 1))
    if slice_scan_order == 'ascending':
        bids_slice_timing = slice_timing
    # bids slice timing for descending acquistion
    elif slice_scan_order == 'descending':
        bids_slice_timing = slice_timing[::-1]
        order = order[::-1]
    # bids slice timing for default acquistion
    elif slice_scan_order == 'default':
        order_odd = list(range(1, numb_slices + 1, 2))
        order_even = list(range(2, numb_slices + 1, 2))
        order = order_odd + order_even
        bids_slice_timing =  [x for _,x in sorted(zip(order, slice_timing))]
    # bids slice timing for interleaved acquistion
    elif slice_scan_order == 'interleaved':
        increment = round(math.sqrt(numb_slices))
        order = []
        for i in list(range(1, increment + 1)):
            order += list(range(i, numb_slices + 1, increment))
        bids_slice_timing = [x for _,x in sorted(zip(order, slice_timing))]
    else:
        print(f'Slice order {slice_scan_order} not implemented yet')
        return

    return order, bids_slice_timing


def update_json(json_path, new_data):
    '''
        Update a json file
        json_path (string): json path
        new_data (dict): data to add to the json

    '''
    with open(json_path, 'r', encoding='utf-8') as json_file:
        file_data = json.load(json_file)
    # Check if data already in json
    for k in list(new_data.keys()):
        if k in list(file_data.keys()):
            print(f"\nKey already in json file for {json_path}")
            return 0
    file_data.update(new_data)
    with open(json_path, 'w', encoding='utf-8') as json_file:
        json.dump(file_data, json_file, indent=4)
    return 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compute SliceTiming metadata BIDS and add it to json'
    )
    parser.add_argument('-bids_path', required=True, help='bids dataset path')
    parser.add_argument(
        '-task', required=True, help='task of the functional MRI to update'
    )
    parser.add_argument(
        '-tr', required=True, help='repetion time in seconds'
    )
    parser.add_argument('-slices', required=True, help='number of slice')
    parser.add_argument(
        '-order',
        required=False,
        default='ascending',
        help='slice order (either ascending, descending, default or interleaved)'
    )
    parser.add_argument(
        '-multiband',
        required=False,
        default=1,
        help='multiband_factor (1 = no multiband)'
    )
    parser.add_argument(
        '-delay_between_vol_sec',
        required=False,
        default=0,
        help='pause between final slice of volume and start of next in seconds (sparse data)'
    )
    parser.add_argument(
        '-package',
        required=False,
        default=1,
        help='number of package'
    )

    args = parser.parse_args()
    slice_scan_order = args.order.lower()

    if slice_scan_order not in ['ascending',
                           'descending',
                           'default',
                           'interleaved']:
        print('Slice order should either ascending,'
              'descending, default or interleaved')
    else:
        # Compute slice order and SliceTiming
        slice_order, bids_slice_timing = compute_bids_slice_timing(
            float(args.tr),
            int(args.slices),
            slice_scan_order=slice_scan_order,
            multiband_factor=int(args.multiband),
            delay_between_vol_sec=float(args.delay_between_vol_sec),
            package=int(args.package)
        )

        print('Slice order: \n', slice_order)
        print('\n')
        for i in list(range(len(slice_order))):
            print(f'Slice : {i + 1} / Timing: {bids_slice_timing[i]}')

        slice_timing_data = {
            'SliceTiming': bids_slice_timing
        }

        # Add SliceTiming to BIDS json
        bids_directory_path = args.bids_path
        layout = BIDSLayout(bids_directory_path)
        subjects = layout.get_subjects()
        sessions = layout.get_sessions()
        task = args.task

        print(
            f'\nAdd SliceTiming metadata for task {task} '
            f'functional MRI in dataset {bids_directory_path}'
        )
        done = []
        not_done = []
        for sub in subjects:
            for sess in sessions:
                f = layout.get(subject=sub, session=sess,
                            task=task, extension='nii.gz')
                if len(f) == 1:
                    json_file = os.path.join(
                        bids_directory_path,
                        'sub-' + sub,
                        'ses-' + sess,
                        'func',
                        layout.get(subject=sub, session=sess, task=task,
                                extension='json')[0].filename
                    )
                    update = update_json(json_file, slice_timing_data)
                    if update == 1:
                        done.append(sub + '_' + sess)
                    else:
                        not_done.append(sub + '_' + sess)
                else:
                    not_done.append(sub + '_' + sess)
        print(f'\nMetadata added for : {done}')
        print(f'\n Metadata NOT added for : {not_done} (error or no task for this subject)')
