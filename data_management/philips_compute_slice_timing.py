""" 
Compute SliceTiming metadata for json BIDS (Philips system)

/!\ For multiband data it is experimental. Please check the validity for your data

Ressources: 

- https://neurostars.org/t/warning-1-code-13-slice-timing-not-defined/5273/3 
- https://neurostars.org/t/heudiconv-no-extraction-of-slice-timing-data-based-on-philips-dicoms/2201/12 
- https://mblab.si/blog/empirical-slice-timing-order-of-philips-achieva-multiband-bold-images/


"""

import argparse
import math
import numpy as np


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

    if package != 1:
        print('Not implemented for several packages')
        return

    # If multibande
    if multiband_factor != 1:
        print('For data with multiband, the slice scan order computed is experimental. '
              'The BIDS SliceTiming is not given because it is not possible to be sure abour it.'
              'Please check with the MR physicist the values')
        # Check that number of slice / MB > 9
        test = numb_slices / multiband_factor
        if test < 9:
            print('Implemented only for numb_slices / multiband_factor>9')
            return
        slices_per_band = int(numb_slices / multiband_factor)
        order_temp = np.zeros([multiband_factor, slices_per_band], dtype='int')
        group_temp = np.zeros(slices_per_band, dtype='int')

        if slice_scan_order == 'ascending':
            for k in range(multiband_factor):
                order_temp[k, :] = range(
                    slices_per_band * k + 1, slices_per_band * (k + 1) + 1, 1)
        elif slice_scan_order == 'descending':
            for k in range(multiband_factor):
                order_temp[k, :] = range(
                    slices_per_band * (k + 1), slices_per_band * k, -1)
        elif slice_scan_order == 'default':
            k = 0
            step = 2
            if math.fmod(slices_per_band, 2) == 0:
                # half_ceil = math.ceil((slices_per_band+1)/step)
                half_locs_per_package = math.ceil((slices_per_band + 1)/step)
                half_floor = math.floor((slices_per_band + 1)/step)
                half_locs_per_package = half_floor
                low_part_start_loc = 0
                low_part_act_loc = 0
                high_part_start_loc = slices_per_band - 1
                high_part_act_loc = high_part_start_loc

                for i in range(slices_per_band + 1):
                    if low_part_act_loc < half_locs_per_package:
                        group_temp[k] = low_part_act_loc
                        low_part_act_loc = low_part_act_loc+step
                        k = k + 1
                    elif high_part_act_loc >= half_locs_per_package:
                        group_temp[k] = high_part_act_loc
                        high_part_act_loc = high_part_act_loc-step
                        k = k + 1
                    else:
                        low_part_act_loc = low_part_start_loc + 1
                        high_part_act_loc = high_part_start_loc - 1

                for k in range(multiband_factor):
                    order_temp[k, :] = group_temp + 1 + k * slices_per_band
        elif slice_scan_order == 'interleaved':
            step = round(math.sqrt(slices_per_band))

            slice_group = 1
            group_temp[0] = 1

            for l in range(slices_per_band - 1):
                current = group_temp[l] + step

                if current > slices_per_band:
                    slice_group = slice_group + 1
                    current = slice_group

                group_temp[l + 1] = current

                for k in range(multiband_factor):
                    order_temp[k, :] = group_temp + k * slices_per_band
        else:
            print(f'Slice order {slice_scan_order} not implemented yet')
            return
        order = []
        for i in range(slices_per_band):
            order.append(list(order_temp[:, i]))

        # Not possible to be sure for SliceTiming for MB data
        # ta = (tr_sec - delay_between_vol_sec) / slices_per_band
        # slice_timing = [ta * i for i in range(slices_per_band)]
        bids_slice_timing = None

    else:
        # Without multiband
        # ta = acquistion time of one slice
        ta = (tr_sec - delay_between_vol_sec) / numb_slices
        slice_timing = [ta * i for i in range(numb_slices)]
        order = list(range(1, numb_slices + 1))
        if slice_scan_order == 'ascending':
            bids_slice_timing = slice_timing
        elif slice_scan_order == 'descending':
            bids_slice_timing = slice_timing[::-1]
            order = order[::-1]
        elif slice_scan_order == 'default':
            order_odd = list(range(1, numb_slices + 1, 2))
            order_even = list(range(2, numb_slices + 1, 2))
            order = order_odd + order_even
            bids_slice_timing = [
                x for _, x in sorted(zip(order, slice_timing))]
        elif slice_scan_order == 'interleaved':
            increment = round(math.sqrt(numb_slices))
            order = []
            for i in list(range(1, increment + 1)):
                order += list(range(i, numb_slices + 1, increment))
            bids_slice_timing = [
                x for _, x in sorted(zip(order, slice_timing))]
        else:
            print(f'Slice order {slice_scan_order} not implemented yet')
            return

    return order, bids_slice_timing


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compute SliceTiming metadata BIDS'
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
            delay_between_vol_sec=int(args.delay_between_vol_sec),
            package=int(args.package)
        )

        print('Slice order: \n', slice_order)
        if args.multiband == 1:
            print('\n')
            for i in list(range(len(slice_order))):
                print(f'Slice : {i + 1} / Timing: {bids_slice_timing[i]}')
        else:
            print('SliceTiming not computed for multiband data')
