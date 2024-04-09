""" 
Compute SliceTiming metadata for json BIDS (Philips system)

Ressources: 

- https://neurostars.org/t/warning-1-code-13-slice-timing-not-defined/5273/3 
- https://neurostars.org/t/heudiconv-no-extraction-of-slice-timing-data-based-on-philips-dicoms/2201/12 
- https://mblab.si/blog/empirical-slice-timing-order-of-philips-achieva-multiband-bold-images/


"""

import argparse
import math


def compute_bids_slice_timing(tr_sec,
                              numb_slices,
                              slice_order,
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
    if slice_order == 'ascending':
        bids_slice_timing = slice_timing
    # bids slice timing for descending acquistion
    elif slice_order == 'descending':
        bids_slice_timing = slice_timing[::-1]
        order = order[::-1]
    # bids slice timing for default acquistion
    elif slice_order == 'default':
        order_odd = list(range(1, numb_slices + 1, 2))
        order_even = list(range(2, numb_slices + 1, 2))
        order = order_odd + order_even
        bids_slice_timing =  [x for _,x in sorted(zip(order, slice_timing))]
    # bids slice timing for interleaved acquistion
    elif slice_order == 'interleaved':
        increment = round(math.sqrt(numb_slices))
        order = []
        for i in list(range(1, increment + 1)):
            order += list(range(i, numb_slices + 1, increment))
        bids_slice_timing = [x for _,x in sorted(zip(order, slice_timing))]
    else:
        print(f'Slice order {slice_order} not implemented yet')
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
    slice_order = args.order.lower()

    if slice_order not in ['ascending',
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
            slice_order=slice_order,
            multiband_factor=int(args.multiband),
            delay_between_vol_sec=int(args.delay_between_vol_sec),
            package=int(args.package)
        )

        print('Slice order: \n', slice_order)
        print('\n')
        for i in list(range(len(slice_order))):
            print(f'Slice : {i + 1} / Timing: {bids_slice_timing[i]}')
