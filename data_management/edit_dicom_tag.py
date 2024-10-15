'''
Edit a specific DICOM tag (change the value)
'''
import argparse
import glob
import os
import pydicom
from pydicom.tag import Tag


def find_dicom_tag_value(input_dicom_dataset, tag):
    '''
    Find DICOM tag value
    '''
    try:
        value = input_dicom_dataset[tag].value
    except ValueError:
        value = 'ValueNotFound'
    return value


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Edit a specific DICOM tag (change the value)'
    )
    parser.add_argument(
        '-d', '--directory', required=True, help='input directory'
    )
    parser.add_argument(
        '-o', '--output', required=True, help='output directory'
    )
    parser.add_argument(
        '-t', '--tag', required=True, help='tag to modified'
    )
    parser.add_argument(
        '-v', '--value', required=True, help='new value'
    )
    args = parser.parse_args()
    directory = args.directory
    output_directory = args.output
    tag = str(args.tag)
    new_value= args.value

    # Get all DICOM in directory
    input_dicom_files = glob.glob(os.path.join(directory, '*'))

    sequences = []
    for dicom_file in input_dicom_files:
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)
        new_dicom = os.path.join(output_directory, os.path.basename(dicom_file))
        ds =  pydicom.read_file(dicom_file)
        ds[Tag(tag)].value = new_value
        ds.save_as(new_dicom)

