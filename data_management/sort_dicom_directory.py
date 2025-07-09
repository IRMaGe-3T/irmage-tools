"""
Sort DICOM file using Patient's Name, 
Study date and Series Description tag
"""
import argparse
import os
import pydicom
import unidecode
import shutil

from pydicom.tag import Tag


def find_dicom_tag_value(input_dicom_dataset, tag):
    """
    Find DICOM tag value
    """
    try:
        value = input_dicom_dataset[tag].value
    except ValueError:
        value = "ValueNotFound"
    return value


def get_all_dicom_files(dicom_directory):
    """
    Get all DICOM files from a directory
    """
    if len(dicom_directory) == 0:
        message = "Input folder is empty."
        raise RuntimeError(message)

    raw_storage_methods_not_taken = [
        "1.2.840.10008.5.1.4.1.1.11.1",
        "1.2.840.10008.5.1.4.1.1.66",
        "Secondary Capture Image Storage",
        "1.2.840.10008.5.1.4.1.1.7",
    ]

    raw_input_files = os.listdir(dicom_directory)
    input_files = []

    # Check if each file is a DICOM
    for file_name in raw_input_files:
        file_start = file_name.split("/")[-1].split("_")[0].lower()
        file_ext = file_name.split(".")[-1]
        not_taken_start = [
            "ps",
            "xx",
            "dicomdir"
        ]
        not_taken_ext = [
            "bvecs",
            "bvals",
            "txt",
            'PAR', 
            'REC',
            'nii',
            'xml'
        ]

        if file_start in not_taken_start or file_ext in not_taken_ext:
            print(f"File {file_name} not taken (not a DICOM image).")
            continue
        elif os.path.isdir(os.path.join(dicom_directory, file_name)):
            # Add all files from subdirectory to file list
            sub_directory = os.path.join(dicom_directory, file_name)
            for file_in_directory in os.listdir(sub_directory):
                file_path = os.path.join(file_name, file_in_directory)
                raw_input_files.append(file_path)
            continue

        try:
            input_dicom = pydicom.read_file(
                os.path.join(dicom_directory, file_name)
            )
        except Exception:
            print(f"File {file_name} not taken (coul not read DICOM).")
            continue
        storage_method = find_dicom_tag_value(input_dicom, Tag(0x08, 0x16))[0]
        if storage_method in raw_storage_methods_not_taken:
            print(f"File {file_name} not taken (bad storage method).")
            continue
        input_files.append(os.path.join(dicom_directory, file_name))
    return input_files


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Sort DICOM directory  using Patient s Name,'
        'Study date and Series Description tag'
    )
    parser.add_argument(
        '-d', '--directory', required=True, help='input directory'
    )
    parser.add_argument(
        '-o', '--output', required=True, help='output directory'
    )
    args = parser.parse_args()
    directory = args.directory
    output_directory = args.output

    print(f'Sort DICOM files from {directory}')

    # Get all DICOM in directory
    input_dicom_files = get_all_dicom_files(directory)

    sequences = []
    for dicom_file in input_dicom_files:
        # Get info as patient name , study date , serie description
        dataset = pydicom.read_file(dicom_file)
        patient_name_value = find_dicom_tag_value(dataset, Tag(0x10, 0x0010))
        patient_name = "".join(
            filter(str.isalnum, unidecode.unidecode(patient_name_value))
        )
        study_date = find_dicom_tag_value(
            dataset, Tag(0x08, 0x0020)
        )
        series_description = find_dicom_tag_value(
            dataset, Tag(0x08, 0x103e)
        )
        serie_num = str(find_dicom_tag_value(
            dataset, Tag(0x20, 0x0011)
        ))

        sequence = patient_name + '_' + study_date + '_' + serie_num + '_' + series_description

        if sequence not in sequences:
            sequences.append(sequence)
        new_path = os.path.join(output_directory, sequence)
        if not os.path.exists(new_path):
            os.makedirs(new_path)
        # Copy dicom into sequence folder
        new_dicom_path = os.path.join(new_path, dicom_file.split('/')[-1])
        i = 1
        while os.path.exists(new_dicom_path):
            new_dicom_path = os.path.join(new_path, dicom_file.split('/')[-1] + '_' + str(i))
            i += 1

        shutil.copy2(dicom_file, new_dicom_path)

    print(f'\nList of the sequneces detected: {sequences}')
    print(f'\nDICOM files copied to {output_directory}')
