"""
Obtain an excel file with all the value of the MRI json 
(from a BIDS folder or from a folder with NIfTI/json)
"""
import argparse
import glob
import json
import os
import pandas as pd


def get_value_json(json_path, field_name):
    """
    Get a value from a json for one filed

    json_path (string): json path
    field_name (string): field name (example: "Manufacturer")
    """
    try:
        field = json_path[field_name]
    except Exception:
        field = "na"
    return field


def get_all_labels(files_path):
    """
    Obtain all the possible fields of the json files 

    files_path (a list): list of the json path
    """
    labels = ["File"]
    for path in files_path:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
            for key, value in data.items():
                if key not in labels:
                    labels.append(key)
    return labels


def get_values(files_path, labels):
    """
    Obtain all the values for the selected labels 
    from the json files

    files_path (a list): list of the json path
    labels (a list): list of the label
    """
    dataframe = pd.DataFrame(columns=labels)
    for path in files_path:
        temp = [os.path.basename(path)]
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
            for label in labels:
                if label != "File":
                    temp.append(get_value_json(data, label))
        dataframe.loc[len(dataframe)] = temp
    return dataframe


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="mri_json_to_xlsx",
        description="Create an xlsx file from multiple json MRI files.",
        epilog="Example "
        "mri_json_to_xlsx.py -b path/to/bids/folder -f study.xlsx"
    )

    parser.add_argument(
        "-p", "--path",
        help="Path to a folder with nii/json (exclusive with --bids",
        required=False
    )
    parser.add_argument(
        "-b", "--bids",
        help="Path to BIDS dataset (exclusive with --path)",
        required=False
    )
    parser.add_argument(
        "-f", "--filename",
        help="Output filename (default = params_sequence.xlsx)",
        required=False
    )
    parser.add_argument(
        "-o", "--ouput_directory",
        help="output directory (default save to input folder)",
        required=False
    )
    args = parser.parse_args()

    if args.path and args.bids:
        parser.error("--path and --bids are exclusive")

    if args.path:
        files = glob.glob(os.path.join(args.path, "*.json"))
    elif args.bids:
        files = glob.glob(os.path.join(
            args.bids, "sub-*", "ses-*", "*", "*.json"))
        files += glob.glob(os.path.join(args.bids, "sub-*", "*", "*.json"))

    if args.ouput_directory:
        output_path = args.ouput_directory
    else:
        if args.path:
            output_path = args.path
        elif args.bids:
            output_path = args.bids

    if args.filename:
        filename = args.filename
        if not ".xlsx" in filename:
            filename += ".xlsx"
    else:
       filename = "params_sequence.xlsx"

    labels_columns = get_all_labels(files)
    df = get_values(files, labels_columns)
    df.to_excel(os.path.join(output_path, filename))
