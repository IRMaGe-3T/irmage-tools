"""
Add / remove / update json key for bids folder
"""
import argparse
import json
import os
import glob

def load_json(path):
    """Load a JSON file safely and return a dictionary."""
    if not os.path.exists(path):
        print("Warning: JSON file does not exist. Creating a new one.")
        return {}

    try:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    except json.JSONDecodeError:
        print("Error: JSON file is corrupted or invalid. Using an empty dictionary.")
        return {}
    except Exception as e:
        print(f"Error loading JSON file: {e}")
        return {}

def save_json(path, data):
    """Save the dictionary into the JSON file safely."""
    try:
        with open(path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=4, ensure_ascii=False)
    except Exception as e:
        print(f"Error saving JSON file: {e}")

def validate_key(key):
    """Ensure the key is a valid non-empty string."""
    if not isinstance(key, str) or key.strip() == "":
        print("Error: Key must be a non-empty string.")
        return False
    return True

def add_key(path, key, value):
    if not validate_key(key):
        return

    data = load_json(path)

    if key in data:
        print(f"Error: Key '{key}' already exists.")
        return

    data[key] = value
    save_json(path, data)
    print(f"Key '{key}' added with value: {value}")

def update_key(path, key, new_value):
    if not validate_key(key):
        return

    data = load_json(path)

    if key not in data:
        print(f"Error: Key '{key}' does not exist, cannot update.")
        return

    data[key] = new_value
    save_json(path, data)
    print(f"Key '{key}' updated to: {new_value}")

def delete_key(path, key):
    if not validate_key(key):
        return

    data = load_json(path)

    if key not in data:
        print(f"Error: Key '{key}' does not exist, cannot delete.")
        return

    del data[key]
    save_json(path, data)
    print(f"Key '{key}' deleted.")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="bids_update_json",
        description="Update/Remove/add a key in json files from BIDS folder (using acq and modality)",
        epilog="Example "
        "bids_update_json.py -b path/to/bids/folder --delete key1"
    )

    parser.add_argument(
        "-b", "--bids",
        help="Path to BIDS dataset",
        required=True
    )
    parser.add_argument(
        "-d", "--datatype",
        help="BIDS data type (anat, dwi, fmap...)",
        required=True
    )
    parser.add_argument(
        "-m", "--modality",
        help="BIDS modality (T1w, dwi, epi..)",
        required=True
    )
    parser.add_argument(
        "-a", "--acq",
        help="BIDS acq ",
        required=True
    )
    parser.add_argument(
        "--dir",
        help="BIDS dir (PA or AP) ",
        required=False
    )
    parser.add_argument(
        "--delete",
        help="Key to delete (ex : PatientName)",
        required=False
    )
    parser.add_argument(
        "--add",
        action = type('', (argparse.Action, ), dict(__call__ = lambda a, p, n, v, o: getattr(n, a.dest).update(dict([v.split('=')])))), 
        default = {},
        help="values to add (ex: --add 'PhaseEncodingDirection'='j-' --add 'TaskName'='rs')",
        required=False
    )

    parser.add_argument(
        "--update",
        action = type('', (argparse.Action, ), dict(__call__ = lambda a, p, n, v, o: getattr(n, a.dest).update(dict([v.split('=')])))), 
        default = {},
        help="values to update in a dictionnary  (ex: --update 'PhaseEncodingDirection'='j-' --update 'TaskName'='rs')",
        required=False
    )
    args = parser.parse_args()

    if not args.add and not args.delete and not args.update:
        parser.error("--add or --delete or --update should be used")

    if args.dir:
        files = glob.glob(os.path.join(
                args.bids, "sub-*", "ses-*", args.datatype, "*_acq-" + args.acq +  "_dir-" + args.dir +"*_" + args.modality + ".json"))
        files += glob.glob(os.path.join(args.bids, "sub-*",  args.datatype, "*_acq-" + args.acq +   "_dir-" + args.dir + "*_" + args.modality + ".json"))
    else:
        files = glob.glob(os.path.join(
                args.bids, "sub-*", "ses-*", args.datatype, "*_acq-" + args.acq + "*_" + args.modality + ".json"))
        files += glob.glob(os.path.join(args.bids, "sub-*",  args.datatype, "*_acq-" + args.acq + "*_" + args.modality + ".json"))

    for json_file in files:
        print("\njson file: ", json_file)
        if args.delete:
            delete_key(json_file, args.delete)
        if args.add:
            for key in args.add.keys():
                add_key(json_file, key, args.add[key])
        if args.update:
            for key in args.update.keys():
                update_key(json_file, key, args.update[key])
