
import argparse
import glob
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Organise a DICOM dataset in BIDS. '
        'Folders should be organized in one folder with one folder '
        'by subject with inside one folder by session'
    )
    parser.add_argument(
        '-s', '--sourcedata', required=True, help='directory with DICOM folders'
    )
    parser.add_argument(
        '-o', '--output', required=True, help='output directory'
    )
    parser.add_argument(
        '-c', '--config_file', required=True, help='dcm2bids config file'
    )
    args = parser.parse_args()
    directory = args.sourcedata
    output_directory = args.output

    config_file = args.config_file
    all_dicom_folders = glob.glob(os.path.join(directory, '*'))


    for folder in all_dicom_folders:
        print(f"\n {folder}")
        if os.path.isdir(folder):
            id_sub = folder.split('/')[-1].split('_')[-2]
            session = folder.split('/')[-1].split('_')[-1]

            print('\nSubject ', id_sub)
            print('Session ', session)
            
            cmd = 'dcm2bids -d ' + folder + ' -p ' + id_sub + \
                ' -s' + session + ' -c ' + config_file + ' -o ' + output_directory
            print('cmd: ', cmd)
            os.system(cmd)

