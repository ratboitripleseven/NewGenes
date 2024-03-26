import os
import pandas as pd
import argparse
import urllib
import gzip
import shutil
from data_loader.NCBI_data_loader import NCBIDataDownloaderPrep

ROOT_FOLDER = 'partition_file/'
SEQUENCES_FOLDER = 'data/NCBI/sequence_files'
PREPPED_SEQUENCES_FOLDER = 'data/NCBI/prep/'


'''
NOTE:

This script is create to download and set up the files after creating a partition file!
partition file has the dl link of the file!

It seems that partition file neglects if file exist or not... 
This will be quite annoying temp solution perhaps: create a new partition file?

'''
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f',
        '--file',
        type = str,
        default = 'phylum_Fibrobacterota_test',
        help = 'Choose taxonomy level'
    )
    return parser.parse_args()

def assert_args(args):
    # check if ROOT_FOLDER exists
    if not os.path.isdir(ROOT_FOLDER):
        raise NotADirectoryError(f'{ROOT_FOLDER} is not found! please use create a partition first')

    file_path = ROOT_FOLDER+args.file+'.csv'
    if not os.path.isfile(file_path):
        raise NotADirectoryError(f'{file_path} is not found! create a partition!')

    return file_path

def download_file(asm_name, link):
    # Taken from NCBI_DATA_LOADER
    if not os.path.isfile(SEQUENCES_FOLDER+f'/{asm_name}.fna'):
        #download link
        print(f'downloading {asm_name}, at {link} ')
        urllib.request.urlretrieve(link, SEQUENCES_FOLDER+f'/{asm_name}.fna.gz')
        with gzip.open(SEQUENCES_FOLDER+f'/{asm_name}.fna.gz', 'rb') as f_in:
            with open(SEQUENCES_FOLDER+f'/{asm_name}.fna', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(SEQUENCES_FOLDER+f'/{asm_name}.fna.gz')
    else:
        print(f"File: {SEQUENCES_FOLDER+'/'+asm_name+'.fna'} exists")
    
    
def create_csv(asm_name, link):
    ncbi_object = NCBIDataDownloaderPrep(asm_name,link)
    ncbi_object.to_HGTDB('csv')
        
def set_partition_file(file_path):
    dataframe = pd.read_csv(file_path)
    print('Downloading the files ')
    for i in range(len(dataframe)):
        try:
            download_file(dataframe.loc[i,'asm_name'], dataframe.loc[i,'dl_link'])
        except:
            print(f"Unable to download {dataframe.loc[i,'asm_name']}, {dataframe.loc[i,'dl_link']} ")
    print('Files are downloaded')
    
    list_of_unavailable_csv_files = []
    print('Creating csvs for files')
    for i in range(len(dataframe)):
        try:
            create_csv(dataframe.loc[i,'asm_name'], dataframe.loc[i,'dl_link'])
        except:
            list_of_unavailable_csv_files.append(i)
            print(f"Unable to create csv for {dataframe.loc[i,'asm_name']} ")
            
    print('#'*80)
    new_file_path = file_path.replace('.csv', '_available.csv')
    new_dataframe = dataframe.drop(index=list_of_unavailable_csv_files)
    print(f'new partition file created only with succesful downloads')
    print(f"please check {new_file_path}")
    new_dataframe.to_csv(new_file_path, index=False)
















if __name__ == '__main__':
    file_path = assert_args(parse_args())
    set_partition_file(file_path)
    
