import os
import pandas as pd
import argparse
import urllib
import gzip
import shutil
from data_loader.NCBI_data_loader import NCBIDataDownloaderPrep

ROOT_FOLDER = 'partition_file/'

# NCBI
NCBI_SEQUENCES_FOLDER = 'data/NCBI/sequence_files'
NCBI_PREPPED_SEQUENCES_FOLDER = 'data/NCBI/prep/'

# new and better
PREPPED_DATA_FOLDER = 'data/prepped/'


'''
NOTE:

when ncbi
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
        default = 'partition_file/phylum_Myxococcota_tiny.csv',
        help = 'Choose taxonomy level'
    )
    parser.add_argument(
        '-n',
        '--ncbi',
        action = 'store_true',
        help = 'Choose taxonomy level'
    )
    return parser.parse_args()

def assert_args(args):
    # check if ROOT_FOLDER exists
    if not os.path.isdir(ROOT_FOLDER):
        raise NotADirectoryError(f'{ROOT_FOLDER} is not found! please use create a partition first')

    #file_path = ROOT_FOLDER+args.file+'.csv'
    file_path = args.file
    if not os.path.isfile(file_path):
        raise NotADirectoryError(f'{file_path} is not found! create a partition!')
    
    if not os.path.isdir(PREPPED_DATA_FOLDER):
        os.makedirs(PREPPED_DATA_FOLDER)
        print("creating annotation folder : ", PREPPED_DATA_FOLDER)
        
    if args.ncbi:
        print('Downloading and prepping NCBI files!')
    else:
        pass
        


    return args.ncbi, file_path

def download_file(asm_name, link):
    # Taken from NCBI_DATA_LOADER
    if not os.path.isfile(NCBI_SEQUENCES_FOLDER+f'/{asm_name}.fna'):
        #download link
        print(f'downloading {asm_name}, at {link} ')
        urllib.request.urlretrieve(link, NCBI_SEQUENCES_FOLDER+f'/{asm_name}.fna.gz')
        with gzip.open(NCBI_SEQUENCES_FOLDER+f'/{asm_name}.fna.gz', 'rb') as f_in:
            with open(NCBI_SEQUENCES_FOLDER+f'/{asm_name}.fna', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(NCBI_SEQUENCES_FOLDER+f'/{asm_name}.fna.gz')
    else:
        print(f"File: {NCBI_SEQUENCES_FOLDER+'/'+asm_name+'.fna'} exists")
    
    
def create_csv(ncbi_status, asm_name, link):
    if ncbi_status:
        ncbi_object = NCBIDataDownloaderPrep(asm_name,link)
        ncbi_object.to_HGTDB('csv', PREPPED_DATA_FOLDER+asm_name+'.csv')
    else:
        non_ncbi_object = NCBIDataDownloaderPrep(asm_name,link)
        non_ncbi_object.to_HGTDB('csv', asm_name.replace('.fasta','.csv'))
        
        
def set_partition_file(ncbi_status, file_path):
    dataframe = pd.read_csv(file_path)
    print('Downloading the files ')
    
    if ncbi_status:
        for i in range(len(dataframe)):
            try:
                download_file(ncbi_status,dataframe.loc[i,'asm_name'], dataframe.loc[i,'dl_link'])
            except:
                print(f"Unable to download {dataframe.loc[i,'asm_name']}, {dataframe.loc[i,'dl_link']} ")
        print('Files are downloaded')
        
        list_of_unavailable_csv_files = []
        print('Creating csvs for NCBI files')
        for i in range(len(dataframe)):
            try:
                create_csv(ncbi_status,dataframe.loc[i,'asm_name'], dataframe.loc[i,'dl_link'])
            except:
                list_of_unavailable_csv_files.append(i)
                print(f"Unable to create csv for {dataframe.loc[i,'asm_name']} ")
    else:
        list_of_unavailable_csv_files = []
        print('Creating csvs for NON-NCBI files')
        for i in range(len(dataframe)):
            create_csv(ncbi_status, dataframe.loc[i,'file_path'], None)
    print('#'*80)
    new_file_path = file_path.replace('.csv', '_RTR.csv')
    new_dataframe = dataframe.drop(index=list_of_unavailable_csv_files)
    if not ncbi_status:
        new_dataframe['file_path'] = [i.replace('.fasta', '.csv') for i in new_dataframe['file_path']]
    print(f'new partition file created only with succesful downloads')
    print(f"please check {new_file_path}")
    new_dataframe.to_csv(new_file_path, index=False)
















if __name__ == '__main__':
    ncbi_status, file_path = assert_args(parse_args())
    set_partition_file(ncbi_status,file_path)
    
