import os
import pandas as pd
import argparse
import urllib
import gzip
import shutil
from data_loader.NCBI_data_loader import NCBIDataDownloaderPrep
import unittest

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
        default = 'test/CP002542.fna',
        help = 'Input file fna or fasta'
    )
    return parser.parse_args()

def assert_args(args):
    file_path = args.file
    
    # check if file ends with .fasta or .fna
    if file_path.split('.')[-1] not in ['fasta','fna']:
        raise ValueError(f'{file_path} is not a fasta or fna file!')
    
    if not os.path.isfile(file_path):
        raise NotADirectoryError(f'{file_path} Does not exist!')

    return file_path
    
    
def create_csv(file_path):
    non_ncbi_object = NCBIDataDownloaderPrep(file_path,None)
    extension = file_path.split('.')[-1]
    non_ncbi_object.to_HGTDB('csv', file_path.replace(extension,'csv'))
        
        


class TestPrepFile(unittest.TestCase):
    def test_create_csv(self):
        create_csv('test/CP002542.fna')
        
        assert os.path.isfile('test/CP002542.csv'), "error target file not created!"
        







if __name__ == '__main__':
    file_path = assert_args(parse_args())
    create_csv(file_path)
    
