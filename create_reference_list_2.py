from Bio import Entrez
import pandas as pd
import urllib
import gzip
import shutil
import os

REFERENCE_ROOT = 'reference/'

def download_files():
    if not os.path.isdir(REFERENCE_ROOT):
        os.makedirs(REFERENCE_ROOT)
        print("creating folder : ", REFERENCE_ROOT)
        
    if not os.path.isfile(REFERENCE_ROOT+'assembly_summary_genbank.txt'):
        urllib.request.urlretrieve('ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt', REFERENCE_ROOT+'assembly_summary_genbank.txt')
        
    if not os.path.isfile(REFERENCE_ROOT+'bac120_metadata.tsv.gz'):
        urllib.request.urlretrieve('https://data.gtdb.ecogenomic.org/releases/latest/bac120_metadata.tsv.gz', REFERENCE_ROOT+'bac120_metadata.tsv.gz')
    
    if not os.path.isfile(REFERENCE_ROOT+'bac120_metadata.tsv'):
        with gzip.open(REFERENCE_ROOT+'bac120_metadata.tsv.gz', 'rb') as f_in:
            with open(REFERENCE_ROOT+'bac120_metadata.tsv', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)



def create_partition():
    gtdb = pd.read_csv(REFERENCE_ROOT+'bac120_metadata.tsv', sep='\t', usecols=['gtdb_taxonomy', 'ncbi_genbank_assembly_accession'])
    assembly_ref = pd.read_csv(REFERENCE_ROOT+'assembly_summary_genbank.txt', sep='\t', skiprows=1, usecols=['#assembly_accession', 'ftp_path','organism_name', 'asm_name'])
    #assembly_ref = assembly_ref[assembly_ref['asm_name'].str.contains('ASM', na=False)].copy()
    gtdb_and_links = pd.merge(gtdb,assembly_ref,how='left', left_on='ncbi_genbank_assembly_accession', right_on='#assembly_accession')
    print(f'the length of gtdb {len(gtdb)}')
    print(f'the length of assembly_ref {len(assembly_ref)}')
    print(f'the length of gtdb_and_links {len(gtdb_and_links)}')
    
    #rename
    # TODO: need to change everything form gdtb to gtdb
    gtdb_and_links.rename(columns = {'gtdb_taxonomy':'GDTB_taxonomy', 'ftp_path': 'link'}, inplace = True)
    # delete anything not ASM 
    gtdb_and_links = gtdb_and_links[gtdb_and_links['asm_name'].str.contains('ASM', na=False)].copy()
    gtdb_and_links.to_csv(REFERENCE_ROOT+'ref_gtdb.csv')
    print(gtdb_and_links)
    
    
        

if __name__ == '__main__':
    download_files()
    create_partition()