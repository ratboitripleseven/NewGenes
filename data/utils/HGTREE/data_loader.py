import os
import sys
sys.path.append("../../../../NewGenes")
from Bio import Entrez
import urllib
import gzip
import shutil

from data.utils.NCBI.utils.prep_genome import prep_genome
from data.utils.NCBI.data_loader import NCBIDataLoader

SEQUENCES_FOLDER = 'data/HGTREE/sequence_files'
CODE = {
    'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
    'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C',
    'tta': 'L', 'tca': 'S', 'taa': '*', 'tga': '*',
    'ttg': 'L', 'tcg': 'S', 'tag': '*', 'tgg': 'W',
    'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R',
    'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
    'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R',
    'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R',
    'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S',
    'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S',
    'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R',
    'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',
    'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G',
    'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G',
    'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
    'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G'
}

class HGTREEDataLoader(NCBIDataLoader):
    def __init__(self, name):
        super().__init__(name)
        
    def _prep_genome(self, name)->dict:
        # print('Prepping')
        return prep_genome(SEQUENCES_FOLDER,name)
    
    def _download_genome(self, name):
        #def get_assemblies(term, download=True, path='sequence_files'):
        """Download genbank assemblies for a given search term.
        https://dmnfarrell.github.io/bioinformatics/assemblies-genbank-python
        Args:
            term: search term, usually organism name
            download: whether to download the results
            path: folder to save to
        """

        #provide your own mail here
        Entrez.email = "adjie.salman@stud.uni-due.de"
        Entrez.api_key = "726c2709d1827c981a38403d4a7d99e5cf08"
        handle = Entrez.esearch(db="assembly", term=name, retmax='200')
        record = Entrez.read(handle)
        ids = record['IdList']
        print (f'found {len(ids)} ids')
        links = []
        for id in ids:
            #get summary
            summary = self._get_assembly_summary(id)
            #get ftp link
            # so FtpPath_RefSeq gives you refseq
            # FtpPath_GenBank gives you genbank
            url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
            
            if url == '':
                continue
            label = os.path.basename(url)
            #get the fasta link - change this to get other formats
            link = os.path.join(url,label+'_cds_from_genomic.fna.gz')
            print (link)
            links.append(link)
            if not os.path.isfile(SEQUENCES_FOLDER+f'/{name}.fna'):
                #download link
                urllib.request.urlretrieve(link, SEQUENCES_FOLDER+f'/{name}.fna.gz')
                with gzip.open(SEQUENCES_FOLDER+f'/{name}.fna.gz', 'rb') as f_in:
                    with open(SEQUENCES_FOLDER+f'/{name}.fna', 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                os.remove(SEQUENCES_FOLDER+f'/{name}.fna.gz')
        