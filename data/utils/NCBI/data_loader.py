'''
This Script/module is supposed to recreate what the paper
'Horizontal Gene Transfer in Bacterial and Archaeal Complete Genomes'
calculated

TODO: currently cub is sum of all cubs in genes (is this right?)
'''
import os
import sys
sys.path.append("../../../../NewGenes")
import argparse
import pandas as pd
from data.utils.NCBI.utils.prep_genome import prep_genome
from data.utils.NCBI.metrics.calc_gc_content import calc_gc_content
from data.utils.NCBI.metrics.calc_relative_freq import calc_relative_freq
from data.utils.NCBI.metrics.calc_12_symbols import calc_12_symbols
from data.utils.NCBI.metrics.calc_48_symbols import calc_48_symbols
from data.utils.NCBI.metrics.calc_cub import calc_cub
import unittest
import math 
from Bio.SeqUtils import GC123
from Bio import Entrez
import urllib
import gzip
import shutil
SEQUENCES_FOLDER = 'data/NCBI/sequence_files'
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

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--f',
        type = str,
        default = 'bsub',
        help = 'ta'
    )


class NCBIDataLoader():
    def __init__(self, name):
        
        self.name = name
        
        # keys are locust_tag
        # not really a good way to handle error!
        self.genes = self._prep_genome(name=name)
        if self.genes == 0:
            self._download_genome(name=name)
            self.genes = self._prep_genome(name=name)
        
        
        #self.genes = pd.DataFrame.from_dict(self._prep_genome(name=name), orient='index')
        
        #init metrics to 0
        self.complete_sequence = ''
        
        # hgtdb
        self.mean_GCT = 0
        self.mean_GC1 = 0
        self.mean_GC2 = 0
        self.mean_GC3 = 0
        self.std_GCT = 0
        self.std_GC1 = 0
        self.std_GC2 = 0
        self.std_GC3 = 0
        
        self.mean_cub = {}
        self.std_cub = {}
        
        self.cub = {}
        for cds in CODE:
            self.std_cub[cds.upper()] = 0
            self.cub[cds.upper()] = 0
        
        # https://doi.org/10.1093/nar/gkm204
        self.rel_freq = None
        self.nucleutide_identity = None
        self.dinucleutide_identity = None
        
        
        
        
        ####
        self.rscu = 0
        self.aa_comp = 0
        
        #fill metrics
        self._fill_metrics_genome()
        #self._fill_metrics_genes()
        
    def _prep_genome(self, name)->dict:
        # print('Prepping')
        return prep_genome(SEQUENCES_FOLDER,name)
    
    def _get_assembly_summary(self, id):
        """Get esummary for an entrez id"""
        esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
        esummary_record = Entrez.read(esummary_handle)
        return esummary_record
    
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
            url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
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
        # now no need to return!
        # return links
        
        
    def _fill_metrics_genome(self):
        # fill metrics based on https://doi.org/10.1093/nar/gkm204 and HGTdb
        for tag in self.genes:
            # hgtdb
            # acess gene
            gene = self.__getitem__(key=tag)
            gene['g_count'], gene['a_count'], gene['c_count'], gene['t_count'], _ = calc_gc_content(gene['sequence'])
            gene['GCT'], gene['GC1'], gene['GC2'], gene['GC3'] = GC123(gene['sequence'])
            
            gene['cub'] = calc_cub(gene['sequence'])
            
            # add to mean calcilattion
            self.mean_GCT += gene['GCT']
            self.mean_GC1 += gene['GC1']
            self.mean_GC2 += gene['GC2']
            self.mean_GC3 += gene['GC3']
            
            for cds in CODE:
                self.cub[cds.upper()] += gene['cub'][cds.upper()]
            # self.mean_cub += gene['cub']
            
            # get complete sequence
            self.complete_sequence += gene['sequence']
            
            # https://doi.org/10.1093/nar/gkm204
            gene['rel_freq'] = calc_relative_freq(gene['sequence'])
            gene['12_symbols'] = calc_12_symbols(gene['sequence'])
            gene['48_symbols'] = calc_48_symbols(gene['sequence'])
            
        
        # https://doi.org/10.1093/nar/gkm204
        self.rel_freq = calc_relative_freq(self.complete_sequence)
        self.nucleutide_identity = calc_12_symbols(self.complete_sequence)
        self.dinucleutide_identity = calc_48_symbols(self.complete_sequence)
        # self.cub = calc_cub(self.complete_sequence)
        
        #finalize mean
        self.mean_GCT = self.mean_GCT/len(self.genes)
        self.mean_GC1 = self.mean_GC1/len(self.genes)
        self.mean_GC2 = self.mean_GC2/len(self.genes)
        self.mean_GC3 = self.mean_GC3/len(self.genes)
        
        # mean for cub (TODO:recheck this!!!)
        for cds in CODE:
            self.mean_cub[cds.upper()]=self.cub[cds.upper()]/len(self.genes)
        # self.mean_cub = self.mean_cub/len(self.genes)
        
        # calculate standard devation of whole genome
        numeratorT = 0
        numerator1 = 0
        numerator2 = 0
        numerator3 = 0
        numeratorCub = {}
        for cds in CODE:
            numeratorCub[cds.upper()] = 0
        for tag in self.genes:
            gene = self.__getitem__(key=tag)
            deltaT = (gene['GCT'] - self.mean_GCT)
            delta1 = (gene['GC1'] - self.mean_GC1)
            delta2 = (gene['GC2'] - self.mean_GC2)
            delta3 = (gene['GC3'] - self.mean_GC3)
            
            deltaCub = {}
            for cds in CODE:
                deltaCub[cds.upper()] = (gene['cub'][cds.upper()] - self.mean_cub[cds.upper()])
            # deltaCub = (gene['cub'] - self.mean_cub)
            
            deltaT *= deltaT
            delta1 *= delta1
            delta2 *= delta2
            delta3 *= delta3
            for cds in CODE:
                squared = deltaCub[cds.upper()]**2
                deltaCub[cds.upper()] = squared
            # deltaCub *= deltaCub
            numeratorT += deltaT
            numerator1 += delta1
            numerator2 += delta2
            numerator3 += delta3
            for cds in CODE:
                numeratorCub[cds.upper()] += deltaCub[cds.upper()]
            
            # numeratorCub += numeratorCub
        
        self.std_GCT = math.sqrt(numeratorT/len(self.genes))
        self.std_GC1 = math.sqrt(numerator1/len(self.genes))
        self.std_GC2 = math.sqrt(numerator2/len(self.genes))
        self.std_GC3 = math.sqrt(numerator3/len(self.genes))
        
        for cds in CODE:
            self.std_cub[cds.upper()] = math.sqrt(numeratorCub[cds.upper()]/len(self.genes))
        # self.std_cub = math.sqrt(numeratorCub/len(self.genes))
        
        # add std calculation for each gene
        for tag in self.genes:
            gene = self.__getitem__(key=tag)
            gene['SDT'] = (gene['GCT'] - self.mean_GCT)/self.std_GCT
            gene['SD1'] = (gene['GC1'] - self.mean_GC1)/self.std_GC1
            gene['SD2'] = (gene['GC2'] - self.mean_GC2)/self.std_GC2
            gene['SD3'] = (gene['GC3'] - self.mean_GC3)/self.std_GC3
            
            for cds in CODE:
                gene['std_cub'][cds.upper()] = (gene['cub'][cds.upper()] - self.mean_cub[cds.upper()])/self.std_cub[cds.upper()]
            
        
        
            
            
            
            
        
        
        
        
    def _fill_metrics_genes(self):
        #second loop to calc SD1,2,3,T for each gene
        print('here')
        
            
        
            
        
        
        
        
        
        
    def __len__(self):
        return len(self.genes)
    
    def __getitem__(self, key):
        return self.genes[key]
    

        
    def list_genes(self):
        for genes in self.genes:
            print(genes)
            
    def print_genome_summary(self):
        print(f'Mean GC Content-> T:{self.mean_GCT}, 1:{self.mean_GC1}, 2:{self.mean_GC2}, 3:{self.mean_GC3}')
        print(f'Std GC content-> T:{self.std_GCT}, 1:{self.std_GC1}, 2:{self.std_GC2}, 3:{self.std_GC3}')
        print(f'Relative nucleotide frequency: {self.rel_freq}')
        print(f'Nucleotide Identity: {self.nucleutide_identity}')
        print(f'Dinucleotide Identity: {self.dinucleutide_identity}')
        print(f'Codon Usage Bias: {self.cub}')
        
    def print_gene_summary(self, tag):
        gene = self.__getitem__(key=tag)
        print(f'Mean GC Content-> T: {gene["GCT"]}, 1:{gene["GC1"]}, 2:{gene["GC2"]}, 3:{gene["GC3"]}')
        print(f'Std GC content-> T:{gene["SDT"]}, 1:{gene["SD1"]}, 2:{gene["SD2"]}, 3:{gene["SD3"]}')
        print(f'Relative nucleotide frequency: {gene["rel_freq"]}')
        print(f'Nucleotide Identity: {gene["12_symbols"]}')
        print(f'Dinucleotide Identity: {gene["48_symbols"]}')
        print(f'Codon Usage Bias: {gene["cub"]}')
        
        
        
        
        

## for testing

class TestNCBIDataLoader(unittest.TestCase):
    
    
    def test_init_positive(self):
        genome = NCBIDataLoader('ecoli')
        assert genome.genes is not None, "Somthing is wrong when reading file"
    
    def test_prep_genome(self):
        genome = NCBIDataLoader('ecoli')
        assert len(genome) != 0, "cannot access length"
        
    




if __name__ == '__main__':
    #args = parse_args()
    unittest.main()