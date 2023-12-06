'''
This Script/module is supposed to recreate what the paper
'Horizontal Gene Transfer in Bacterial and Archaeal Complete Genomes'
calculated

TODO: currently cub is sum of all cubs in genes (is this right?)
'''
import os
import sys
sys.path.append("../NewGenes")
import argparse
import pandas as pd
from data_loader.utils.prep_genome import prep_genome
from data_loader.utils.annotate_HGTs import *
from data_loader.metrics.calc_gc_content import calc_gc_content
from data_loader.metrics.calc_relative_freq import calc_relative_freq
from data_loader.metrics.calc_12_symbols import calc_12_symbols
from data_loader.metrics.calc_48_symbols import calc_48_symbols
from data_loader.metrics.calc_cub import calc_cub
from data_loader.metrics.calc_RSCU_and_RFC import calc_RSCU_and_RFC
import unittest
import math 
from Bio.SeqUtils import GC123
from Bio import Entrez
import urllib
import gzip
import shutil
SEQUENCES_FOLDER = 'data/NCBI/sequence_files'
PREPPED_SEQUENCES_FOLDER = 'data/NCBI/prep/'
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
    
class NCBIDataLoader:
    def __init__(self, partition_file, data_type = 'F'):
        self.partition_file = partition_file
        self.partition_frame = None
        self._check_prepped_file_availability()
        self.data_type = data_type
        
        
        
        
    
    def _check_prepped_file_availability(self):
        self.partition_frame  = pd.read_csv(self.partition_file)
        for i in range(len(self.partition_frame)):
            identifier = self.partition_frame.loc[i,'GenBank accession number']
            # check if file is prepped
            if not os.path.isfile(PREPPED_SEQUENCES_FOLDER+identifier+'.csv'):
                temp_downloader = NCBIDataDownloaderPrep(identifier)
                temp_downloader.to_HGTDB('csv')
    
    def dataset_prep(self):
        if self.data_type == 'A':
            raise NotImplementedError('Not yet implemented')
        elif self.data_type == 'B':
            raise NotImplementedError('Not yet implemented')
        elif self.data_type == 'C':
            raise NotImplementedError('Not yet implemented')
        elif self.data_type == 'D':
            raise NotImplementedError('Not yet implemented')
        elif self.data_type == 'E':
            raise NotImplementedError('Not yet implemented')
        elif self.data_type == 'F':
            columns_to_drop = ["gene", "protein_id","protein", "location", "GC1","GC2","GC3", "GCT", "Sim1", "Sim2", "Sim3", "SimT", "SimGC", "SimMah"]
            columns = 6
        else:
            raise ValueError(f'No data type {self.data_type}')
        
        X_train = np.array([]).reshape(0,columns)
        y_train = np.array([]).reshape(0,)
        X_test = np.array([]).reshape(0,columns)
        y_test = np.array([]).reshape(0,)
        
        # TODO: Add sleep here!
        for i in range(len(self.partition_frame)):
            identifier = self.partition_frame.loc[i,'GenBank accession number']
            temp_data = pd.read_csv(PREPPED_SEQUENCES_FOLDER+identifier+'.csv', index_col=0)
            temp_data['HGT'] = temp_data['HGT'].fillna(0)
            temp_data['HGT'] = temp_data['HGT'].replace('H',1)
            # TODO:think about this better
            temp_data[ temp_data['Dev.AA'] != '[]'] = 1
            temp_data['Dev.AA'] = temp_data['Dev.AA'].replace('[]', 0)
            
            temp_data = temp_data.drop(columns=columns_to_drop)
            array = temp_data.values
            
            X,y = array[:,0:-1], array[:,-1]
            if self.partition_frame.loc[i,'partition'] == 'train':
                X_train = np.concatenate([X, X_train], axis = 0)
                y_train = np.concatenate([y, y_train], axis = 0)
            else:
                X_test = np.concatenate([X, X_test], axis = 0)
                y_test = np.concatenate([y, y_test], axis = 0)
        
        return X_train, y_train, X_test, y_test
            
        


class NCBIDataDownloaderPrep:
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
        
        self.mean_RSCU = {}
        self.std_RSCU = {}
        
        self.mean_RFC = {}
        self.std_RFC = {} 
        
        self.cub = {}
        for cds in CODE:
            self.std_cub[cds.upper()] = 0
            self.cub[cds.upper()] = 0
            
            self.mean_RSCU[cds.upper()] = 0
            self.std_RSCU[cds.upper()] = 0
            
            self.mean_RFC[cds.upper()] = 0
            self.std_RFC[cds.upper()] = 0
            
        # https://doi.org/10.1093/nar/gkm204
        self.rel_freq = None
        self.nucleutide_identity = None
        self.dinucleutide_identity = None
        
        
        
        
        ####
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
            # for now going forward no refseq files are going to be used!
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
            
            # add codon usage bias per gene
            gene['cub'] = calc_cub(gene['sequence'])
            
            # add RSCU and Relative Frequency of Codon Usage
            gene['RSCU'], gene['RFC'] = calc_RSCU_and_RFC(gene['cub'])
            
            
            # add to mean calcilattion
            self.mean_GCT += gene['GCT']
            self.mean_GC1 += gene['GC1']
            self.mean_GC2 += gene['GC2']
            self.mean_GC3 += gene['GC3']
            
            # calculate cub
            # the mean RSCU and RFC is only intermediate!
            for cds in CODE:
                self.cub[cds.upper()] += gene['cub'][cds.upper()]
                self.mean_RSCU[cds.upper()] += gene['RSCU'][cds.upper()]
                self.mean_RFC[cds.upper()] += gene['RFC'][cds.upper()]
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
        
        # mean for cub, RSCU and RFC
        for cds in CODE:
            self.mean_cub[cds.upper()]=self.cub[cds.upper()]/len(self.genes)
            self.mean_RSCU[cds.upper()]=self.mean_RSCU[cds.upper()]/len(self.genes)
            self.mean_RFC[cds.upper()]=self.mean_RFC[cds.upper()]/len(self.genes)
        # self.mean_cub = self.mean_cub/len(self.genes)
        
        # https://math.stackexchange.com/questions/1433374/difference-between-these-two-standard-deviation-formulas
        # calculate standard devation of whole genome (population)
        numeratorT = 0
        numerator1 = 0
        numerator2 = 0
        numerator3 = 0
        numeratorCub = {}
        numeratorRSCU = {}
        numeratorRFC = {}
        for cds in CODE:
            numeratorCub[cds.upper()] = 0
            numeratorRSCU[cds.upper()] = 0
            numeratorRFC[cds.upper()] = 0
        for tag in self.genes:
            gene = self.__getitem__(key=tag)
            deltaT = (gene['GCT'] - self.mean_GCT)
            delta1 = (gene['GC1'] - self.mean_GC1)
            delta2 = (gene['GC2'] - self.mean_GC2)
            delta3 = (gene['GC3'] - self.mean_GC3)
            
            deltaCub = {}
            deltaRSCU = {}
            deltaRFC = {}
            for cds in CODE:
                deltaCub[cds.upper()] = (gene['cub'][cds.upper()] - self.mean_cub[cds.upper()])
                deltaRSCU[cds.upper()] = (gene['RSCU'][cds.upper()] - self.mean_RSCU[cds.upper()])
                deltaRFC[cds.upper()] = (gene['RFC'][cds.upper()] - self.mean_RFC[cds.upper()])
            # deltaCub = (gene['cub'] - self.mean_cub)
            
            deltaT *= deltaT
            delta1 *= delta1
            delta2 *= delta2
            delta3 *= delta3
            for cds in CODE:
                squared_cub = deltaCub[cds.upper()]**2
                squared_RSCU = deltaRSCU[cds.upper()]**2
                squared_RFC = deltaRFC[cds.upper()]**2
                
                deltaCub[cds.upper()] = squared_cub
                deltaRSCU[cds.upper()] = squared_RSCU
                deltaRFC[cds.upper()] = squared_RFC
            # deltaCub *= deltaCub
            numeratorT += deltaT
            numerator1 += delta1
            numerator2 += delta2
            numerator3 += delta3
            for cds in CODE:
                numeratorCub[cds.upper()] += deltaCub[cds.upper()]
                numeratorRSCU[cds.upper()] += deltaRSCU[cds.upper()]
                numeratorRFC[cds.upper()] += deltaRFC[cds.upper()]
            
            # numeratorCub += numeratorCub
        
        self.std_GCT = math.sqrt(numeratorT/len(self.genes))
        self.std_GC1 = math.sqrt(numerator1/len(self.genes))
        self.std_GC2 = math.sqrt(numerator2/len(self.genes))
        self.std_GC3 = math.sqrt(numerator3/len(self.genes))
        
        for cds in CODE:
            self.std_cub[cds.upper()] = math.sqrt(numeratorCub[cds.upper()]/len(self.genes))
            self.std_RSCU[cds.upper()] = math.sqrt(numeratorRSCU[cds.upper()]/len(self.genes))
            self.std_RFC[cds.upper()] = math.sqrt(numeratorRFC[cds.upper()]/len(self.genes))
        # self.std_cub = math.sqrt(numeratorCub/len(self.genes))
        
        # add std calculation for each gene (sample)
        for tag in self.genes:
            gene = self.__getitem__(key=tag)
            gene['SDT'] = (gene['GCT'] - self.mean_GCT)/self.std_GCT
            if abs(gene['SDT']) >= 1.5:
                if abs(gene['SDT']) >= 2:
                    if gene['SDT'] > 0:
                        gene['SimT'] = +2
                    else:
                        gene['SimT'] = -2
                else:
                    if gene['SDT'] > 0:
                        gene['SimT'] = +1
                    else:
                        gene['SimT'] = -1
            else:
                gene['SimT'] = 0
            
            gene['SD1'] = (gene['GC1'] - self.mean_GC1)/self.std_GC1
            if abs(gene['SD1']) >= 1.5:
                if abs(gene['SD1']) >= 2:
                    if gene['SD1'] > 0:
                        gene['Sim1'] = +2
                    else:
                        gene['Sim1'] = -2
                else:
                    if gene['SD1'] > 0:
                        gene['Sim1'] = +1
                    else:
                        gene['Sim1'] = -1
            else:
                gene['Sim1'] = 0
                
            gene['SD2'] = (gene['GC2'] - self.mean_GC2)/self.std_GC2
            if abs(gene['SD2']) >= 1.5:
                if abs(gene['SD2']) >= 2:
                    if gene['SD2'] > 0:
                        gene['Sim2'] = +2
                    else:
                        gene['Sim2'] = -2
                else:
                    if gene['SD2'] > 0:
                        gene['Sim2'] = +1
                    else:
                        gene['Sim2'] = -1
            else:
                gene['Sim2'] = 0
            
            gene['SD3'] = (gene['GC3'] - self.mean_GC3)/self.std_GC3
            if abs(gene['SD3']) >= 1.5:
                if abs(gene['SD3']) >= 2:
                    if gene['SD3'] > 0:
                        gene['Sim3'] = +2
                    else:
                        gene['Sim3'] = -2
                else:
                    if gene['SD3'] > 0:
                        gene['Sim3'] = +1
                    else:
                        gene['Sim3'] = -1
            else:
                gene['Sim3'] = 0
                
            # simGC
            count = 0
            # if simt is more than 1.5
            if gene['SimT'] >=1:
                count+=1
                if gene['SimT'] >=2:
                    count+=1
                # if sim1 and sim3 is equal sign and if 
            elif (gene['SD1'] * gene['SD3']) > 0:
                if gene['Sim1'] >=1:
                    count+=1
                elif gene['Sim3'] >=1:
                    count+=1
            gene['SimGC']=count
            
            for cds in CODE:
                gene['std_cub'][cds.upper()] = (gene['cub'][cds.upper()] - self.mean_cub[cds.upper()])/self.std_cub[cds.upper()]
                
        
            
        
        
            
            
            
            
        
        
        
        
    def _fill_metrics_genes(self):
        #second loop to calc SD1,2,3,T for each gene
        print('here')
        
            
        
            
        
        
        
        
        
        
    def __len__(self):
        return len(self.genes)
    
    def __getitem__(self, key):
        return self.genes[key]
    

    def to_HGTDB(self, return_type='pd'):
        '''
        Return the whole object into a dataframe style like HGTDB
        or print a csv
        
        '''
        if return_type not in ['pd', 'csv']:
            raise ValueError('argument unknown! choose pd or csv')
        
            
        
        # get HGT candidates based on GC content
        hgt_candidates_GC = GC_Content_Deviation(self)
        
        # get list of to exclude in list above based on Amino Acid content
        calculate_amino_acid_content_genome_mean(self)
        calculate_amino_acid_content_gene_mean(self)
        calculate_amino_acid_content_genome_std(self)
        list_of_non_extraneous_genes_AA = check_amino_acid_deviation(self)
        
        # coombine information from above
        extraneous_but_non_HGT = []
        deemed_HGT =[]
        for i in hgt_candidates_GC:
            if i in list_of_non_extraneous_genes_AA:
                extraneous_but_non_HGT.append(i)
            else:
                deemed_HGT.append(i)
        
        
        # calculate mahalanobis distances
        calculate_mahalanobis_distances(self)
        hgt_candidates_Mah = get_potential_HGT_Mah(self)

        
        print(len(extraneous_but_non_HGT))
        print(len(deemed_HGT))
        combined_list_hgt = deemed_HGT+hgt_candidates_Mah
        
        print('Total hgt by gc and mah: {}'.format(len(combined_list_hgt)))
        for tag in self.genes:
            gene = self.__getitem__(key=tag)
            if tag in combined_list_hgt:
                gene['HGT']='H'
            else:
                gene['HGT']=None
        
        # output to as pandas
        if return_type == 'pd':
            # return a pandas dataframe object
            return pd.DataFrame.from_dict(self.genes, orient='index')
        elif return_type == 'csv':
            # return 0 
            # print a hgdtb like csv file
            if not os.path.isdir(PREPPED_SEQUENCES_FOLDER):
                os.makedirs(PREPPED_SEQUENCES_FOLDER)
                print("creating folder : ", PREPPED_SEQUENCES_FOLDER)
            df = pd.DataFrame.from_dict(self.genes, orient='index')
            print('printing csv... saving as {}.csv'.format(self.name))
            df.drop(columns=['g_count','a_count','c_count','t_count', 'rel_freq', '12_symbols', '48_symbols', 'cub','sequence','std_cub', 'RSCU', 'RFC', 'AA_Content_mean'], inplace=True)
            df['Mah'] = [float(x) for x in df['Mah']]
            df.to_csv(PREPPED_SEQUENCES_FOLDER+'{}.csv'.format(self.name))
            print('done printing')
            return 0
        
        
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

class TestNCBIDataLoaderPrep(unittest.TestCase):
    
    def test_dataloader(self):
        dataloader = NCBIDataLoader('partition_file/class_chlorobia.csv')
        x1,y1,x2,y2 = dataloader.dataset_prep()
        assert len(x1)!=0, "error!"
    
    
    #def test_init_positive(self):
    #    genome = NCBIDataDownloaderPrep('AE000657')
    #    assert genome.genes is not None, "Somthing is wrong when reading file"
    
    #def test_prep_genome(self):
    #    genome = NCBIDataDownloaderPrep('AL009126')
    #    assert len(genome) != 0, "cannot access length"
        
    #def test_prep_genome_csv_out(self):
    #    genome = NCBIDataDownloaderPrep('AL009126')
    #    test = genome.to_HGTDB('csv')
    #    assert test == 0, 'something went wrong in creating csv'
        
        
    




if __name__ == '__main__':
    #args = parse_args()
    unittest.main()