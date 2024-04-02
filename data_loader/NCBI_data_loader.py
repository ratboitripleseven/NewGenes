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
import torch
from sklearn import preprocessing
from torch.nn.utils.rnn import pad_sequence
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
    '''
    This is for classical (?)
    '''
    def __init__(self, partition_file, data_type = 'F'):
        self.partition_file = partition_file
        self.partition_frame = None
        self._check_prepped_file_availability()
        self.data_type = data_type
        
    
    def _check_prepped_file_availability(self):
        self.partition_frame  = pd.read_csv(self.partition_file)
        for i in range(len(self.partition_frame)):
            identifier = self.partition_frame.loc[i,'asm_name']
            link = self.partition_frame.loc[i,'link']
            # check if file is prepped
            if not os.path.isfile(PREPPED_SEQUENCES_FOLDER+identifier+'.csv'):
                print(f'Getting {identifier}')
                temp_downloader = NCBIDataDownloaderPrep(identifier, link)
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
        
        for i in range(len(self.partition_frame)):
            identifier = self.partition_frame.loc[i,'asm_name']
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
            
        
class NCBIDatasetSequential(torch.utils.data.Dataset):
    def __init__(self, data_type, partition_file, partition):
        '''
        Note: Annotate is not suppose to have any labls as it is supposed to be used to annotate them!
        
        '''
        # Annotate bypasses any partition and just us th whole partition file to load
        if partition not in ['train','test','valid','annotate']:
            raise ValueError('not partition or train, test or valid!')
        
        
        
        self.partition = partition
        self.data_type = data_type
        self.partition_file = partition_file
        self.min_max_scaler = preprocessing.MinMaxScaler()
        # self.one_hot_encoder = preprocessing.OneHotEncoder()
        # self.drop_na = drop_na
        self.null_count = 0
        self.na_count = 0
        
        self.max_sequence = 0
        
        # load all in ram perhaps?
        self.data_x = []
        self.data_y = []
        self.data_seq_length = []
        self.data_id = []
        
        
        self._init_dataset()
        
    def _load_single_file(self, species, data_type):
        '''
        I took this from HGTDB squential data laoder
        '''

        preprocessed_path = PREPPED_SEQUENCES_FOLDER
        
        csv_list = os.listdir(preprocessed_path)
        
    
        if species is not None:
            csv_file = None
            for path in csv_list:
                # check if current path is a file
                if os.path.basename(path).replace(".csv", "") == species:
                    csv_file =os.path.join(preprocessed_path, path)
                    df = pd.read_csv(csv_file)
            if csv_file is None:
                raise ValueError(f'{species} not found')
        else:
            raise ValueError(f'{species} not found')
        
        # need to drop the first unnamed column!
        df = df.drop(df.columns[0],axis=1)
            
        # set dataset according to data type
        if data_type == 'A':
            df = df.drop(columns=["gene","protein","protein_id","location","SD1","SD2","SD3","SDT","Mah","Sim1","Sim2","Sim3","SimT","SimMah","Dev.AA","SimGC", "HGT"]) # only GC1,GC2,GC3,GCT
        elif data_type == 'B':
            df = df.drop(columns=["FunctionCode","Strand","Length","SD1","SD2","SD3","SDT"]) # only GC1,GC2,GC3,GCT,Mah,AADev
        elif data_type == 'C':
            # C is basically A but padding is with 0.666 instead of 0
            df = df.drop(columns=["FunctionCode","Strand","AADev","Length","SD1","SD2","SD3","SDT","Mah"]) # only GC1,GC2,GC3,GCT
        elif data_type == 'D':
            # So z score is basically what is calculated for sd1,2,3,t
            df = df.drop(columns=["FunctionCode","Strand","AADev","Length","GC1","GC2","GC3","GCT","Mah"]) # only SD1,SD2,SD3,SDT
        elif data_type == 'E':
            # So z score is basically what is calculated for sd1,2,3,t
            # this z score is cutoff at two and scaled down by two so data ranges from -1 to 1
            df = df.drop(columns=["gene","protein","protein_id","location","GC1","GC2","GC3","GCT","Mah","Sim1","Sim2","Sim3","SimT","SimMah","Dev.AA","SimGC", "HGT"]) # only SD1,SD2,SD3,SDT
        elif data_type == 'F':
            # So z score is basically what is calculated for sd1,2,3,t
            # z score for mah is added here
            # this z score is cutoff at two and scaled down by two so data ranges from -1 to 1
            df = df.drop(columns=["FunctionCode","Strand","AADev","Length","GC1","GC2","GC3","GCT"]) # only SD1,SD2,SD3,SDT, Mah
        
        #print(df)
        # count nulls!
        #df.bfill(inplace=True)
        self.null_count +=df.isnull().sum().sum()
        self.na_count +=df.isna().sum().sum()
        if df.isna().sum().sum() >0:
            print(df[df.isna().any(axis=1)])
            
        df = df.bfill(axis='columns')
        #for column in df.columns:
        #    df[column] = df[column].fillna(0)
        self.null_count +=df.isnull().sum().sum()
        self.na_count +=df.isna().sum().sum()
        if df.isna().sum().sum() >0:
            print(df[df.isna().any(axis=1)])
            

        
        
        if data_type == 'D':
            pass
        elif data_type == 'E':
            df['SD1'] = [ (2*(abs(x)/x))/2 if abs(x)>2 else x/2 for x in df['SD1']]
            df['SD2'] = [ (2*(abs(x)/x))/2 if abs(x)>2 else x/2 for x in df['SD2']]
            df['SD3'] = [ (2*(abs(x)/x))/2 if abs(x)>2 else x/2 for x in df['SD3']]
            df['SDT'] = [ (2*(abs(x)/x))/2 if abs(x)>2 else x/2 for x in df['SDT']]
        elif data_type == 'F':
            # need to reindex!
            df=df.reset_index(drop=True)
            df['SD1'] = [ (2*(abs(x)/x))/2 if abs(x)>2 else x/2 for x in df['SD1']]
            df['SD2'] = [ (2*(abs(x)/x))/2 if abs(x)>2 else x/2 for x in df['SD2']]
            df['SD3'] = [ (2*(abs(x)/x))/2 if abs(x)>2 else x/2 for x in df['SD3']]
            df['SDT'] = [ (2*(abs(x)/x))/2 if abs(x)>2 else x/2 for x in df['SDT']]
            
            # calculate mean mah
            mean_Mah = df['Mah'].sum()/len(df)
            # calculate sdmah
            sum_of_mah_diffs = 0
            for gene_idx in range(len(df)):
                mah_diff = df.loc[gene_idx, 'Mah'] - mean_Mah
                squared_mah_diff = mah_diff * mah_diff
                sum_of_mah_diffs = sum_of_mah_diffs + squared_mah_diff
                
            sd_mah = math.sqrt(sum_of_mah_diffs/len(df))
            
            # calculate deviation for each val
            standard_deviation = []
            for x in df['Mah']:
                standard_deviation.append( (x - mean_Mah)/ sd_mah )
                
            # cut off and divide by 2
            standard_deviation_standardized = [ (2*(abs(x)/x))/2 if abs(x)>2 else x/2 for x in standard_deviation]
            # replace mah values (i think this is fine)
            df['Mah'] = standard_deviation_standardized
        else:
            df=(df-df.min())/(df.max()-df.min())
        array = df.values
        
        if self.partition != 'annotate':
            x = array[:,0:-1]
            y = array[:,-1]
            y = np.expand_dims(y, axis=1)
        else:
            x = array
            y = np.zeros(len(x))


        return x,y 
        
    
    def _init_dataset(self):
        '''
        '''
        partition_frame = pd.read_csv(self.partition_file)
        if self.partition != 'annotate':
            # if not annotate partition
            partition_frame = partition_frame[partition_frame['partition']==self.partition].reset_index(drop=True)
        
        
        for i in range(len(partition_frame)):
            x,y = self._load_single_file(partition_frame.loc[i,'asm_name'], self.data_type)
            
            if self.max_sequence<len(x):
                self.max_sequence = len(x)
                
            self.data_id.append(partition_frame.loc[i,'asm_name'])
            self.data_x.append(torch.from_numpy(np.float32(x)))
            self.data_y.append(torch.from_numpy(np.float32(y)))
            self.data_seq_length.append(torch.tensor(len(x)))
            
        # pad 'sequences'
        # padded stuff are tagged as 0!
        if self.data_type == 'C':
            self.data_x = pad_sequence(self.data_x, batch_first=True, padding_value=0.666)
        else:
            self.data_x = pad_sequence(self.data_x, batch_first=True)
        self.data_y = pad_sequence(self.data_y, batch_first=True)

    def __getitem__(self, ind):
        datum = self.data_x[ind]
        label = self.data_y[ind]
        seq_length = self.data_seq_length[ind]
        data_id = self.data_id[ind]
        
        output = {
            "datum" : datum,
            "seq_length" : seq_length,
            "data_id" : data_id,
            "label" : label
        }
        return output
    
    def __len__(self):
        return len(self.data_x)
        
        
        
        

class NCBIDataDownloaderPrep:
    def __init__(self, name, link=None):
        
        self.name = name
        self.link = link
        # keys are locust_tag
        # not really a good way to handle error!
        self.genes = self._prep_genome(name=name)
        if self.genes == 0:
            # self._download_genome(name=name)
            self._download_genome(self.name, self.link)
            self.genes = self._prep_genome(name=name)
        
        
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
        file_path = SEQUENCES_FOLDER+'/'+ name+'.fna'
        print(f'file path : {file_path}')
        return prep_genome(file_path)
    
    def _download_genome(self, name, link):
        #def get_assemblies(term, download=True, path='sequence_files'):
        """Download genbank assemblies for a given search term.
        https://dmnfarrell.github.io/bioinformatics/assemblies-genbank-python
        Args:
            term: search term, usually organism name
            download: whether to download the results
            path: folder to save to
        """
        
        if link:
            #download link
            print(f'downloading {name}, at {link} ')
            urllib.request.urlretrieve(link, SEQUENCES_FOLDER+f'/{name}.fna.gz')
            with gzip.open(SEQUENCES_FOLDER+f'/{name}.fna.gz', 'rb') as f_in:
                with open(SEQUENCES_FOLDER+f'/{name}.fna', 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(SEQUENCES_FOLDER+f'/{name}.fna.gz')
        else:
            raise ValueError('Link is Unspecified!')
        

        
        
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
        dataloader = NCBIDataLoader('partition_file/phylum_Fibrobacterota_test_available.csv')
        x1,y1,x2,y2 = dataloader.dataset_prep()
        assert len(x1)!=0, "error!"
    
    
    def test_NCBI_Sequential(self):
        dataset = NCBIDatasetSequential('A','partition_file/phylum_Fibrobacterota_test_available.csv','annotate' )
        print(dataset[0])
        assert len(dataset) ==2, "Something wrong wtih dataset sequential"
        
        
    
    def test_prep_genome(self):
        genome = NCBIDataDownloaderPrep('AL009126')
        assert len(genome) != 0, "cannot access length"
        
    #def test_prep_genome_csv_out(self):
    #    genome = NCBIDataDownloaderPrep('AL009126')
    #    test = genome.to_HGTDB('csv')
    #    assert test == 0, 'something went wrong in creating csv'
        
        
    




if __name__ == '__main__':
    #args = parse_args()
    #test = NCBIDataDownloaderPrep('AAA.csv')
    unittest.main()