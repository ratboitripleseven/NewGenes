import pandas as pd
import os
import sys
import numpy as np
import unittest
import torch
from sklearn import preprocessing
from torch.nn.utils.rnn import pad_sequence
import math

# TODO:Find better way to do the sys path thing.. this is bound to be problematic!!!
sys.path.append("../../../../NewGenes")



class HGTDBDataLoader():
    def __init__(self, data_type, partition_file):
        self.data_type = data_type
        self.partition_file = partition_file
        self.null_count = 0
        self.na_count = 0

    
    def _load_single_file(self, species, data_type, drop_na=False, genome_info = False):
        '''
        species: 
            choose specific species or use all
        
        data_type: 
            choose from [A,B,C,D,E] 
                
        '''
        
        #if cross_partition is not None and return_folds is True: raise ValueError('Choose either cross partition OR folds!')

        preprocessed_path = "data/HGTDB/preprocessed_data"
        
        
        csv_list = os.listdir(preprocessed_path)
        
        
        if species == 'all':
            df = pd.DataFrame()
            for path in csv_list:
                # check if current path is a file
                if os.path.isfile(os.path.join(preprocessed_path, path)):
                    #list_csv_files.append(os.path.join(preprocessed_path, path))
                    df_temp = pd.read_csv(os.path.join(preprocessed_path, path), index_col='ID')
                    if genome_info:
                        genome = os.path.basename(path).replace(".csv", "")
                        col_genom = [genome for i in range(len(df_temp))]
                        df_temp['Genome'] = col_genom
                    df = pd.concat([df,df_temp])
        elif species is not None:
            csv_file = None
            for path in csv_list:
                # check if current path is a file
                if os.path.basename(path).replace(".csv", "") == species:
                    csv_file =os.path.join(preprocessed_path, path)
                    df = pd.read_csv(csv_file, index_col='ID')
            if csv_file is None:
                raise ValueError(f'{species} not found')
        else:
            raise ValueError(f'{species} not found')
            
        
        #df = read_csv(path, index_col='ID')
        # drop col
        #df.drop(columns=["Sim1","Sim2","Sim3","SimT","SimGC","SimMah", "Genename","COG","Function"],inplace=True)
        # Josefa drop this but in theory a gct = 0 can exist i.e. when both g and c counts are and only A and T exist. But maybe in biology?
        # from wikipedia GC content 100% or 0% is virtually impossible!
        #df.drop(df[df["GCT"] == 0].index,inplace=True)
        #df["HGT"].replace(2, 0, inplace=True)
        
        
        # set dataset according to data type
        if data_type == 'A':
            df = df.drop(columns=["FunctionCode","Strand","AADev","Length","SD1","SD2","SD3","SDT","Mah"]) # only GC1,GC2,GC3,GCT
        elif data_type == 'B':
            df = df.drop(columns="FunctionCode") # without Function code
        elif data_type == 'C':
            df = df.drop(columns=["FunctionCode","Strand","AADev"]) # without function code, Strand and AAdev
        elif data_type == 'D':
            df = df.drop(columns=["FunctionCode","Strand","AADev","Length","GC1","GC2","GC3","GCT", "Mah"]) # only SD 
        elif data_type == 'E':
            df = df.drop(columns=["FunctionCode","GC1","GC2","GC3","GCT"]) 
        elif data_type == 'F':
            df = df.drop(columns=["FunctionCode","Length","Strand","GC1","GC2","GC3","GCT"]) 
        else:
            raise ValueError(f'No data type {data_type}')
        
        # drop available nans and null
        if drop_na:
            df.dropna(inplace=True)
        
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

        
        # return as numpy array
        array = df.values
        # return X, y
        return array[:,0:-1], array[:,-1]
        
    
    def dataset_prep(self):
        '''
        TODO: Think about the design better!
        this abstract method should set X, Y, train and test
        '''
        
        partition_frame = pd.read_csv(self.partition_file)
        if self.data_type == 'A':
            columns = 4
        elif self.data_type == 'B':
            columns = 12
        elif self.data_type == 'C':
            columns = 10
        elif self.data_type == 'D':
            columns = 4
        elif self.data_type == 'E':
            columns = 8
        elif self.data_type == 'F':
            columns = 6
        else:
            raise ValueError(f'Unknown data type {self.data_type}')
        X_train = np.array([]).reshape(0,columns)
        y_train = np.array([]).reshape(0,)
        X_test = np.array([]).reshape(0,columns)
        y_test = np.array([]).reshape(0,)
        
        for i in range(len(partition_frame)):
            X,y = self._load_single_file(partition_frame.loc[i,'file'], self.data_type)
            if partition_frame.loc[i,'partition'] == 'train':
                X_train = np.concatenate([X, X_train], axis = 0)
                y_train = np.concatenate([y, y_train], axis = 0)
            else:
                X_test = np.concatenate([X, X_test], axis = 0)
                y_test = np.concatenate([y, y_test], axis = 0)
        
        return X_train, y_train, X_test, y_test
        
        
class HGTDBDatasetSequential(torch.utils.data.Dataset):
    def __init__(self, data_type, partition_file, partition):
        
        if partition not in ['train','test','valid']:
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
        
        
        self._init_dataset()
        

    
    def _load_single_file(self, species, data_type):
        '''
        replace this shit
        some of this shit is still legacy!!!
        '''

        preprocessed_path = "data/HGTDB/preprocessed_data"
        
        csv_list = os.listdir(preprocessed_path)
        
    
        if species is not None:
            csv_file = None
            for path in csv_list:
                # check if current path is a file
                if os.path.basename(path).replace(".csv", "") == species:
                    csv_file =os.path.join(preprocessed_path, path)
                    df = pd.read_csv(csv_file, index_col='ID')
            if csv_file is None:
                raise ValueError(f'{species} not found')
        else:
            raise ValueError(f'{species} not found')
            
        # set dataset according to data type
        if data_type == 'A':
            df = df.drop(columns=["FunctionCode","Strand","AADev","Length","SD1","SD2","SD3","SDT","Mah"]) # only GC1,GC2,GC3,GCT
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
            df = df.drop(columns=["FunctionCode","Strand","AADev","Length","GC1","GC2","GC3","GCT","Mah"]) # only SD1,SD2,SD3,SDT
        elif data_type == 'F':
            # So z score is basically what is calculated for sd1,2,3,t
            # z score for mah is added here
            # this z score is cutoff at two and scaled down by two so data ranges from -1 to 1
            df = df.drop(columns=["FunctionCode","Strand","AADev","Length","GC1","GC2","GC3","GCT"]) # only SD1,SD2,SD3,SDT, Mah
        
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
            

        
        
        #df.dropna(inplace=True)
        
        #after replacing nan
        #null_count = 0
        #na_count = 0
        #null_count +=df.isnull().sum().sum()
        #na_count +=df.isna().sum().sum()
        #print(null_count)
        #print(na_count)
        
        # preprocess i.e. normalize and shit
        # labels are not affected since there is only two options 0,1
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
        x = array[:,0:-1]
        y = array[:,-1]
        y = np.expand_dims(y, axis=1)
        #OHE = preprocessing.OneHotEncoder()
        #OHE.fit(y)
        #y = OHE.transform(y).toarray()
        return x,y 
        
    
    def _init_dataset(self):
        '''
        '''
        partition_frame = pd.read_csv(self.partition_file)
        partition_frame = partition_frame[partition_frame['partition']==self.partition].reset_index(drop=True)
        
        for i in range(len(partition_frame)):
            x,y = self._load_single_file(partition_frame.loc[i,'file'], self.data_type)
            
            if self.max_sequence<len(x):
                self.max_sequence = len(x)
                
            
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
        
        output = {
            "datum" : datum,
            "seq_length" : seq_length,
            "label" : label
        }
        return output
    
    def __len__(self):
        return len(self.data_x)
        
        
class HGTDBDatasetSequential_v2(HGTDBDatasetSequential):
    '''
    Inherits from HGTDBDatasetSequential
    This is for paritioning based on genes
    does not require info on how long the genomes are since no padding
    i.e. not processed by pad packed
    '''
    def __init__(self, data_type, partition_file, partition):
        super().__init__(data_type, partition_file, partition)
        
    def _init_dataset(self):
        '''
        '''
        partition_frame = pd.read_csv(self.partition_file)
        partition_frame = partition_frame[partition_frame['partition']==self.partition].reset_index(drop=True)
        
        x_set= []
        y_set = []
        
        self.max_sequence = 0
        
        for i in range(len(partition_frame)):
            x,y = self._load_single_file(partition_frame.loc[i,'file'], self.data_type)
        
                
            self.data_x.extend(torch.from_numpy(np.float32(x)))
            self.data_y.extend(torch.from_numpy(np.float32(y)))
            #self.data_seq_length.extend(torch.tensor(len(x)))
            
    def __getitem__(self, ind):
        '''
        torch data loader does not accept nones!
        '''
        datum = self.data_x[ind]
        label = self.data_y[ind]
        seq_length = None
        
        output = {
            "datum" : datum,
            #"seq_length" : seq_length,
            "label" : label
        }
        return output

            
        
                    

class TestHGTDBDataLoaderPrep(unittest.TestCase):
    
    def test_dataloader(self):
        dataloader = HGTDBDataLoader('F','partition_file/HGTDB_firmicutes.csv')
        x1,y1,x2,y2 = dataloader.dataset_prep()
        assert len(x1)!=0, "error!"
    
    def test_sequential_dataloader(self):
        hgtdb_train = HGTDBDatasetSequential('C','partition_file/HGTDB_firmicutes_trisplit.csv', 'train')
        #print(hgtdb_train[0])
        assert len(hgtdb_train)!=0, "error"
        
    def test_sequential_dataloader(self):
        hgtdb_train = HGTDBDatasetSequential('C','partition_file/HGTDB_firmicutes_trisplit.csv', 'train')
        dataloader = torch.utils.data.DataLoader(dataset=hgtdb_train,batch_size=5,shuffle=True)
        print(next(iter(dataloader)))
        assert len(hgtdb_train)!=0, "error"
        
    def test_sequential_dataloader_v2(self):
        hgtdb_train = HGTDBDatasetSequential_v2('C','partition_file/HGTDB_firmicutes_trisplit.csv', 'train')
        #print(hgtdb_train[0])
        print(f' Training length {len(hgtdb_train)}')
        assert len(hgtdb_train)!=0, "error"
    
    def test_sequential_dataloader_v2_2(self):
        hgtdb_train = HGTDBDatasetSequential_v2('C','partition_file/HGTDB_firmicutes_trisplit.csv', 'train')
        dataloader = torch.utils.data.DataLoader(dataset=hgtdb_train,batch_size=5,shuffle=True)
        
        print(next(iter(dataloader)))
        assert len(hgtdb_train)!=0, "error"
        
    def test_lengths(self):
        dataloader_1 = HGTDBDataLoader('F','partition_file/HGTDB_firmicutes.csv')
        x1,y1,x2,y2 = dataloader_1.dataset_prep()
        hgtdb_train = HGTDBDatasetSequential_v2('C','partition_file/HGTDB_firmicutes_trisplit.csv', 'train')
        #print(hgtdb_train[0])
        assert len(x1)==len(hgtdb_train), "error"
    
    
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