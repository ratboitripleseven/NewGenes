import torch
import pandas as pd
import os
from pandas import read_csv
from torch.utils.data import Dataset
import sys
sys.path.append("../../../../NewGenes")




## Single species
class SingleSpecies(Dataset):
    """Single species dataset"""

    def __init__(self, species, data_type,normalize=False):
        """
        Arguments:
            species: checks file string to get required species
            type: A,B,C
        """
        if data_type in ["A","B","C", "D", "E"]:
            self.data_type = data_type
        else:
            raise ValueError("Unknown data type. options = A, B, C ")
        
        self.csv_file = None
        self.preprocessed_path = "data/HGTDB/preprocessed_data"
        self.get_csv(species)
        
        
        self.normalize = normalize
        
        # load data
        self.load_csv()
        
        
    def get_csv(self, species):
        for path in os.listdir(self.preprocessed_path):
            # check if current path is a file
            if os.path.basename(path).replace(".csv", "") == species:
                self.csv_file=os.path.join(self.preprocessed_path, path)
            #if os.path.isfile(os.path.join(self.preprocessed_path, path)):
            #    self.list_csv_files.append(os.path.join(self.preprocessed_path, path))
        if self.csv_file is None:
            raise ValueError(f'{species} not found')
        

    def load_csv(self):
        #set dataframe
        df = read_csv(self.csv_file, index_col="ID")
        #print(ecoli_dataset_A.shape)
        if self.data_type == "B":
            df.drop(columns="FunctionCode", inplace=True) # without Function code
        elif self.data_type == "C":
            df.drop(columns=["FunctionCode","Strand","AADev"], inplace =True) # without function code, Strand and AAdev
        elif self.data_type == "D":
            df.drop(columns=["FunctionCode","Strand","AADev","Length","GC1","GC2","GC3","GCT", "Mah"], inplace=True) # just SD
        elif self.data_type == "E":
            df.drop(columns=["FunctionCode","GC1","GC2","GC3","GCT"], inplace = True) 
            
        # min max normalize
        if self.normalize:
            df=(df-df.min())/(df.max()-df.min())
        
        self.genes=df.to_numpy()

    
    def __len__(self):
        return len(self.genes)

    def __getitem__(self, idx):

        sample = self.genes[idx]
        sample = {'gc_signature': sample[0:-1], 'hgt': sample[-1]}

        return sample
    
    

## Single species
class AllSpecies(Dataset):
    """
    This dataset class take all in EDA/prepocessed_data
    and turn in into one single large dataset
    """

    def __init__(self, data_type, normalize=False, partition_type=None, partition=None):
        """
        Arguments:
            data_type: A,B,C,D,E
            
            normalize: def-> False; normalize value (recommended)
            
            partition_type: def-> None; holdout, fold1, fold2, fold3, fold4, fold5, fold6
            
            partition: partition -> None; train, test
            
            if both partition_type and partition is set to None -> gives back all data
        """
        if partition_type in [None, 'holdout', 'fold1','fold2','fold3','fold4','fold5','fold6']:
            self.partition_type = partition_type
        else:
            raise ValueError(f'partition type {partition_type} is unknown. Choose holdout or sixfold. Default None')
        
        if partition in [None, 'train', 'test']:
            self.partition = partition
        else:
            raise ValueError(f'partition type {partition} is unknown. train or test. Default is None')
        
        if self.partition_type is None and partition is not None:
            raise ValueError('Partition is given but not partition type. Please define Partition type')
        elif self.partition_type is not None and partition is None:
            raise ValueError('Partition_type is given but not partition. Please select which partition')
        
        if data_type in ["B","C","D", "E"]:
            self.data_type = data_type
        elif data_type in ["A"]:
            raise ValueError("A is still in progress")
        else:
            raise ValueError("Unknown data type. options = B, C, D")
        
        self.normalize = normalize
        
        self.preprocessed_path = "data/HGTDB/preprocessed_data"
        self.list_csv_files = []
        #self.list_csv_to_read = []
        
        # get all csvs
        self.create_csv_list()
        # load data
        self.load_csv()
        

    def create_csv_list(self):
        for path in os.listdir(self.preprocessed_path):
            # check if current path is a file
            if os.path.isfile(os.path.join(self.preprocessed_path, path)):
                self.list_csv_files.append(os.path.join(self.preprocessed_path, path))
        
        
        if self.partition_type == 'holdout':
            if self.partition == 'train':
                self.list_csv_files = self.list_csv_files[:int(len(self.list_csv_files)*(2/3))]
            elif self.partition == 'test':
                self.list_csv_files = self.list_csv_files[int(len(self.list_csv_files)*(2/3)):]
        elif self.partition_type == 'fold1':
            if self.partition == 'train':
                self.list_csv_files = self.list_csv_files[int(len(self.list_csv_files)*(1/6)):]
            elif self.partition == 'test':
                self.list_csv_files = self.list_csv_files[:int(len(self.list_csv_files)*(1/6))]
        elif self.partition_type == 'fold2':
            if self.partition == 'train':
                self.list_csv_files = self.list_csv_files[:int(len(self.list_csv_files)*(1/6))] + self.list_csv_files[int(len(self.list_csv_files)*(2/6)):]
            elif self.partition == 'test':
                self.list_csv_files = self.list_csv_files[int(len(self.list_csv_files)*(1/6)):int(len(self.list_csv_files)*(2/6))]
        elif self.partition_type == 'fold3':
            if self.partition == 'train':
                self.list_csv_files = self.list_csv_files[:int(len(self.list_csv_files)*(2/6))] + self.list_csv_files[int(len(self.list_csv_files)*(3/6)):]
            elif self.partition == 'test':
                self.list_csv_files = self.list_csv_files[int(len(self.list_csv_files)*(2/6)):int(len(self.list_csv_files)*(3/6))]
        elif self.partition_type == 'fold4':
            if self.partition == 'train':
                self.list_csv_files = self.list_csv_files[:int(len(self.list_csv_files)*(3/6))] + self.list_csv_files[int(len(self.list_csv_files)*(4/6)):]
            elif self.partition == 'test':
                self.list_csv_files = self.list_csv_files[int(len(self.list_csv_files)*(3/6)):int(len(self.list_csv_files)*(4/6))]
        elif self.partition_type == 'fold5':
            if self.partition == 'train':
                self.list_csv_files = self.list_csv_files[:int(len(self.list_csv_files)*(4/6))] + self.list_csv_files[int(len(self.list_csv_files)*(5/6)):]
            elif self.partition == 'test':
                self.list_csv_files = self.list_csv_files[int(len(self.list_csv_files)*(4/6)):int(len(self.list_csv_files)*(5/6))]
        elif self.partition_type == 'fold6':
            if self.partition == 'train':
                self.list_csv_files = self.list_csv_files[:int(len(self.list_csv_files)*(5/6))]
            elif self.partition == 'test':
                self.list_csv_files = self.list_csv_files[int(len(self.list_csv_files)*(5/6)):]


    def load_csv(self):
        # load all into one dataframe
        # print(self.list_csv_files)
        df = pd.DataFrame()
        for csv_file in self.list_csv_files:
            df_temp = read_csv(csv_file, index_col="ID")
            df = pd.concat([df, df_temp])
        
        
        #print(ecoli_dataset_A.shape)
        if self.data_type == "B":
            df.drop(columns="FunctionCode", inplace=True) # without Function code
        elif self.data_type == "C":
            df.drop(columns=["FunctionCode","Strand","AADev"], inplace =True) # without function code, Strand and AAdev
        elif self.data_type == "D":
            df.drop(columns=["FunctionCode","Strand","AADev","Length","GC1","GC2","GC3","GCT", "Mah"], inplace=True)
        elif self.data_type == 'E':
            df.drop(columns=["FunctionCode","GC1","GC2","GC3","GCT"], inplace = True) 
            
            
        #print(df.info(verbose=True))
        #print(df.describe())
        #print(df.isnull().values.any())
        #print(df.loc[df.isnull().any(axis=1)])
        df.dropna(inplace=True)
        #print(df.min())
        #print(df.max())
        # min max normalize
        if self.normalize:
            df=(df-df.min())/(df.max()-df.min())
        
        
        #it seems that there is one row that has nan value and messes up my damn network
        #df.dropna(inplace=True)
        self.genes=torch.from_numpy(df.to_numpy())

    
    def __len__(self):
        return len(self.genes)

    def __getitem__(self, idx):

        sample = self.genes[idx]
        sample = {'gc_signature': sample[0:-1], 'hgt': sample[-1]}

        return sample
    
    
    

if __name__ == '__main__':
    test_single = SingleSpecies(species="ecoli",data_type="B")
