import pandas as pd
import os
import sys
import numpy as np
import unittest

# TODO:Find better way to do the sys path thing.. this is bound to be problematic!!!
sys.path.append("../../../../NewGenes")



class HGTDBDataLoader():
    def __init__(self, data_type, partition_file):
        self.data_type = data_type
        self.partition_file = partition_file

    
    def _load_single_file(self, species, data_type, cross_partition = None, drop_na=False,return_df = False, genome_info = False):
        '''
        species: 
            choose specific species or use all
        
        data_type: 
            choose from [A,B,C,D,E] 
        
        cross_partition: 
            choose either train or test. This does 1 fold partition.
            first 2/3 is for training
            last 1/3 is for testing
        
        '''
        
        #if cross_partition is not None and return_folds is True: raise ValueError('Choose either cross partition OR folds!')

        preprocessed_path = "data/HGTDB/preprocessed_data"
        
        
        csv_list = os.listdir(preprocessed_path)
        
        if cross_partition is not None:
            #assert species == 'all', 'cross partitioning is only available when loading all genomes'
            if species != 'all': raise ValueError('cross partitioning is only available when loading all genomes')
            
            # partition 2:1 train, test
            if cross_partition == 'train':
                csv_list = csv_list[:int(len(csv_list)*(2/3))]
            elif cross_partition == 'test':
                csv_list = csv_list[int(len(csv_list)*(2/3)):]
            else:
                raise ValueError(f'No cross partition for {cross_partition}. Choose train or test')
        
        
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
            pass
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
        
        if return_df:
            # return as pandas
            return df
        else:
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
            if self.data_type == 'B':
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
        
        
        

class TestHGTDBDataLoaderPrep(unittest.TestCase):
    
    def test_dataloader(self):
        dataloader = HGTDBDataLoader('F','partition_file/HGTDB_firmicutes.csv')
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