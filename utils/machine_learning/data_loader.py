import pandas as pd
import os

def ml_load_species(species, data_type, cross_partition = None, drop_na=False,return_df = False, genome_info = False):
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

    preprocessed_path = "../EDA/preprocessed_data"
    
    
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
        
        
def ml_load_species_folds(data_type, drop_na=False,return_df = False, genome_info = False):
    preprocessed_path = "../EDA/preprocessed_data"
    csv_list = os.listdir(preprocessed_path)
    
    fold_01 = csv_list[:int(len(csv_list)*(1/6))]
    fold_02 = csv_list[int(len(csv_list)*(1/6)):int(len(csv_list)*(2/6))]
    fold_03 = csv_list[int(len(csv_list)*(2/6)):int(len(csv_list)*(3/6))]
    fold_04 = csv_list[int(len(csv_list)*(3/6)):int(len(csv_list)*(4/6))]
    fold_05 = csv_list[int(len(csv_list)*(4/6)):int(len(csv_list)*(5/6))]
    fold_06 = csv_list[int(len(csv_list)*(5/6)):]
    
    df_dict = {}
    counter = 1
    for fold in [ fold_01, fold_02, fold_03, fold_04, fold_05, fold_06]:
        df = pd.DataFrame()
        for path in fold:
            # check if current path is a file
            if os.path.isfile(os.path.join(preprocessed_path, path)):
                #list_csv_files.append(os.path.join(preprocessed_path, path))
                df_temp = pd.read_csv(os.path.join(preprocessed_path, path), index_col='ID')
                if genome_info:
                    genome = os.path.basename(path).replace(".csv", "")
                    col_genom = [genome for i in range(len(df_temp))]
                    df_temp['Genome'] = col_genom
                df = pd.concat([df,df_temp])
                
                
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
        else:
            raise ValueError(f'No data type {data_type}')
        
        if drop_na:
            df.dropna(inplace=True)
        
        df_dict[f'fold_{counter}'] = df
        counter+=1
        
    
    if return_df: 
        return df_dict['fold_1'], df_dict['fold_2'], df_dict['fold_3'], df_dict['fold_4'], df_dict['fold_5'], df_dict['fold_6']
    else:
        array_1 = df_dict['fold_1'].values
        array_2 = df_dict['fold_2'].values
        array_3 = df_dict['fold_3'].values
        array_4 = df_dict['fold_4'].values
        array_5 = df_dict['fold_5'].values
        array_6 = df_dict['fold_6'].values
        
        
        return (array_1[:,0:-1], array_1[:,-1]), (array_2[:,0:-1], array_2[:,-1]), (array_3[:,0:-1], array_3[:,-1]), (array_4[:,0:-1], array_4[:,-1]), (array_5[:,0:-1], array_5[:,-1]), (array_6[:,0:-1], array_6[:,-1])
        
    