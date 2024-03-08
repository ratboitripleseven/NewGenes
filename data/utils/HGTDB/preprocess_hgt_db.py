import os
import pandas as pd
import numpy as np
import sys
sys.path.append("../../../../NewGenes")

ROOT_FOLDER = "data/HGTDB/temp_data"
OUTPUT_FOLDER = "data/HGTDB/preprocessed_data"

'''
TODO: Ask Josefa regarding Function Code.. It seems she does some checking.. Or check it yourself
NOTE: For now Type A data (with Function Code) is not possible
NOTE: Josefa's bsub has one less entry (resolved by adding Josefa's testing function)
NOTE: have equal entries with hinf, ecoli2, ecoli

'''

def preprocess(path_to_file):
    df = pd.read_csv(path_to_file, sep='\t', index_col="Synonym")
    
    to_drop = ["Coordinates", "GCRegion", "PID", "Gene name", "SimMah", "SimT", "Sim1", "Sim2", "Sim3","Function", "COG","SimGC"]
    df.drop(columns=to_drop, inplace=True)
    
    #rename
    #df.rename(columns={"Gene name": "Genename", "Dev.AA": "AADev"}, inplace=True)
    # rename some columns & index
    df.rename(columns={ "Dev.AA": "AADev"}, inplace=True)
    df.index.names = ["ID"]
    
    df_obj = df.select_dtypes(['object'])
    df[df_obj.columns] = df_obj.apply(lambda x: x.str.strip())

    # HGT column - target column
        # reorder HGT column to the end
    cols = df.columns.tolist()
    cols.append(cols.pop(cols.index("HGT")))
    df = df[cols]
    # binary encoding HGT: HGT: 1, non-HGT: 0, excluded from analysis: 0
    # TODO: ask Josefa whether to delete f or not
    hgt_dict = {"":0, "H":1, "f": 0}
    # Trouble when no f is found
    # print('{} genes were excluded from analysis'.format(df["HGT"].value_counts()["f"]))
    print('{} genes were excluded from analysis'.format(df[df['HGT'] == 'f']['HGT'].count()))
    df.replace({"HGT":hgt_dict}, inplace=True)

    # binary encoding strand
    strand_dict = {'+':0, '-':1}
    df.replace({"Strand":strand_dict}, inplace=True)

    # change missing values to naN
    #df["Genename"] = np.where(df["Genename"] == "-", df.index, df["Genename"])
    #df.replace("", np.nan, inplace=True)
    #df.replace("-", np.nan, inplace=True)
    df["FunctionCode"].replace('-','',inplace=True)
    # binary encoding AADev
    df["AADev"].replace(np.nan, 0, inplace=True)
    df["AADev"].replace('', 0, inplace=True)
    df["AADev"] = np.where(df["AADev"] != 0, 1, df["AADev"])
    
    
    # taken from Josefa's correct strain
    gc_df = df.loc[df["GCT"] == 0]
    print("GCT == 0:")
    if len(gc_df.index) >= 1:
        print(gc_df.index.to_list())
        # deletes all non-HGT genes with GCT = 0 & prints a warning otherwise
        for gene in gc_df.index.to_list():
            if gc_df.loc[gene]["HGT"] == 0:
                df.drop(index=gene, inplace=True)
                print("deleted {}".format(gene))
            else:
                print("Warning: a HGT gene as an abnormal GCT of 0!")
    else:
        print("none")

    # if an index is duplicate
    print("duplicated index?")
    idx=pd.Index(df.index)#.values)
    duplicate_lst = idx.duplicated(keep=False)
    if True in duplicate_lst:
        print("There are duplicates! \nchanged to:")
            # list of rownumber whose indices are doubled
        duplicate_row_lst = np.where(duplicate_lst == True)[0].tolist()
        #print(duplicate_row_lst)
        duplicate_ind_lst = idx.to_list()
            # correct indices:
        first = True
        unique_ID = 0
        for row in duplicate_row_lst:
            if first:
                new_ind = "{}a_{}".format(duplicate_ind_lst[row], unique_ID)
                idx = idx.delete(row).insert(row, new_ind)
                first = False
            else:
                new_ind = "{}b_{}".format(duplicate_ind_lst[row], unique_ID)
                idx = idx.delete(row).insert(row, new_ind)
                first = True
            print(new_ind)
            unique_ID +=1
            # set corrected indices as new indices for the df
        df.set_index(idx, inplace=True)#pd.Series(duplicate_ind_lst), inplace=True).rename_axis("ID") #set_axis(duplicate_ind_lst,axis=0, inplace=True)
    else:
        print("no")
    
    return df
    
    

if __name__ == '__main__':
    # get downloaded files
    downloaded_files = []
    # Iterate directory
    for path in os.listdir(ROOT_FOLDER):
        # check if current path is a file
        if os.path.isfile(os.path.join(ROOT_FOLDER, path)):
            downloaded_files.append(os.path.join(ROOT_FOLDER, path))
            
    
    # create output folder
    if not os.path.isdir(OUTPUT_FOLDER):
        os.makedirs(OUTPUT_FOLDER)
        print("created folder : ", OUTPUT_FOLDER)
    else:
        print(OUTPUT_FOLDER, "folder already exists.")
    
    for files in downloaded_files:
        print(f"Proccessing {files}")
        df = preprocess(files)
        df.to_csv(files.replace(".tsv", ".csv").replace(ROOT_FOLDER,OUTPUT_FOLDER), sep=",")
        print("Completed")
