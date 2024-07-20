'''
This script only to be used for reporting
'''
import unittest
import os
import sys
import json
import csv
import pandas as pd
import argparse

ANNOTATED_PREFIX = 'annotated_file_HGT/'
SNAKEMAKE_RESULT_PREFIX = 'results/HGT_DB/analysis/card/'


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-m',
        '--model',
        type = str,
        default = None,
        help = 'model that was used to annotate'
    )
    parser.add_argument(
        '-p',
        '--partition',
        type = str,
        default = None,
        help = 'path to partition file'
    )
    return parser.parse_args()

def assert_config(args):
    if args.model is None:
        raise ValueError('Please input model!')
    if args.partition is None:
        raise ValueError('Please add path to partition file! make sure it is RTR')
    
    if not os.path.isdir(ANNOTATED_PREFIX + args.model):
        raise ValueError(f' model {args.model} not found')
    if not os.path.exists(args.partition):
        raise ValueError(f' partition {args.partition} not found')
    
    return args


def open_annotated_fasta(file_path)-> dict:
    file_type = None
    if file_path.endswith('.fasta'):
        try:
            with open(os.path.join(file_path)) as f:
                lines = f.readlines()
        except IOError:
            print(f"ERROR: The file {os.path.join(file_path)} does not exist")
            return 0
        
    genes = dict()
    for i in range(len(lines)):
        s = lines[i].strip()
        if s[0] == '>':
            gene = s.replace('>','')
            gene = gene.split(' ')[0]
            HGT_start_index = s.find('HGT:')
            HGT_end_index = s.find(']', HGT_start_index)
            HGT = s[HGT_start_index+4:HGT_end_index]
            
            genes[gene] = {
                'HGT' : int(HGT),
                'sequence':''
            }
        else:
            genes[gene]['sequence'] += s
            
    return genes
            
   

            
def open_result_json(file_path):
    with open(file_path) as payload:
        json_payload = json.load(payload)
    return json_payload
        
    
def create_annotation_statistics(model, name, csv_file):
    
    fasta_file = ANNOTATED_PREFIX +model+'/'+ name + '.fasta'
    json_file = SNAKEMAKE_RESULT_PREFIX + name +'/'+ name+'.json'
    
    #print(fasta_file)
    #print(json_file)
    genes_dict = open_annotated_fasta(fasta_file)
    arg_genes_json = open_result_json(json_file)
    hgtdb_csv = pd.read_csv(csv_file,index_col=0)
    
    
    available_genes = list(genes_dict.keys())
    arg_genes = list(arg_genes_json.keys())
    
    available_hgts = [i for i in available_genes if genes_dict[i]['HGT'] == 1]
    
    # this could be redundant but I am making sure the ARG pipelone is alright
    list_of_arg_as_hgt = []
    for i in arg_genes:
        if i in available_hgts:
            list_of_arg_as_hgt.append(i)
            
    
    print(f'Total available genes {len(available_genes)}')
    print(f'Total amount of HGTs {len(available_hgts)}')
    print(f'Total amount of ARG: {len(arg_genes_json)}')
    print(f'Total amount of HGTs that are hgt: {len(list_of_arg_as_hgt)}')
    
    #return 0
    #return genome, available_genes, available_hgts, arg_genes, list_of_arg_as_hgt
    temp_dict = {
        'available_genes':len(available_genes),
        'available_hgts':available_hgts,
        'arg_genes':arg_genes,
        'list_of_arg_as_hgt':list_of_arg_as_hgt
    }
    
    #stats = stats_dict[genome]
    #list_of_hgts = stats['available_hgts']
    list_of_hgts = available_hgts
    #list_of_args = stats['arg_genes']
    list_of_args = arg_genes
    #list_of_args_hgts = stats['list_of_arg_as_hgt']
    list_of_args_hgts = list_of_arg_as_hgt
    
    to_csv = []
    for genes in hgtdb_csv.index:
        sim1,sim2,sim3,simt = hgtdb_csv.loc[genes][['Sim1','Sim2','Sim3','SimT']].values
        if sim1 * sim3 > 0 and (abs(sim1)>= 1 or abs(sim3)>= 1):
                extraneous = True
        elif abs(simt) >= 1:
                extraneous = True
        else:
                extraneous = False
        if genes not in list_of_hgts and genes not in list_of_args:
            to_csv.append(('None',genes,sim1,sim2,sim3,simt, extraneous))
    if len(list_of_hgts)>0:
        for hgt in list_of_hgts:
            #print(hgt)
            #sd1,sd2,sd3,sdt = hgtdb_csv.loc[hgt][['SD1','SD2','SD3','SDT']].values
            sim1,sim2,sim3,simt = hgtdb_csv.loc[hgt][['Sim1','Sim2','Sim3','SimT']].values
            if sim1 * sim3 > 0 and (abs(sim1)>= 1 or abs(sim3)>= 1):
                extraneous = True
            elif abs(simt) >= 1:
                extraneous = True
            else:
                extraneous = False
            if hgt not in list_of_args_hgts:
                to_csv.append(('HGT',hgt,sim1,sim2,sim3,simt, extraneous))
    if len(list_of_args)>0:
        for arg in list_of_args:
            #print(arg)
            #sd1,sd2,sd3,sdt = hgtdb_csv.loc[arg][['SD1','SD2','SD3','SDT']].values
            sim1,sim2,sim3,simt = hgtdb_csv.loc[arg][['Sim1','Sim2','Sim3','SimT']].values
            if sim1 * sim3 > 0 and (abs(sim1)>= 1 or abs(sim3)>= 1):
                extraneous = True
            elif abs(simt) >= 1:
                extraneous = True
            else:
                extraneous = False
            if arg not in list_of_args_hgts:
                to_csv.append(('ARG',arg,sim1,sim2,sim3,simt, extraneous))
    if len(list_of_args_hgts)>0:
        for arg_hgt in list_of_args_hgts:
            #print(arg)
            #sd1,sd2,sd3,sdt = hgtdb_csv.loc[arg][['SD1','SD2','SD3','SDT']].values
            sim1,sim2,sim3,simt = hgtdb_csv.loc[arg_hgt][['Sim1','Sim2','Sim3','SimT']].values
            if sim1 * sim3 > 0 and (abs(sim1)>= 1 or abs(sim3)>= 1):
                extraneous = True
            elif abs(simt) >= 1:
                extraneous = True
            else:
                extraneous = False
            to_csv.append(('ARG_HGT',arg_hgt,sim1,sim2,sim3,simt, extraneous))
            
    with open('report_'+model+'_'+name+'.csv','w') as out:
        csv_out=csv.writer(out)
        #csv_out.writerow(['genome','type','gene','SD1','SD2','SD3','SDT'])
        csv_out.writerow(['type','gene','Sim1','Sim2','Sim3','SimT','extraneous'])
        csv_out.writerows(to_csv)
    
    

        

def open_partition_file(partition_file):
    dataframe = pd.read_csv(partition_file)
    return dataframe['file_path'].tolist()


if __name__ == '__main__':
    #unittest.main()
    #create_annotation_statistics('annotated_file_HGT/000000_release_HGBC/ASM2013877v1.fasta','results/HGT_DB/analysis/card/ASM2013877v1/ASM2013877v1.json', 'sample_sequences/test/ASM2013877v1.csv')
    #csv_files = open_partition_file('sample_sequences/test/test_partition_RTR.csv')
    #model = '000000_release_HGBC'
    args = assert_config(parse_args())
    csv_files = open_partition_file(args.partition)
    model = args.model
    for i in csv_files:
        sample_name = i.split('/')[-1]
        sample_name = sample_name.replace('.csv','')
        print(sample_name)
        create_annotation_statistics(model, sample_name, i)
        #print(i.split('/')[-1])
    