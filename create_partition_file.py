import os
import argparse
import pandas as pd
import random

COMPLETE_DATASET = 'reference/ref_gtdb.csv'
#COMPLETE_DATASET = 'newgenes_complete_reference.csv'
ROOT_FOLDER = 'partition_file/'



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-l',
        '--level',
        type = str,
        default = 'phylum',
        help = 'Choose taxonomy level'
    )
    parser.add_argument(
        '-n',
        '--name',
        type = str,
        default = 'bacteroidota',
        help = 'name of chosen taxon'
    )
    return parser.parse_args()

def assert_args(args):
    try:
        taxonomy_table = pd.read_csv(COMPLETE_DATASET)
        no_link_or_taxon = taxonomy_table[taxonomy_table['GDTB_taxonomy'].isnull() | taxonomy_table['link'].isnull() ].index
        taxonomy_table.drop(no_link_or_taxon, inplace=True)
        taxonomy_table["phylum"]= [x.split(';')[1] for x in taxonomy_table["GDTB_taxonomy"].values]
        taxonomy_table["class"]= [x.split(';')[2] for x in taxonomy_table["GDTB_taxonomy"].values]
        taxonomy_table["order"]= [x.split(';')[3] for x in taxonomy_table["GDTB_taxonomy"].values]
        taxonomy_table["family"]= [x.split(';')[4] for x in taxonomy_table["GDTB_taxonomy"].values]
        taxonomy_table["genus"]= [x.split(';')[5] for x in taxonomy_table["GDTB_taxonomy"].values]
        taxonomy_table["species"]= [x.split(';')[6] for x in taxonomy_table["GDTB_taxonomy"].values]
        phylum = list(set( [phylum for phylum in taxonomy_table['phylum'].values ]))
        kelas = list(set( [kelas for kelas in taxonomy_table['class'].values ]))
        order = list(set( [order for order in taxonomy_table['order'].values ]))
        family = list(set( [family for family in taxonomy_table['family'].values ]))
        genus = list(set( [genus for genus in taxonomy_table['genus'].values ]))
        species = list(set( [species for species in taxonomy_table['species'].values ]))
    except FileNotFoundError:
        print(f"File {COMPLETE_DATASET} not found!")
        exit()
    
    taxonomy_levels = ['phylum','class','order','family','genus','species']
    if args.level not in taxonomy_levels:
        raise ValueError('taxonomy level unknown!')
    
    if args.level == 'phylum':
        query = 'p__'+str(args.name).capitalize()
        if query not in phylum:
            raise ValueError(f'given argument {args.name} is not found in phylum')
    elif args.level == 'class':
        query = 'c__'+str(args.name).capitalize()
        if query not in kelas :
            raise ValueError(f'given argument {args.name} is not found in class')
    elif args.level == 'order':
        query = 'o__'+str(args.name).capitalize()
        if query not in order:
            raise ValueError(f'given argument {args.name} is not found in order')
    elif args.level == 'family':
        query = 'f__'+str(args.name).capitalize()
        if query not in family:
            raise ValueError(f'given argument {args.name} is not found in family')
    elif args.level == 'genus':
        query = 'g__'+str(args.name).capitalize()
        if query not in genus:
            raise ValueError(f'given argument {args.name} is not found in class')
    elif args.level == 'species':
        query = 's__'+str(args.name).capitalize()
        if query not in species:
            raise ValueError(f'given argument {args.name} is not found in species')
    
    
    # return filtered taxonomy
    filtered_table = taxonomy_table[taxonomy_table[args.level] == query].copy()
    filtered_table.reset_index(inplace=True,drop=True)
    
    # filterd_table should be at least 4 to split 3:1
    if len(filtered_table)<4:
        raise ValueError('Total filtered query is less than 4! Not able to split train:test!')
    else:
        partition = ['train' if i<int(3/4 * len(filtered_table)) else 'test' for i in range(len(filtered_table)) ]
        random.shuffle(partition)
        filtered_table['partition'] = partition
        
    
    return args, filtered_table

    

    
if __name__ == '__main__':
    args, filtered_table = assert_args(parse_args())
    if not os.path.isdir(ROOT_FOLDER):
        os.makedirs(ROOT_FOLDER)
        print("creating folder : ", ROOT_FOLDER)
    
    # added another column called dl_link
    dl_link = []
    for i in range(len(filtered_table)):
        link = filtered_table.loc[i, 'link']
        link = link.replace('https:','ftp:')
        link = link+'/'+link.split('/')[-1]+'_cds_from_genomic.fna.gz'
        dl_link.append(link)
    filtered_table['dl_link'] = dl_link
    
    # print out filtered table as partition file
    filtered_table.to_csv(ROOT_FOLDER+f'{args.level}_{args.name}.csv', index=False)
    print('Done')