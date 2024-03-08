import os
import random
from collections import deque
import csv
import argparse

ROOT_FOLDER = 'data/HGTDB/preprocessed_data'
SAVE_FOLDER = 'partition_file/'

'''
This script partition the whole hgtdb into 6 parts

partition into train:valid:test (4:1:1)
'''

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-n',
        '--name',
        type = str,
        default = 'HGTDB_CV',
        help = 'Name of cv'
    )
    parser.add_argument(
        '--valid', 
        dest='flag', 
        action='store_true',
        help='Set the flag value to True.'
    )
    parser.add_argument(
        '--no-valid', 
        dest='flag', 
        action='store_false',
        help='Set the flag value to False.'
    )
    parser.set_defaults(flag=True)
    return parser.parse_args()

def main(name, valid):
    if not os.path.isdir(ROOT_FOLDER):
        print('ROOT folder not found')
        print('Please download HGTDB')
        exit()
        
    list_of_files = os.listdir(ROOT_FOLDER)
    
    block_length = int(len(list_of_files)/6)
    
    print(f'Partitioning into blocks of {6}')
    print(f'Each block has {block_length}')
    
    # shuffles list
    random.shuffle(list_of_files)
    
    # total blocks = 6
    block_indices = [0,1,2,3,4,5]
    
    temp_holder = []
    
    count = 0 
    count_block = 0
    for i in list_of_files:
        temp_holder.append((i,count_block))
        if count < block_length-1:
            # print(count)
            count+=1
        else:
            count = 0
            count_block+=1

    # check_dist(temp_holder)
    print(name)
    print(valid)
    tag_distribution(temp_holder, block_indices,name,valid)
    

def tag_distribution(holder, block_indices, name, valid):
    indices = deque(block_indices)
    for i in range(len(block_indices)):
        partition_lod = []
        for files,block in holder:
            if indices[0] == block:
                temp_dict = {
                    'file': files.replace('.csv', ''),
                    'partition': 'test'
                }
            elif indices[1] == block:
                if valid:
                    partition='valid'
                else:
                    partition='train'
                temp_dict = {
                    'file': files.replace('.csv', ''),
                    'partition': partition
                }
            else:
                temp_dict = {
                    'file': files.replace('.csv', ''),
                    'partition': 'train'
                }
            partition_lod.append(temp_dict)
            
        with open(f'{SAVE_FOLDER}{name}_{i}.csv', 'w', newline='') as output_file:
            dict_writer = csv.DictWriter(output_file, ['file','partition'])
            dict_writer.writeheader()
            dict_writer.writerows(partition_lod)
            
        #print(partition_lod)
        indices.rotate(1)
        #print(indices[0])
        #break
        
    

    
    
        
        
    
    
    
def check_dist(holder):
    block_0 =0
    block_1 =0
    block_2 =0
    block_3 =0
    block_4 =0
    block_5 =0
    
    for item in holder:
        if item[1] == 0:
            block_0 += 1
        if item[1] == 1:
            block_1 += 1
        if item[1] == 2:
            block_2 += 1
        if item[1] == 3:
            block_3 += 1
        if item[1] == 4:
            block_4 += 1
        if item[1] == 5:
            block_5 += 1
            
    print(f'block_0 {block_0}')
    print(f'block_1 {block_1}')
    print(f'block_2 {block_2}')
    print(f'block_3 {block_3}')
    print(f'block_4 {block_4}')
    print(f'block_0 {block_5}')
    
            
    print(holder)

if __name__ == '__main__':
    args = parse_args()
    print(args.flag)
    main(args.name, args.flag)
    