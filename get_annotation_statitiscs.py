'''
This script only to be used for reporting
'''
import unittest
import os
import sys
import json



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
        
    
def create_annotation_statistics(annotated_file, result_file):
    genes_dict = open_annotated_fasta(annotated_file)
    arg_genes_json = open_result_json(result_file)
    
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
    
    return 0
    
    
    
    

        
class TestAnnotationStatistics(unittest.TestCase):
    def test_open_annotated_fasta(self):
        test = open_annotated_fasta('annotated_file_HGT/HGTDB_ALL/ASM221032v1.fasta')
        #print(test['B7988_00005'])
        assert test['B7988_00005']!=0, "error!"
        
    def test_open_result_json(self):
        test = open_result_json('results/HGT_DB/analysis/card/ASM221032v1/ASM221032v1.json')
        
        #print(test['B7988_12605'])
        assert True, "error!"
        
    def test_create_annotation_statistics(self):
        test = create_annotation_statistics('annotated_file_HGT/HGTDB_ALL/ASM221032v1.fasta','results/HGT_DB/analysis/card/ASM221032v1/ASM221032v1.json')
        
        assert test == 0, 'error!'
        



if __name__ == '__main__':
    #unittest.main()
    create_annotation_statistics('annotated_file_HGT/HGTDB_ALL/ASM221032v1.fasta','results/HGT_DB/analysis/card/ASM221032v1/ASM221032v1.json')