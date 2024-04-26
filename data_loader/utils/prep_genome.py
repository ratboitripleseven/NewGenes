'''
This script takes in .fna file and then return genes in dict
'''
import unittest
import os
import sys
sys.path.append("../../../../../NewGenes")


def prep_genome(file_path)-> dict:
    file_type = None
    if file_path.endswith('.fna'):
        file_type = 'fna'
        try:
            with open(os.path.join(file_path)) as f:
                lines = f.readlines()
        except IOError:
            print(f"ERROR: The file {os.path.join(file_path)} does not exist")
            return 0
    elif file_path.endswith('.fasta'):
        file_type = 'fasta'
        try:
            with open(os.path.join(file_path)) as f:
                lines = f.readlines()
        except IOError:
            print(f"ERROR: The file {os.path.join(file_path)} does not exist")
            return 0
    else:
        raise ValueError('file does not have the extension .fna or .fasta')
        
    genes = dict()
    for i in range(len(lines)):
        s = lines[i].strip()
        if s[0] == '>':
            if file_type == 'fna':
                # get locus_tag (gene)
                locus_tag_index = s.find('locus_tag=')
                if locus_tag_index == -1:
                    locus_tag = "No_Locust_Tag_{}".format(i)
                else:
                    locus_tag_end_index = s.find(']', locus_tag_index)
                    locus_tag = s[locus_tag_index+10:locus_tag_end_index]
                
                # get gene (synonym name)
                gene_index = s.find('gene=')
                if gene_index == -1:
                    gene = None
                else:
                    gene_end_index = s.find(']', gene_index)
                    gene = s[gene_index+5:gene_end_index]
                
                # get protein
                protein_index = s.find('protein=')
                protein_end_index = s.find(']', protein_index)
                protein = s[protein_index+8:protein_end_index]
                
                # get protein id
                protein_id_index = s.find('protein_id=')
                protein_id_end_index = s.find(']', protein_id_index)
                protein_id = s[protein_id_index+11:protein_id_end_index]
                
                # get location
                location_index = s.find('location=')
                location_end_index = s.find(']', location_index)
                location = s[location_index+9:location_end_index]
            elif file_type == 'fasta':
                locus_tag = "Tag_{}".format(i)
                gene = s
                protein = None
                protein_id = None
                location = None
            
            
            genes[locus_tag] = {
                'gene' : gene,
                'protein' : protein,
                'protein_id' : protein_id,
                'location': location,
                'g_count': None,
                'a_count': None,
                'c_count': None,
                't_count': None,
                'GC1': None,
                'SD1': None,
                'Sim1': None,
                'GC2': None,
                'SD2': None,
                'Sim2': None,
                'GC3': None,
                'SD3': None,
                'Sim3': None,
                'GCT': None,
                'SDT': None,
                'SimT': None,
                'SimGC': None,
                'Mah': None,
                'SimMah': None,
                'Dev.AA': None,
                'HGT':None,
                'rel_freq': None,
                '12_symbols': None,
                '48_symbols': None,
                'cub': None,
                'RSCU': None,
                'std_cub':{},
                'sequence':''
            }
        else:
            genes[locus_tag]['sequence'] += s
    return genes
            
   
            
    
        
    
        
    
        
class TestNCBIDataLoaderPrep(unittest.TestCase):
    def test_prep_fast_1(self):
        test = prep_genome('annotated_file_HGT/HGTDB_ALL/ASM221032v1.fasta')
        #print(test)
        assert test!=0, "error!"

    def test_prep_fast_2(self):
        test = prep_genome('data/NCBI/sequence_files/ASM668v1.fna')
        print(test)
        assert test!=0, "error!"



if __name__ == '__main__':
    unittest.main()