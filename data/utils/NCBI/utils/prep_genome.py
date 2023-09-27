'''
This script takes in .fna file and then return genes in dict
'''
import os
import sys
sys.path.append("../../../../../NewGenes")
SEQUENCES_FOLDER = 'data/NCBI/sequence_files'

def prep_genome(sequence_folder,genome)-> dict:
    if not os.path.isdir(sequence_folder):
        sys.exit("The folder sequence_files does not exist!")
    
    try:
        with open(os.path.join(sequence_folder,genome+'.fna')) as f:
            lines = f.readlines()
    except IOError:
        print(f"ERROR: The file {os.path.join(sequence_folder,genome)}.fna does not exist")
        return 0
        
    genes = dict()
    for i in range(len(lines)):
        s = lines[i].strip()
        if s[0] == '>':
            # get locus_tag (gene)
            locus_tag_index = s.find('locus_tag=')
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
                'GC2': None,
                'SD2': None,
                'GC3': None,
                'SD3': None,
                'GCT': None,
                'SDT': None,
                'rel_freq': None,
                '12_symbols': None,
                '48_symbols': None,
                'cub': None,
                'std_cub':{},
                'sequence':''
            }
        else:
            genes[locus_tag]['sequence'] += s
    return genes
            
   
            
    
        
    
        
    
        
    




if __name__ == '__main__':
    genome = prep_genome(SEQUENCES_FOLDER,'bsub')
    for i in genome:
        print(i)