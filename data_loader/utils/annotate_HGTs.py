import sys
sys.path.append("../../NewGenes")
import math 
from Bio.SeqUtils import GC123
import numpy as np
from random import choices,choice
from data_loader.metrics.calc_cub import calc_cub

######### 

# GC CONTENT

#########
def calc_GC_deviation(genome, gene):
    # DEPRECATE THIS
    dev_GCT = genome[gene]['GCT'] - genome.mean_GCT
    dev_GC1 = genome[gene]['GC1'] - genome.mean_GC1
    dev_GC2 = genome[gene]['GC2'] - genome.mean_GC2
    dev_GC3 = genome[gene]['GC3'] - genome.mean_GC3
    equal_sign_check = dev_GC1*dev_GC3
    #populate stats
    if abs(dev_GC1) > (1.5*genome.std_GC1):
        if abs(dev_GC1) > (2*genome.std_GC1):
            if dev_GC1 > 0:
                genome[gene]['Sim1'] = 2
            else:
                genome[gene]['Sim1'] = -2
        else:
            if dev_GC1 > 0:
                genome[gene]['Sim1'] = +1
            else:
                genome[gene]['Sim1'] = -1
    else:
        genome[gene]['Sim1'] = 0
    
    if abs(dev_GC2) > (1.5*genome.std_GC2):
        if abs(dev_GC2) > (2*genome.std_GC2):
            if dev_GC2 > 0:
                genome[gene]['Sim2'] = +2
            else:
                genome[gene]['Sim2'] = -2
        else:
            if dev_GC2 > 0:
                genome[gene]['Sim2'] = +1
            else:
                genome[gene]['Sim2'] = -1
    else:
        genome[gene]['Sim2'] = 0
        
    if abs(dev_GC3) > (1.5*genome.std_GC3):
        if abs(dev_GC3) > (2*genome.std_GC3):
            if dev_GC3 > 0:
                genome[gene]['Sim3'] = +2
            else:
                genome[gene]['Sim3'] = -2
        else:
            if dev_GC3 > 0:
                genome[gene]['Sim3'] = +1
            else:
                genome[gene]['Sim3'] = -1
    else:
        genome[gene]['Sim3'] = 0
    
    if abs(dev_GCT) > (1.5*genome.std_GCT):
        if abs(dev_GCT) > (2*genome.std_GCT):
            if dev_GCT > 0:
                genome[gene]['SimT'] = +2
            else:
                genome[gene]['SimT'] = -2
        else:
            if dev_GC2 > 0:
                genome[gene]['SimT'] = +1
            else:
                genome[gene]['SimT'] = -1
    else:
        genome[gene]['SimT'] = 0
    
    return dev_GCT, dev_GC1, dev_GC3, equal_sign_check
    
def GC_Content_Deviation(genome, return_genomic_strips = False):
    # Init list of extraneous genes
    list_of_extraneous_genes = []
    
    # part 1: GCT, GC1, GC3 deviation
    for gene in genome.genes:        
        # consider only more than 300bp
        # either gct > 1.5
        # or sign gc1 and gc3 equal and one of is > 1.5
        if len(genome[gene]['sequence']) > 300:
            if genome[gene]['SimGC'] >0:
                list_of_extraneous_genes.append(gene)
        else:
            # annotate filtered genes
            genome[gene]['HGT']='f'
    
    # part two: 11 genes window
    list_of_genes = list(genome.genes.keys())
    genomic_strips = []
    
    for k in range(len(list_of_genes)-10):
        window = {}
        j = 0
        while j < 11:
            # get window
            locust_tag = list_of_genes[k + j]
            
            # take genes that are more than 300bp
            if len(genome[locust_tag]['sequence'])>300:
                data = genome[locust_tag]
                window[locust_tag] = data
            # iterate        
            j+=1
            
        # count total extraneous genes in window
        extraneous_counter = 0
        for l in window.keys():
            if l in list_of_extraneous_genes:
                extraneous_counter+=1
        
        # check windows with more than or
        # equal to 5 extraneous genes
        if extraneous_counter >=5:
            # add their sequences together
            window_sequences = ''
            for m in window.keys():
                window_sequences += window[m]['sequence']
            
            # get standard deviation of strip
            GCT, GC1, GC2, GC3 = GC123(window_sequences)
            SDT = (GCT - genome.mean_GCT)/genome.std_GCT
            SD1 = (GC1- genome.mean_GC1)/genome.std_GC1
            SD2 = (GC2 - genome.mean_GC2)/genome.std_GC2
            SD3 = (GC3 - genome.mean_GC3)/genome.std_GC3
            
            # tag genes as extraneous if they have equal deviation to the its strip
            for n in window.keys():
                # check only genes not in current list
                if n not in list_of_extraneous_genes:
                    check_SDT = window[n]['SDT']*SDT
                    check_SD1 = window[n]['SD1']*SD1
                    check_SD2 = window[n]['SD2']*SD2
                    check_SD3 = window[n]['SD3']*SD3
                    
                    if (check_SDT > 0) and (check_SD1 > 0) and (check_SD2 > 0) and (check_SD3 > 0):
                        list_of_extraneous_genes.append(n)
                
            
            genomic_strips.append(window)
    
    if return_genomic_strips:
        return list_of_extraneous_genes, genomic_strips
    else:
        return list_of_extraneous_genes
    



######### 

# AA CONTENT

#########
CODE = {
    'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
    'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C',
    'tta': 'L', 'tca': 'S', 'taa': '*', 'tga': '*',
    'ttg': 'L', 'tcg': 'S', 'tag': '*', 'tgg': 'W',
    'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R',
    'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
    'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R',
    'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R',
    'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S',
    'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S',
    'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R',
    'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',
    'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G',
    'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G',
    'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
    'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G'
}
LIST_OF_AMINO_ACID = set([CODE[cds] for cds in CODE])
LIST_OF_AMINO_ACID = {v: [cds.upper() for cds in CODE if CODE[cds]==v] for v in LIST_OF_AMINO_ACID}

def calculate_amino_acid_content_genome_mean(genome):
    amino_acid_count_genome = {amino_acid: 0 for amino_acid in LIST_OF_AMINO_ACID}
    for locust_tag in genome.genes:
        for amino_acid in LIST_OF_AMINO_ACID:
            for codons in LIST_OF_AMINO_ACID[amino_acid]:
                amino_acid_count_genome[amino_acid] += genome[locust_tag]['cub'][codons]

    total_amino_acid_count = sum([amino_acid_count_genome[amino_acid] for amino_acid in amino_acid_count_genome ])
    for i in amino_acid_count_genome:
        amino_acid_count_genome[i]=amino_acid_count_genome[i]*100/total_amino_acid_count

    genome.mean_AA = amino_acid_count_genome
    
def calculate_amino_acid_content_gene_mean(genome):
    # calculate mean
    for locust_tag in genome.genes:
        amino_acid_count_gene = {amino_acid: 0 for amino_acid in LIST_OF_AMINO_ACID}
        for amino_acid in LIST_OF_AMINO_ACID:
            for codons in LIST_OF_AMINO_ACID[amino_acid]:
                amino_acid_count_gene[amino_acid]  += genome[locust_tag]['cub'][codons]
                
        total_amino_acid_count = sum([amino_acid_count_gene[amino_acid] for amino_acid in amino_acid_count_gene ])
        for AA in amino_acid_count_gene:
            amino_acid_count_gene[AA] = amino_acid_count_gene[AA]*100/total_amino_acid_count

        genome[locust_tag]['AA_Content_mean']=amino_acid_count_gene
        
def calculate_amino_acid_content_genome_std(genome):
    amino_acid_std = {amino_acid: 0 for amino_acid in LIST_OF_AMINO_ACID}
    
    for amino_acid in LIST_OF_AMINO_ACID:
        AA_mean = genome.mean_AA[amino_acid]
        sum_diff = 0
        for locust_tag in genome.genes:
            diff = (genome[locust_tag]['AA_Content_mean'][amino_acid]) - AA_mean
            diff = diff*diff
            sum_diff+=diff
        
        amino_acid_std[amino_acid]=math.sqrt(sum_diff/len(genome))
            
    genome.std_AA = amino_acid_std
    
def check_amino_acid_deviation(genome):
    list_of_non_extraneous_genes_AA = []
    for locust_tag in genome.genes:
        list_of_dev_AA = []
        for amino_acid in LIST_OF_AMINO_ACID:
            devAA = genome[locust_tag]['AA_Content_mean'][amino_acid] - genome.mean_AA[amino_acid]
            if devAA > 3*genome.std_AA[amino_acid]:
                if locust_tag in list_of_non_extraneous_genes_AA:
                    pass
                else:
                    list_of_non_extraneous_genes_AA.append(locust_tag)
                list_of_dev_AA.append((amino_acid, devAA/genome.std_AA[amino_acid]))
                #break
        if list_of_dev_AA is None:
            genome[locust_tag]['Dev.AA'] = None
        else:
            genome[locust_tag]['Dev.AA'] = list_of_dev_AA
    return list_of_non_extraneous_genes_AA


######### 

# MAHALANOBIS DISTANCE

#########
CODE_COVARMAT = {
    'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
    'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C',
    'tta': 'L', 'tca': 'S', 'ttg': 'L', 'tcg': 'S',
    'tgg': 'W', 'ctt': 'L', 'cct': 'P', 'cat': 'H',
    'cgt': 'R', 'ctc': 'L', 'ccc': 'P', 'cac': 'H',
    'cgc': 'R', 'cta': 'L', 'cca': 'P', 'caa': 'Q',
    'cga': 'R', 'ctg': 'L', 'ccg': 'P', 'cag': 'Q',
    'cgg': 'R', 'att': 'I', 'act': 'T', 'aat': 'N',
    'agt': 'S', 'atc': 'I', 'acc': 'T', 'aac': 'N',
    'agc': 'S', 'ata': 'I', 'aca': 'T', 'aaa': 'K',
    'aga': 'R', 'atg': 'M', 'acg': 'T', 'aag': 'K',
    'agg': 'R', 'gtt': 'V', 'gct': 'A', 'gat': 'D',
    'ggt': 'G', 'gtc': 'V', 'gcc': 'A', 'gac': 'D',
    'ggc': 'G', 'gta': 'V', 'gca': 'A', 'gaa': 'E',
    'gga': 'G', 'gtg': 'V', 'gcg': 'A', 'gag': 'E',
    'ggg': 'G'
}
'''
###### Monte Carlo
# current random dna generator
def random_dna_sequence_v2(genome, length):
  if length%3!=0:
    raise ValueError('length needs to be disvisible by 3!')
  
  weights = []
  for i in CODE_COVARMAT:    
      mu, sigma = genome.mean_cub[i.upper()], genome.std_cub[i.upper()]
      s = np.random.normal(mu, sigma, 1)
      weights.append(s)
  list_of_cds = choices([cds.upper() for cds in CODE_COVARMAT.keys()], weights=weights, k=length)
  DNA=""
  for cds in list_of_cds:
    DNA+=cds
  return DNA

## generate sequences and create normed vector
def calculate_mean_cub_generated_normed(random_dnas):
    # calculate cub for each generated sequences
    list_of_cubs = []
    for i in range(len(random_dnas)):
        list_of_cubs.append(calc_cub(random_dnas[i]))
    # calculate cub means of 10 000 random sequences
    mean_cub_gen = {}
    for cds in CODE_COVARMAT:
        mean_cub_gen[cds.upper()] = 0
        
    for i in range(len(random_dnas)):
        curr_cub = list_of_cubs[i]
        for cds in CODE_COVARMAT:
            mean_cub_gen[cds.upper()] += curr_cub[cds.upper()]
    
    # norm
    total_codon_count = 0
    for cds in mean_cub_gen:
        total_codon_count+=mean_cub_gen[cds.upper()]
    
    for cds in mean_cub_gen:
        mean_cub_gen[cds] = mean_cub_gen[cds]*1000/total_codon_count
    
    return mean_cub_gen
######
'''



def calculate_mahalanobis_distance(covariance_matrix, gene, genome, total_cds_count_in_genome):
    X = np.zeros((1,61))
    Xhat = np.zeros((1,61))
    for i,cds in enumerate(CODE_COVARMAT):
        X[0][i] = genome[gene]['cub'][cds.upper()]
        Xhat[0][i] = (genome.cub[cds.upper()]/total_cds_count_in_genome)*1000
    
    # the X is normalised!
    substraction = (X*1000/np.sum(X))-Xhat
    inverse_covariance_matrix = np.linalg.inv(covariance_matrix)
    
    mahalonobis_distance = np.matmul(np.matmul(substraction,inverse_covariance_matrix),np.transpose(substraction))
    return mahalonobis_distance*1000
# calculate mahalanobis distances

def calculate_mahalanobis_distances(genome):
    ## MONTE CARLO PROCEDURE
    
    '''
    list_of_random_dnas = []
    for i in range(20000):
        list_of_random_dnas.append(random_dna_sequence_v2(genome, 300))
    
    mean_cub_gen_normed = calculate_mean_cub_generated_normed(list_of_random_dnas)
    '''
    
    
    ## CALCULATE COVARIANCE MATRIX
    # normed using cub + normed
    total_cds_count_in_genome = 0
    for locust_tag in genome.genes:
        for cds in CODE_COVARMAT:
            total_cds_count_in_genome += genome[locust_tag]['cub'][cds.upper()]
    covarmat_v1 = np.zeros((61,61))
    # fill covariance matrix
    for i,i_tag in enumerate(CODE_COVARMAT):
        for j,j_tag in enumerate(CODE_COVARMAT):
            intermediate_sum = 0
            for locust_tag in genome.genes:
                total_codon_count_locust_tag = 0
                
                for tag in genome[locust_tag]['cub']:
                    total_codon_count_locust_tag += genome[locust_tag]['cub'][tag]
                
                diff_A = ((genome[locust_tag]['cub'][i_tag.upper()]/total_codon_count_locust_tag)*1000) - ((genome.cub[i_tag.upper()]/total_cds_count_in_genome)*1000)
                diff_B = ((genome[locust_tag]['cub'][j_tag.upper()]/total_codon_count_locust_tag)*1000) - ((genome.cub[j_tag.upper()]/total_cds_count_in_genome)*1000)
                intermediate_sum+=diff_A * diff_B

            covarmat_v1[i][j]= intermediate_sum
            
    covarmat_v1

    ## CALCULATE DISTANCES
    for locust_tag in genome.genes:
        genome[locust_tag]['Mah'] = calculate_mahalanobis_distance(covarmat_v1, locust_tag,  genome, total_cds_count_in_genome)
        

def get_potential_HGT_Mah(genome):
    genome_mean_Mah = 0
    # Calculate mean mah distance
    for genes in genome.genes:
        genome_mean_Mah += genome[genes]['Mah']
    genome_mean_Mah = genome_mean_Mah/len(genome)
    
    # calculate sum difference for std
    sum_diff = 0
    for genes in genome.genes:
        diff = genome[genes]['Mah'] - genome_mean_Mah
        diff = diff*diff
        sum_diff += diff
        
    # calc std
    genome_std_Mah = math.sqrt(sum_diff/len(genome))
    
    # extraneous is defined as:
    # whe difference from mean is 2 times more than sigma
    list_of_extraneous_Mah = []
    for genes in genome.genes:
        diff = genome[genes]['Mah'] - genome_mean_Mah 
        if diff > 1.5 * genome_std_Mah:
            if diff > 2 * genome_std_Mah:
                # odd / extraneous
                list_of_extraneous_Mah.append(genes)
                genome[genes]['SimMah'] = 2
            else:
                # potential odd
                genome[genes]['SimMah'] = 1
        else:
            genome[genes]['SimMah'] = 0
                
    
    return list_of_extraneous_Mah