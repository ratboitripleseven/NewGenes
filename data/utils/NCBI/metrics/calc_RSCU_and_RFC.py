'''
Calculate RSCU and RFC together
RSCU
http://genomes.urv.es/CAIcal/tutorial.pdf
https://www.biostars.org/p/437939/

RFC
https://www.biotite-python.org/examples/gallery/sequence/codon_usage.html

RSCU is defined as S × Nc/Na, where 
S represents the number of synonymous codons encoding the same amino acid,
Nc is the frequency of the codon in the genome, and 
Na is the relative frequency of the codon for that amino acid

RFC is defined as Nc/Na

NOTE: This script requires codon usage bias!


'''
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

def calc_RSCU_and_RFC(gene_cub):
    RSCU = {}
    RFC = {} 
    for cds in CODE:
        RSCU[cds.upper()]=0
        RFC[cds.upper()]=0
    for cds in CODE:
        # check amino acid of current codon
        amino_acid = CODE[cds]
        # get total amount of synonymous codon
        S = list(CODE.values()).count(amino_acid) 
        # get total count of synonymous codn
        Na = 0
        for codon in CODE:
            if CODE[codon]==amino_acid:
                Na+= gene_cub[codon.upper()]
        
        Nc = gene_cub[cds.upper()]
        
        if Na == 0:
            rfc = 0
        else:
            rfc = (Nc / Na)
        
        rscu = S * rfc
        
        RSCU[cds.upper()] = rscu
        RFC[cds.upper()] = rfc
    # return rscu and rf
    return RSCU, RFC
    
    