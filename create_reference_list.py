from Bio import Entrez
import pandas as pd
import time
import os



def get_assembly_summary(id):
    """Get esummary for an entrez id"""
    esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
    esummary_record = Entrez.read(esummary_handle)
    return esummary_record

def check_genome(name, email, api_key):
    #def get_assemblies(term, download=True, path='sequence_files'):
    """Download genbank assemblies for a given search term.
    https://dmnfarrell.github.io/bioinformatics/assemblies-genbank-python
    Args:
        term: search term, usually organism name
        download: whether to download the results
        path: folder to save to
    """

    #provide your own mail here
    Entrez.email = email
    Entrez.api_key = api_key
    handle = Entrez.esearch(db="assembly", term=name, retmax='200')
    record = Entrez.read(handle)
    ids = record['IdList']
    print (f'found {len(ids)} ids')
    links = []
    genbank_accession = []
    refseq_accession = []
    for id in ids:
        #get summary
        summary = get_assembly_summary(id)
        #print(summary)
        #get ftp link
        # for now going forward no refseq files are going to be used!
        url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
        gb_accession =summary['DocumentSummarySet']['DocumentSummary'][0]['Synonym']['Genbank']
        rs_accession =summary['DocumentSummarySet']['DocumentSummary'][0]['Synonym']['RefSeq']
        #print(summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq'])
        if url == '':
            continue
        label = os.path.basename(url)
        #get the fasta link - change this to get other formats
        link = os.path.join(url,label+'_cds_from_genomic.fna.gz')
        links.append(link)
        genbank_accession.append(gb_accession)
        refseq_accession.append(rs_accession)

    if len(links)>=1:
        return links[0], genbank_accession[0], refseq_accession[0]
    elif len(links) == 0:
        return '', '', ''
        
        
def string_separator(string_test):
    comma_index = string_test.find(',')
    if comma_index != -1:
        result = string_test[:comma_index]
    else:
        result = string_test
    return result
        
def main():
    hgtree_genomes = pd.read_csv('data/HGTREE/version_1/Genome_Information_20150618.txt', sep='\t')
    gban = hgtree_genomes['GenBank accession number'].tolist()
    hgtree_genomes['gban'] = [string_separator(acession) for acession in gban]
    
    spam_counter =0
    links, gbs, rss = [],[],[]
    for i in gban:
        try:
            link, gb, rs = check_genome(i, "adjie.salman@stud.uni-due.de", "726c2709d1827c981a38403d4a7d99e5cf08")
        except:
            link, gb, rs = 'fetch_error','fetch_error','fetch_error'
        links.append(link)
        gbs.append(gb)
        rss.append(rs)
        spam_counter +=1
        if spam_counter ==3:
            print('sleeping')
            time.sleep(3)
            spam_counter = 0
            print('continuing')
            
        
    hgtree_genomes['link_to_file'] = links
    hgtree_genomes['genbank_accession'] = gbs
    hgtree_genomes['refseq_acession'] = rss
    
    hgtree_genomes.to_csv('complete_reference_redo.csv')
        
        
    
        
if __name__ == '__main__':
    main()