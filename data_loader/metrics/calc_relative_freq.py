'''
Calculate Relative Frequency of a sequence

Relative Frequency is described as fk(i) = ck(i)/lk
ck(i) is the ocunt of nucleotide i in sequence sk
lk is length of sequence

Definition from:
Detecting laterally transferred genes: use of entropic clustering methods and genome position
(Azad and Lawrence, Nucleic Acids Research, 2007)
https://doi.org/10.1093/nar/gkm204

'''

def calc_relative_freq(sequence):
    g_count=0 
    a_count=0 
    c_count=0 
    t_count = 0
    for base in sequence:
            if base == 'G':
                g_count +=1
            elif base == 'A':
                a_count +=1
            elif base == 'C':
                c_count +=1
            elif base == 'T':
                t_count +=1

    rep = {
       'G': g_count/len(sequence),
        'A': a_count/len(sequence),
        'C': c_count/len(sequence),
        'T': t_count/len(sequence)
    }
    return rep