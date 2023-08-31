'''
Calculate 12 symbol representation of a sequence

12 symbol is described as nucleotide identity and the three codon positions

Definition from:
Detecting laterally transferred genes: use of entropic clustering methods and genome position
(Azad and Lawrence, Nucleic Acids Research, 2007)
https://doi.org/10.1093/nar/gkm204

'''

def calc_12_symbols(sequence):
    d = {}
    for nt in ["A", "T", "G", "C"]:
        d[nt] = [0, 0, 0]
    for i in range(0, len(sequence), 3):
        codon = sequence[i : i + 3]
        if len(codon) < 3:
            codon += "  "
        for pos in range(0, 3):
            for nt in ["A", "T", "G", "C"]:
                if codon[pos] == nt or codon[pos] == nt.lower():
                    d[nt][pos] += 1

    #print(d)
    # calculate length and divide
    for codon_pos in range(0,3):
        total_counts_per_pos = 0
        for nt in ["A", "T", "G", "C"]:
            total_counts_per_pos += d[nt][codon_pos]
        for nt in ["A", "T", "G", "C"]:
            d[nt][codon_pos] =d[nt][codon_pos]/total_counts_per_pos
            
    return d