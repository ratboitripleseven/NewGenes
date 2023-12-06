'''
Calculate 48 symbol representation of a sequence

48 symbol is described as dinucleotide identity and the three codon positions

Definition from:
Detecting laterally transferred genes: use of entropic clustering methods and genome position
(Azad and Lawrence, Nucleic Acids Research, 2007)
https://doi.org/10.1093/nar/gkm204

'''


def calc_48_symbols(sequence):
    dn = {}
    for nt in ["AA", "AG", "AT", "AC", "TA", "TG", "TT", "TC","GA", "GG", "GT", "GC","CA", "CG", "CT", "CC"]:
        dn[nt] = [0, 0, 0]
    for i in range(0, len(sequence), 3):
        codon = sequence[i : i + 4]
        if len(codon) < 4:
            codon += "  "
        for pos in range(0, 3):
            for nt in ["AA", "AG", "AT", "AC", "TA", "TG", "TT", "TC","GA", "GG", "GT", "GC","CA", "CG", "CT", "CC"]:
                if codon[pos] == nt[0] and codon[pos+1] == nt[1] :
                    dn[nt][pos] += 1
                elif codon[pos] == nt[0].lower() and codon[pos+1] == nt[1].lower():
                    dn[nt][pos] += 1
                    
    # calculate length and divide
    for codon_pos in range(0,3):
        total_counts_per_pos = 0
        for nt in ["AA", "AG", "AT", "AC", "TA", "TG", "TT", "TC","GA", "GG", "GT", "GC","CA", "CG", "CT", "CC"]:
            total_counts_per_pos += dn[nt][codon_pos]
        for nt in ["AA", "AG", "AT", "AC", "TA", "TG", "TT", "TC","GA", "GG", "GT", "GC","CA", "CG", "CT", "CC"]:
            dn[nt][codon_pos] =dn[nt][codon_pos]/total_counts_per_pos
    return dn