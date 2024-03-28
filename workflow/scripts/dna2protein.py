from Bio import SeqIO
from Bio.Seq import Seq
#from Bio.Alphabet import IUPAC
import sys

sys.stderr = open(snakemake.log[0], "w")

input_file = snakemake.input.fasta[0]
output_file = snakemake.output.faa


def translate_dna_to_protein(input_file, output_file):

    # Open the input FASTA file
    with open(input_file, "r") as handle:
        # Open output file to write protein sequences
        with open(output_file, "w") as output_handle:
            # Parse each sequence in the FASTA file
            for record in SeqIO.parse(handle, "fasta"):#, alphabet=generic_dna):
                # Translate DNA sequence to protein sequence
                protein_sequence = record.seq.translate()
                # Write protein sequence to the output file
                output_handle.write(f">{record.id}\n{protein_sequence}\n")


translate_dna_to_protein(input_file, output_file)