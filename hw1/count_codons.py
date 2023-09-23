#--------------
# Imports
#--------------
# !pip install -r requirements.txt
import fastaparser
from collections import Counter # dictionary subclass
import pandas as pd
import itertools
import argparse

#--------------
# Functions
#--------------

def make_argparser():
    '''
    Parameters
    -----------
        Technically, this function takes no parameters. Note, though, that the parser expects two strings (separated by a space) after the command count_codons.py: the input_filename, and the output_filename. See the count_codons() docstring for more info.

    Return
    --------
        This function returns the an argparse parser, which is then used as the input for the count_codons() function. See the docs for more info on how the argparse methods work: https://docs.python.org/3/library/argparse.html
    '''
    parser = argparse.ArgumentParser(prog = 'count_codons.py', 
                                     description = 'A program to count 3-base codons from an input FASTA (.fna) file')
    
    # Path to input file
    parser.add_argument('input_filename',
                        help = 'Path to input fasta file')
    
    # Name of output file
    parser.add_argument('output_filename',
                        help = 'Desired output file name, with .csv extension included')
    return parser

def count_codons(input_filename: str, output_filename: str):
    '''
    Parameters
    -----------
        input_filename: str 
            The name of the input fasta file to be used. The '.fna' extension must be included in the string.
            Example: 'SARS-CoV-2_separate_genes.fna'

        output_filename: str 
            The desired name for the output csv file. The '.csv' extension must be included in the string.
            Example: 'separate_genes_output.csv'

    Return
    --------
        This function writes the output csv file to the root folder; it returns nothing.
    '''
    # Read in FASTA file
    # https://pypi.org/project/fastaparser/
    with open(input_filename) as fasta_file:
            parser = fastaparser.Reader(fasta_file)
            seqs = [seq.sequence_as_string() for seq in parser]

    # For each length of sequence, split into codons, starting at first position (assumes no introns, no promoter region, etc.)
    # Note that this works for sequences that are separated into their component genes as well as whole-genome sequences
    codons_ls = list(map(lambda x: [seqs[x][i:i+3] for i in range(0, len(seqs[x]), 3)], range(len(seqs))))
    codons = list(itertools.chain.from_iterable(codons_ls))

    # Count codons
    codon_dict = Counter(codons)

    # Convert to dataframe which can be written to csv & rename cols
    codon_df = pd.DataFrame.from_dict(codon_dict, orient = 'index')
    codon_df = codon_df.reset_index()
    codon_df.columns = ['codon', 'count']
    codons_full = codon_df.query('codon.str.len() == 3') # drop rows with incomplete sequence (not having 3 bases)
    codons_full.to_csv(output_filename, index = False)

if __name__ == '__main__':
    parser = make_argparser()
    args = parser.parse_args()
    count_codons(args.input_filename, args.output_filename)