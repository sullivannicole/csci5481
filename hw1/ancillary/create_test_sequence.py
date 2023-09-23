import fastaparser
import pandas as pd

# Create test FASTA file
# Generated random 20-base seq here: https://www.bioinformatics.org/sms2/random_dna.html
with open("fake_genome.fna", 'w') as fasta_file:
        writer = fastaparser.Writer(fasta_file)
        writer.writefasta(('id20 test sequence', 'CCGCACTTTGTCAGAGTCGG'))


# Expected output
expected_output = pd.DataFrame({'codon': ['CCG', 'CAC', 'TTT', 'GTC', 'AGA'],
                                'count': [1, 1, 1, 2, 1]})

expected_output.to_csv('fake_genome_expected_output.csv', index = False)