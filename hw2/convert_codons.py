# Databricks notebook source
!pip install fastaparser
import pandas as pd
import fastaparser
import re

# COMMAND ----------

def parse_header(header: list) -> dict:
    """
    Parameters
    -----------
        header:

    Return
    --------
    """

    header_split = re.findall("\[(.*?)\]", header)
    header_ls = [{x.split("=")[0]: x.split("=")[1]} for x in header_split]

    return {k: v for d in header_ls for k, v in d.items()}


def read_seqs(fna_file: str, 
              seq_only: bool = True) -> pd.DataFrame:
    """ 
    """
    with open(fna_file) as fasta_file:
        rparser = fastaparser.Reader(fasta_file)
        seqs = [seq.sequence_as_string() for seq in rparser]

        if seq_only:
            return seqs

        qparser = fastaparser.Reader(fasta_file, parse_method="quick")
        headers = [seq.header for seq in qparser]

    seq_df = pd.DataFrame(
        {
            "sequence": seqs,
            "lcl": [x.split("[")[0] for x in headers],
            "specimen": fna_file[:-4],
        }
    )

    headers_df = pd.DataFrame(list(map(parse_header, headers)))
    genome_df = pd.concat([seq_df, headers_df], axis=1)

    return genome_df

# COMMAND ----------

# Read in files
aa_xwalk = pd.read_csv('data/codon_aa_xwalk.csv')
spike = read_seqs('data/sars_spike_protein.fna', seq_only = False)
pfizer = read_seqs('data/pfizer_mrna.fna', seq_only = False)

# COMMAND ----------

# Split into codons
def split_into_codons(seq):
    codon_ls = [seq[i:i+3] for i in range(0, len(seq), 3)]
    return codon_ls

spike_codons = split_into_codons(spike['sequence'][0])
pfizer_codons = split_into_codons(pfizer['sequence'][0][54:-299]) # Exclude start/end gaps

# COMMAND ----------

# Create dfs and connect to aa crosswalk
spike_df = pd.DataFrame({'V1': spike_codons})
pfizer_df = pd.DataFrame({'V1': pfizer_codons})

spike_aas = pd.merge(spike_df, aa_xwalk, on = 'V1', how = 'left')
pfizer_aas = pd.merge(pfizer_df, aa_xwalk, on = 'V1', how = 'left')

# Convert aa's to one long string
spike_aa_str = spike['lcl'][0] + '\n' + ''.join(list(spike_aas.V3.values))
pfizer_aa_str = pfizer['lcl'][0] + '\n' + ''.join(list(pfizer_aas.V3.values))

# Write out results
with open('outputs/sars_spike_protein.aa', 'w') as f:
    f.write(spike_aa_str)

with open('outputs/pfizer_mrna.aa', 'w') as f:
    f.write(pfizer_aa_str)

# COMMAND ----------

spike_2 = pd.read_csv('sars_spike_protein.aa')