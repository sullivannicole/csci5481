# Databricks notebook source
# ------------
# Imports
# ------------

# !pip install fastaparser
import fastaparser
import pandas as pd
import numpy as np
import re
import itertools   

def parse_header(header: list) -> dict:
    """
    Parameters
    -----------
        header: a list produced of sequence headers produced by the fastaparser library reading a .fna file

    Return
    --------
        A dictionary containing a key-value pair for each of the items in the header list
    """

    header_split = re.findall("\[(.*?)\]", header)
    header_ls = [{x.split("=")[0]: x.split("=")[1]} for x in header_split]

    return {k: v for d in header_ls for k, v in d.items()}


def read_seqs(fna_file: str, 
              seq_only: bool = True) -> pd.DataFrame:
    """
    Parameters
    -----------
    fna_file: str
        The path to the fna file to read in, including the '.fna' extension. E.g. 'pfizer_mrna.fna' for a local fasta file containing Pfizer mRNA sequence.
    seq_only: Boolean
        Whether or not to return just the sequences contained in the file, or to include all the metadata contained in the fasta file (headers, tags, locii, etc.). Set to True by default, since it's often useful to extract just the sequence.

    Return
    --------
    pd.DataFrame 
        A dataframe with one row per sequence in the input fasta file. At the very least (if seq_only = True), the dataframe will contain the columns "sequence", "lcl" (the sequence headers), and "specimen", which is the file name, excluding the ".fna" extension. If seq_only = False, the output dataframe will contain all additional columns that can be extracted from the input fasta file. These additional columns vary depending on what's contained in the input, but can include, tags, locii of the start/end of the sequence within the entire genome, and other useful metadata.
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

def calc_max_conserved(seq_matrix_col):
    return np.max(list(map(lambda x: np.sum(np.char.count(x, seq_matrix_col)), ['A', 'T', 'G', 'C'])))

def sum_bases(seq_matrix_col):
    return np.sum(list(map(lambda x: np.sum(np.char.count(x, seq_matrix_col)), ['A', 'T', 'G', 'C'])))

# COMMAND ----------

# MAGIC %md
# MAGIC
# MAGIC # Problem 1

# COMMAND ----------

# Create numpy array with one base per cell, one gene per row
genes = read_seqs('inputs/Homework4-seqs-with-primers.fna', seq_only = True)
seq_array = [np.array(list(x), dtype = str) for x in genes]
seq_mat = np.vstack(seq_array)

# Calculate % conserved at each position
conserved_counts = list(map(calc_max_conserved, seq_mat.T)) # transpose so it iterates over columns (which will be rows)
# base_sums = list(map(sum_bases, seq_mat.T))
# pc_conserved = np.divide(np.array(conserved_counts), np.array(base_sums)) # i thought gaps **weren't** supposed to count at first
pc_conserved = np.array(conserved_counts)/len(genes)

# Create dataframe and write out
conservation_df = pd.DataFrame(pc_conserved, columns = ['pc_conserved']).reset_index(inplace = False)
conservation_df.to_csv('outputs/solution-problem-1.txt', sep = '\t', index = False)
spark.createDataFrame(conservation_df).createOrReplaceTempView('vw_pc_conserved')

# COMMAND ----------

# MAGIC %md
# MAGIC
# MAGIC # Problem 2

# COMMAND ----------

# MAGIC %r
# MAGIC library(SparkR)
# MAGIC library(tidyverse)
# MAGIC library(svglite)
# MAGIC
# MAGIC options(repr.plot.width = 1500, repr.plot.height = 700)
# MAGIC
# MAGIC lt_blue <- "#96C1FF"
# MAGIC neon_blue <- "#0061FF"
# MAGIC bgrd <- "black"
# MAGIC txt_color <- "white"
# MAGIC
# MAGIC theme_bar_narrow <- theme(panel.background = element_rect(fill = bgrd, color = bgrd),
# MAGIC                           plot.background = element_rect(fill = bgrd, color = bgrd),
# MAGIC                           panel.grid.major = element_blank(),
# MAGIC                           panel.grid.minor = element_blank(),
# MAGIC                             axis.ticks = element_blank(),
# MAGIC                             plot.title = element_text(color = txt_color, face = "bold", size = 25),
# MAGIC                             axis.title = element_text(color = txt_color, size = 20),
# MAGIC                             axis.text = element_text(size = 13, color = txt_color),
# MAGIC                             legend.position = "bottom", 
# MAGIC                             legend.title = element_text(color = txt_color, size = 15), 
# MAGIC                             legend.text = element_text(color = "white", size = 13),
# MAGIC                             legend.key = element_rect(fill = "transparent"),
# MAGIC                             legend.background = element_rect(fill="transparent"),
# MAGIC                             plot.margin = margin(2.2,2.2,2.2,2.2, "cm"))
# MAGIC
# MAGIC %r
# MAGIC # Visualize rolling averages - previous 5 dicodons and next 5 dicodons
# MAGIC pc_conserved <- SparkR::collect(SparkR::sql('select index, pc_conserved, avg(pc_conserved) over (partition by "none" order by index rows between 30 preceding and 30 following) as roll_avg
# MAGIC from vw_pc_conserved;'))
# MAGIC
# MAGIC pc_conserved %>%
# MAGIC ggplot(aes(index, roll_avg)) +
# MAGIC geom_line(color = neon_blue, size = 0.6) +
# MAGIC labs(x = '\nPosition in 16S rRNA gene',
# MAGIC y = '% conserved\n') +
# MAGIC theme_bar_narrow
# MAGIC
# MAGIC ggsave("outputs/solution-problem-2.svg")

# COMMAND ----------

# %r

# pc_conserved %>%
# mutate(lo_hi = ifelse(roll_avg < 0.70, 'lo', 'hi')) %>%
# ggplot(aes(index, roll_avg, color = lo_hi, group = 1)) +
# geom_line(size = 0.6) +
# scale_color_manual(values = c(neon_blue, "#00BAB6")) +
# labs(x = '\nPosition in 16S rRNA gene',
# y = '% conserved\n',
# color = '') +
# theme_bar_narrow

# COMMAND ----------

# MAGIC %md
# MAGIC
# MAGIC # Problem 3

# COMMAND ----------

# MAGIC %sql
# MAGIC
# MAGIC -- Use 0.71 as threshold for variable/conserved regions
# MAGIC create or replace temporary view q3 as
# MAGIC
# MAGIC with mov_av as (
# MAGIC   select *, avg(pc_conserved) over (partition by "none" order by index rows between 30 preceding and 30 following) as roll_avg
# MAGIC   from vw_pc_conserved
# MAGIC ),
# MAGIC
# MAGIC var_cat as (select *, if(roll_avg > 0.71, 1, 0) as var_region
# MAGIC from mov_av),
# MAGIC
# MAGIC grouped_vals as (select *, lag(var_region) over (partition by 'nothing' order by index) as grouped_var
# MAGIC from var_cat),
# MAGIC
# MAGIC nums_to_sum as (select *, if(var_region != grouped_var, 1, 0) as num_to_sum
# MAGIC from grouped_vals
# MAGIC ),
# MAGIC
# MAGIC var_reg as (select *, sum(num_to_sum) over (partition by "nothing" order by index) as var_region_grouper
# MAGIC from nums_to_sum),
# MAGIC
# MAGIC regions_summarized as (select var_region_grouper, min(index) as start_position, max(index) as end_position,
# MAGIC if(var_region_grouper % 2 = 0 or var_region_grouper = 0, 'variable', 'conserved') as region_type
# MAGIC from var_reg
# MAGIC group by 1)
# MAGIC
# MAGIC select * except(var_region_grouper)
# MAGIC from regions_summarized
# MAGIC order by start_position;

# COMMAND ----------

q3_df = spark.sql('select * from q3').toPandas()
q3_df.to_csv('outputs/solution-problem-3.txt', sep = '\t', index = False)

# COMMAND ----------

# MAGIC %md
# MAGIC
# MAGIC # Problem 4

# COMMAND ----------

# MAGIC %r
# MAGIC pc_conserved %>%
# MAGIC mutate(lo_hi = ifelse(roll_avg < 0.70, 'lo', 'hi'),
# MAGIC color_values = ifelse(roll_avg < 0.7, 0.7, roll_avg)) %>%
# MAGIC ggplot() +
# MAGIC geom_line(aes(index, roll_avg, color = color_values, group = 1), size = 0.6) +
# MAGIC scale_color_gradient2(mid = neon_blue, low = "#00BAB6", midpoint = 0.8) + # mid =
# MAGIC labs(x = '\nPosition in 16S rRNA gene',
# MAGIC y = '% conserved\n') +
# MAGIC theme_bar_narrow +
# MAGIC theme(legend.position = "none")
# MAGIC
# MAGIC ggsave("outputs/solution-problem04.svg")
