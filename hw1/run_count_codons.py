# Databricks notebook source
# MAGIC %md
# MAGIC
# MAGIC This script runs the shell commands for each of the 3 input FASTA files.

# COMMAND ----------

# MAGIC %sh python count_codons.py 'SARS-CoV-2_whole_genome.fna' 'whole_genome_output.csv'

# COMMAND ----------

# MAGIC %sh python count_codons.py 'SARS-CoV-2_separate_genes.fna' 'separate_genes_output.csv'

# COMMAND ----------

# MAGIC %sh python count_codons.py 'fake_genome.fna' 'fake_genome_actual_output.csv'
