# Databricks notebook source
# MAGIC %md
# MAGIC
# MAGIC This script runs the shell commands for the given example, as well as for tasks 2, 3 and 8.
# MAGIC It takes about ~2min to run the full covid-Pfizer alignment on a 20-core machine.

# COMMAND ----------

# MAGIC %sh python aligner.py -q 'data/pfizer_mrna_truncated_example.fna' -r 'data/sars_spike_protein_truncated_example.fna' -o 'outputs/output_trunc.txt' -g -2 -p -1 -m 1

# COMMAND ----------

# MAGIC %sh python aligner.py -q 'data/pfizer_mrna_truncated_example.fna' -r 'data/sars_spike_protein_truncated_example.fna' -o 'outputs/output_trunc_ignore_outer.txt' -g -2 -p -1 -m 1 --ignore_outer_gaps

# COMMAND ----------

# MAGIC %sh python aligner.py -q 'data/pfizer_mrna.fna' -r 'data/sars_spike_protein.fna' -o 'outputs/question_2_output_test2.txt' -g -2 -p -1 -m 1

# COMMAND ----------

# MAGIC %sh python aligner.py -q 'data/pfizer_mrna.fna' -r 'sars_spike_protein.fna' -o 'outputs/question_3_output.txt' -g -2 -p -1 -m 1 --ignore_outer_gaps

# COMMAND ----------

# MAGIC %sh python aligner.py -q 'pfizer_mrna.aa' -r 'sars_spike_protein.aa' -o 'outputs/question_8_output.txt' -g -2 -p -1 -m 1