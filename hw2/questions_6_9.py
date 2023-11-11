# Databricks notebook source
# MAGIC %md
# MAGIC
# MAGIC # Question 6

# COMMAND ----------

import pandas as pd 

q3 = pd.read_csv("question_3_output.txt")

# COMMAND ----------

# Excl start/end gaps - look only at mismatches in coding portions
coding_region = q3['1711'][2][54:-299]

# COMMAND ----------

# Get a count of matches, mismatches and gaps
coding_region.count('|')
coding_region.count('x')
coding_region.count(' ')

# COMMAND ----------

# MAGIC %md
# MAGIC
# MAGIC # Question 9

# COMMAND ----------

# Get a count of mismatches + gaps
q8 = pd.read_csv("question_8_output.txt")
q8['1270'][2].count('x') + q8['1270'][2].count(' ')

# COMMAND ----------

q8['1270'][1][985:987] +\
q8['1270'][2][985:987] +\
q8['1270'][3][985:987]

# COMMAND ----------

# MAGIC %md
# MAGIC
# MAGIC # Question 11

# COMMAND ----------

q3['1711'][1].count('G') + q3['1711'][1].count('C') # 1426
q3['1711'][3].count('G') + q3['1711'][3].count('C') # 2373