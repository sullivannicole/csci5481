# Databricks notebook source
import pandas as pd
my_edges = pd.read_csv("example2/edges.txt", names = ['anc', 'desc', 'length'], sep = "\t")
canvas_edges = pd.read_csv("example2/edges_canvas.txt", names = ['anc', 'desc', 'length'], sep = "\t")

# COMMAND ----------

# 1-> 4
my_edges.query("anc == 11 and desc == 1").length.values[0] +\
    my_edges.query("anc == 11 and desc == 12").length.values[0] +\
        my_edges.query("anc == 12 and desc == 4").length.values[0]

# COMMAND ----------

# 1-> 4
canvas_edges.query("anc == 11 and desc == 1").length.values[0] +\
    canvas_edges.query("anc == 11 and desc == 12").length.values[0] +\
        canvas_edges.query("anc == 12 and desc == 4").length.values[0]

# COMMAND ----------

# 1-> 6
my_edges.query("anc == 8 and desc == 6").length.values[0] +\
    my_edges.query("anc == 8 and desc == 9").length.values[0] +\
        my_edges.query("anc == 9 and desc == 11").length.values[0] +\
             my_edges.query("anc == 11 and desc == 1").length.values[0]

# COMMAND ----------

# 1-> 6
canvas_edges.query("anc == 9 and desc == 6").length.values[0] +\
    canvas_edges.query("anc == 8 and desc == 9").length.values[0] +\
        canvas_edges.query("anc == 8 and desc == 11").length.values[0] +\
             canvas_edges.query("anc == 11 and desc == 1").length.values[0]

# COMMAND ----------

# 4 -> 7
my_edges.query("anc == 8 and desc == 6").length.values[0] +\
    my_edges.query("anc == 8 and desc == 9").length.values[0] +\
        my_edges.query("anc == 9 and desc == 11").length.values[0] +\
             my_edges.query("anc == 11 and desc == 1").length.values[0]
