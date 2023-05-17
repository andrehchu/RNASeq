import pandas as pd
import numpy as np
from scipy import stats
import statsmodels.stats.multitest as smm

count_data = {}

with open("data.txt") as f:
    header = f.readline().strip().split("\t")
    for line in f:
        fields = line.strip().split("\t")
        gene_id = fields[0]
        counts = [float(x) if x != "NaN" else x for x in fields[1:]]
        if "NaN" in counts:
            print(f"Found 'NaN' value in gene_id: {gene_id}")
        count_data[gene_id] = counts