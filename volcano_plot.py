##Volcano Plot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

deseqresults = pd.read_csv("~/lab4-spring23/chow_vs_hfd_deseq2.csv")
plt.figure(figsize=(10,10))

pvalue_threshold = 0.05
fc_threshold = 0

deseqresults = deseqresults.rename({"Unnamed: 0" : "gene_id"}, axis='columns')


plt.xlabel('log2FoldChange') #x label
plt.ylabel('-log10pvalue') #y label
plt.title('Differential Expression')

#plot scatter plot while labeling significant expression with red based on padj<0.05 AND log2fold > 0 or log2fold < 0
plt.scatter(deseqresults['log2FoldChange'], -np.log10(deseqresults['padj']), c = np.where((deseqresults["padj"] < 0.05)  & ((deseqresults["log2FoldChange"] < 0) | (deseqresults["log2FoldChange"] > 0)), "red", "black"))


#using threshold of p value 5% and taking top 10
deseqresults = deseqresults[deseqresults['padj'] < 0.05]
pSort =  deseqresults.sort_values('padj')
top_ten = pSort.iloc[:10]

for index, row in top_ten.iterrows():
    plt.annotate(row["gene_id"], (row["log2FoldChange"], -np.log10(row["pvalue"])))
    
plt.savefig("VolcanoPlotDeseq2.png")
    
