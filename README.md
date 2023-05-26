# RNASeq
This is our CSE 185 Project that performs differential gene expression analysis on RNA-seq data. Our implementation includes processing and normalizing the RNA-seq count data and statistical testing using the Mann Whitney U test to identify any significant genes with significantly different levels of expression. This tool also generates volcano plots and MA-plots for the resulting data. We modeled this tool around the DESeq2 R package for differential expression analysis using python rather than R.

#Installation instrutions

#Basic usage instructions
The basic usage of [analysis.py] is as follows:

analysis [-c countdata.txt] [-o outputfile.txt]

Small test example command:

analysis -c GSE221626_counts.txt -o -padj bonferroni -pval_thresh 0.5

#Complete usage insturctions

#Credits

#Bonus:badges