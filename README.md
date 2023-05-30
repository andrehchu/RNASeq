# RNASeq
This is our CSE 185 Project that performs differential gene expression analysis on RNA-seq data. Our implementation includes processing and normalizing the RNA-seq count data and statistical testing using the Mann Whitney U test to identify any significant genes with significantly different levels of expression. This tool also generates volcano plots and MA-plots for the resulting data. We modeled this tool around the DESeq2 R package for differential expression analysis using python rather than R.

# Installation instructions

# Basic usage instructions
Our code is currently just a python script and has not been made into a package yet.
The basic usage of [207analysis.py] is as follows:

``` python3 207analysis.py [-c countdata.txt] [-o outputfile.csv] [-padj float] [-pval_thresh float ```

Small test example command:

``` python3 207analysis.py -c GSE221626_counts.txt -o -padj bonferroni -pval_thresh 0.5 ```

# Complete usage instructions
There are four required inputs for [analysis]

* ```-c FILE ```, ```--countdata FILE ```: txt file that contains RNA-seq count data
* ```-o FILE, --output_file```: csv files that contains the output data
* ```-padj PADJ```, ```--pvalue_adjusted PADJ```: options bonferroni or fdr_bh for the padjusted values
* ```-pval_thresh THRESH```, ```--pvalue_threshold THRESH```: decimal value that specifies the pvalue threshold for volcano and MA-plots

# Credits

