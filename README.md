# Introduction
pySeqDiff is a Python package designed to perform differential gene expression analysis on raw RNA-seq count data. The aim of this project was to replicate the method used in the R package, DESeq2, in order to better understand how differential gene expression analysis is performed, and to allow this kind of analyis to be performed using Python. Our tool has a similar output file format as DESeq2, and also generates a volcano plot and MA plot to the output directory.

We appreciate any feedback. This tool was developed by students at UCSD as a part of a project for CSE 185: Advanced Bioinformatics Laboratory in Spring 2023. 

# Usage Instructions:

## Setup:
``` git clone https://github.com/andrehchu/RNASeq ```

``` cd RNASeq ```

``` pip install . ```

Make sure the necessary packages are downloaded before running the program!


Run the following command to install/upgrade the necessary packages:

``` pip install -r prereqs.txt ```

## Basic Usage Instructions:

``` pySeqDiff -c <raw_countdata.txt> -o <output_file.txt> [-pval_thresh <significance threshold>] [-s] ```

Type ``` pySeqDiff --help ``` for more usage information and description of arguments.

## Example test:

Use the count data file provided to run the following command:

``` pySeqDiff -c pySeqDiff/GSE207721_kyse150_raw.txt -o pySeqDiff_kyse150_results.txt -pval_thresh 0.05 -s ```

You may encounter a warning message as such:

/Library/Frameworks/Python.framework/Versions/3.11/lib/python3.11/site-packages/pandas/core/arraylike.py:396: RuntimeWarning: divide by zero encountered 
in log10
  result = getattr(ufunc, method)(*inputs, **kwargs)

This can be ignored and will not affect the output file or plots.

# How to interpret the output:
 The count data file used for testing is obtained from this study 'https://pubmed.ncbi.nlm.nih.gov/35932580/'. The output shows the gene ID, if it is experimental or wildtype, its log 2 fold change, p-value and adjusted p-value, and the corresponding base mean. The data can be found at GEO Accession Viewer ID GSE207721, where we concatenated data from the kyse150 cell line raw read counts (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi).

# References

Fatemeh Seyednasrollah et. al. “Comparison of software packages for detecting differential expression in RNA-seq studies”, Briefings in Bioinformatics, Volume 16, Issue 1, January 2015, Pages 59–70, https://doi.org/10.1093/bib/bbt086.

Khetani, Radhika. “Count Normalization with DEDeq2.” Introduction to DGE - ARCHIVED, 26 Apr. 2017, hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html. 
Love, M.I. et. al. “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2”. Genome Biol 15, 550 (2014). https://doi.org/10.1186/s13059-014-0550-8.

Wang, Jian et. al. “KDM2A plays a dual role in regulating the expression of malignancy-related genes in esophageal squamous cell carcinoma”, Biochemical and Biophysical Research Communications, Volume 624, 2022, Pages 53-58, ISSN 0006-291X, https://doi.org/10.1016/j.bbrc.2022.07.035.


# Acknowledgements
Thank you to professor Melissa Gymrek, teaching assistants Luisa Amaral, and Ryan Eveloff, as well as peers for guidance and support.
