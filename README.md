# Basic Usage Instructions:

# Setup:
``` https://github.com/andrehchu/RNASeq ```

``` cd RNASeq ```

# Usage Instructions:
``` python3 207analysis.py -c [RNA-seq count data file] -o [output_file] -s [sorted output if present]  ```
# Use the count data file I have provided and run the command:
``` python3 207analysis.py -c GSE207721_raw_counts_GRCh38.p13_NCBI.tsv -o [output file name of your choice] -s ```

# How to interpret the output:
 The count data file used for testing is obtained from this study 'https://pubmed.ncbi.nlm.nih.gov/35932580/'. The output shows the gene ID, if it is experimental or wildtype, its log 2 fold change, p-value and adjusted p-value, and the corresponding base mean.