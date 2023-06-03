# Basic Usage Instructions:

# Setup:
``` git clone https://github.com/andrehchu/RNASeq ```

``` pip install . ```

# Usage Instructions:
Make sure the necessary packages are downloaded before running the program!

# If necessary: use the following commands:

Run the following command to install/upgrade the necessary packages:
``` pip install -r prereqs.txt ```

# Use the count data file I have provided and run the command:
``` pySeqDiff -c pySeqDiff/GSE207721_raw_counts_GRCh38.p13_NCBI.tsv -o out.txt -s ```

# How to interpret the output:
 The count data file used for testing is obtained from this study 'https://pubmed.ncbi.nlm.nih.gov/35932580/'. The output shows the gene ID, if it is experimental or wildtype, its log 2 fold change, p-value and adjusted p-value, and the corresponding base mean.
