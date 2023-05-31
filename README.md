# Basic Usage Instructions:

# Setup:
``` git clone https://github.com/andrehchu/RNASeq ```

``` cd RNASeq ```

# Usage Instructions:
Make sure the necessary packages are downloaded before running the program!
If necessary: use the following commands:
    1. pandas: ```pip install pandas```
    2. statsmodels: ```pip install statsmodels```
    3. scikit-learn: ```pip install scikit-learn```
    4. argparse: ```pip install argparse```
    5. numpy: ```pip install numpy```
    6. scipy: ```pip install scipy```
``` python3 207analysis.py -c [RNA-seq count data file] -o [output_file] -s [sorted output if present]  ```
# Use the count data file I have provided and run the command:
``` python3 207analysis.py -c GSE207721_raw_counts_GRCh38.p13_NCBI.tsv -o [output file name of your choice] -s ```

# How to interpret the output:
 The count data file used for testing is obtained from this study 'https://pubmed.ncbi.nlm.nih.gov/35932580/'. The output shows the gene ID, if it is experimental or wildtype, its log 2 fold change, p-value and adjusted p-value, and the corresponding base mean.