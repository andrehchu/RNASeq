#### THE CODE BELOW WAS RAN ON UCSD DATAHUB!

%%R

# Load the libraries we need
library("DESeq2")

# Import count data
countdata <- read.table("~/GSE207721_raw_counts_GRCh38.p13_NCBI.tsv", header=TRUE, row.names=1, sep = '\t')

# Define metadata with gene identifiers as row names
metadata <- data.frame(
  condition = factor(c(rep("WT", 2), rep("Experimental", 2), rep("WT", 3), rep("Experimental", 2)))
)
rownames(metadata) <- colnames(countdata)

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = metadata, design = ~ condition)

# Run DESeq
dds <- DESeq(dds)

# Obtain and save results to a text file
res <- results(dds)
res_sorted <- res[order(res$pvalue), ]

# Save sorted results to a text file
write.table(res_sorted, file = "~/benchmark207_sorted.txt", sep = "\t", quote = FALSE, row.names = TRUE)