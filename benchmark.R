library("DESeq2")

# import count data
countdata <- read.table("~/GSE221626_counts.txt", header=TRUE, row.names=1, sep = '\t')

# remove 'WT_R26-Adar1p150KI' column as it has only one sample
countdata <- countdata[, -ncol(countdata)]
countdata[] <- lapply(countdata, round)

# define metadata
metadata <- data.frame(
  row.names = colnames(countdata),
  condition = factor(c(rep("WT", 3), rep("Adar1KI", 3), rep("Adar1_p110KI", 3), rep("Adar1_p150KI", 3), rep("E861AKI", 3)))
)

# create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = metadata, design = ~ condition)

# run DESeq
dds <- DESeq(dds)

# obtain and save results to text file
res <- results(dds)
write.table(res, file = "~/benchmark.txt", sep = "\t", quote = FALSE, row.names = TRUE)