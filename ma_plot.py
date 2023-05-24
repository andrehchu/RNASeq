import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

result_df = pd.concat(results)

plt.figure(figsize=(10,10))

pvalue_threshold = 0.05


#color points red if p value is less than 0.1 or choose other threshold value

plt.xlabel('Mean of Normalized Counts') #x label
plt.ylabel('log2FoldChange') #y label
plt.title('TBD')

plt.scatter(result_df["log2fold_change"], result_df["mean_normalizedcounts"], c = np.where((result_df["p_value_corrected"] < pvalue_threshold), "red", "black"))
