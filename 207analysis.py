import pandas as pd
import numpy as np
from scipy import stats
import statsmodels.stats.multitest as smm
from sklearn.preprocessing import scale
from sklearn.decomposition import PCA
import argparse
import matplotlib.pyplot as plt

def load_data(file_path):
    # read in count data
    df = pd.read_csv(file_path, sep='\t', index_col=0)
    # for our dataset, we drop the final col because it has only one variant
    df = df.iloc[:, :-1] 
    # add pseudocount for values of zero
    df.replace(0, 0.1, inplace=True) 
    if df.isnull().sum().sum() > 0:
        raise ValueError("Input data contains missing values.")
    #if df.nunique(axis=1).min() == 1:
      #  raise ValueError("Input data contains rows with constant values.")
    if (df < 0).any().any():
        print("Negative values in raw data: ", df[df < 0])
    return df

def define_metadata(df):
    metadata = pd.DataFrame({
        'condition': ['WT' if col in ['GSM6311037', 'GSM6311038', 'GSM6311041', 'GSM6311042', 'GSM6311043'] else 'Experimental' for col in df.columns]
    }, index=df.columns)
    return metadata

def normalize_counts(df, psedudocount=1):
    df = df.replace(0, psedudocount)
    # Find the geometric mean of samples for each gene
    geomeans = []
    for index, row in df.iterrows():
        geomean = 1
        for sample in df.columns[1:].values.tolist():
            geomean *= row[sample]
        geomean = np.sqrt(geomean)
        geomeans.append(geomean)

    # Find the ratios (sample/ref) for each sample
    df_ratios = df.copy()
    for index, row in df_ratios.iterrows():
        for sample in df_ratios.columns[1:].values.tolist():
            df_ratios.at[index,sample] = (df_ratios.at[index,sample]) / geomeans[index]

    # Find size factors for each sample/col
    sfs = []
    for column in df_ratios.columns[1:]:
        sfs.append(np.median(df_ratios[column].tolist()))

    # Divide each sample value by the size factor to normalize
    i = 0
    for column in df.columns[1:].values.tolist():
        df[column] = df[column].div(sfs[i])
        i = i + 1

    return df

def base_means(df):
    baseMeans = []
    for index,row in df.iterrows:
        sum = 0
        numSamples = 0
        for column in df.columns:
            sum += row[column]
            numSamples += 1
        baseMeans.append(sum/numSamples)
    return baseMeans

def batch_correction(df, n_components=2, pseudocount=0.1):
    #Identify constant rows
    constant_rows = df.index[df.nununiqe(axis=1) <= 1]

    #Separate constant and variable rows
    df_variable = df.drop(constant_rows)
    df_constant = df.loc[constant_rows]

    #Perform PCA and batch correction on variable rows
    df_T = df_variable.T
    pca = PCA(n_components=n_components)
    pca_result = pca.fit_transform(scale(df_T))
    corrected_T = df_T - pca_result.dot(pca.components_)
    corrected_variable = corrected_T.T

    corrected_variable = corrected_variable.clip(lower=pseudocount)

    #concatenate variable and constant rows
    corrected = pd.concat([corrected_variable, df_constant])

    return corrected

def calculate_fold_changes(control_data, treatment_data):
    return np.log2(treatment_data.mean(axis=1) / control_data.mean(axis=1))

def calculate_p_values(control_data, treatment_data):
    p_values = []
    for gene_id in control_data.index:
        _, p_value = stats.mannwhitneyu(control_data.loc[gene_id], treatment_data.loc[gene_id], alternative='two-sided')
        p_values.append(p_value)
    return pd.Series(p_values, index=control_data.index)

def correct_p_values(p_values, padj_type):
    if (padj_type == 'bonferroni'):
        _, p_values_corrected, _, _ = smm.multipletests(p_values, method='bonferroni')
        return p_values_corrected
    else:
        _, p_values_corrected, _, _ = smm.multipletests(p_values, method='fdr_bh')
        return p_values_corrected

def differential_expression_analysis(countdata, metadata, padj_type):
    results = []
    for condition in metadata['condition'].unique():
        if condition == "WT":
            continue
        control_data = countdata.loc[:, metadata['condition'] == 'WT']
        treatment_data = countdata.loc[:, metadata['condition'] == condition]
        fold_changes = calculate_fold_changes(control_data, treatment_data)
        p_values = calculate_p_values(control_data, treatment_data)
        p_values_corrected = correct_p_values(p_values.values, padj_type)
        base_means = base_means(countdata)
        result = pd.DataFrame({
            'gene_id': countdata.index,
            'condition': condition,
            'log2fold_change': fold_changes,
            'p_value': p_values,
            'p_value_corrected': p_values_corrected,
            'base_means': base_means,
        })
        results.append(result)
    return results

def save_results(results, output_file):
    result_df = pd.concat(results)
    result_df.to_csv(output_file, sep='\t', index=False)

def volcano_plot(results, pval_thresh):
    result_df = pd.concat(results)

    plt.figure(figsize=(10,10))

    fc_threshold = 0

    plt.xlabel('log2FoldChange') #x label
    plt.ylabel('-log10pvalue') #y label
    plt.title('Differential Expression')

    #plot scatter plot while labeling significant expression with red based on padj<0.05 AND log2fold > 0 or log2fold < 0
    plt.scatter(result_df["log2fold_change"], -np.log10(result_df["p_value_corrected"]), c = np.where((result_df["p_value_corrected"] < pval_thresh)  & ((result_df["log2fold_change"] < 0) | (result_df["log2fold_change"] > 0)), "red", "black"))

    #using threshold of p value 5% and taking top 10
    result_df = result_df[result_df['p_value_corrected'] < 0.05]
    pSort =  result_df.sort_values('p_value_corrected')
    top_ten = pSort.iloc[:10]

    for index, row in top_ten.iterrows():
        plt.annotate(row["gene_id"], (row["log2fold_change"], -np.log10(row["pvalue"])))
        
    plt.savefig("VolcanoPlot.png")
    plt.show()
    
def ma_plot(results, pval_thresh, lower_ybound, upper_ybound):
    result_df = pd.concat(results)

    plt.figure(figsize=(10,10))

    #color points red if p value is less than 0.1 or choose other threshold value
    plt.ylim(lower_ybound, upper_ybound)
    plt.xlabel('Mean of Normalized Counts') #x label
    plt.ylabel('log2FoldChange') #y label
    plt.title('MA-Plot')
    plt.savefig("MAPlot.png")
    plt.scatter(result_df["log2fold_change"], result_df["mean_normalizedcounts"], c = np.where((result_df["p_value_corrected"] < pval_thresh), "red", "black"))


def main():
    myParser = argparse.ArgumentParser(description='Local alignment program')
    myParser.add_argument('-c', '--countdata', help="Input counts txt file", type=str)
    myParser.add_argument('-o', '--output_file', help="Write output to file", type=str)
    myParser.add_argument('-padj', '--pvalue_adjusted', help="P value adjustment type. " "Default: fdr_bh",type=str)
    myParser.add_argument('-pval_thresh', '--pvalue_threshold', help="P value threshold for MA-Plot" ,type=str)
    inputArgs = myParser.parse_args()
    
    #countdata = load_data("~/GSE221626_counts.txt")
    countdata = load_data(inputArgs.countdata)
    outfile = inputArgs.output_file
    pval_adj = inputArgs.pvalue_adjusted
    pval_thresh = inputArgs.pvalue_threshold
    metadata = define_metadata(countdata)
    countdata = normalize_counts(countdata)
    countdata = batch_correction(countdata)
    results = differential_expression_analysis(countdata, metadata, pval_adj)
    save_results(results, outfile)
    volcano_plot(results, pval_thresh)
    ma_plot(results, pval_thresh)

if __name__ == "__main__":
    main()