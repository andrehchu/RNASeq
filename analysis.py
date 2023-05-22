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
    return df

def define_metadata(df):
    metadata = pd.DataFrame({
        'condition': ["WT"] * 3 + ["Adar1KI"] * 3 + ["Adar1_p110KI"] * 3 + ["Adar1_p150KI"] * 3 + ["E861AKI"] * 3
    }, index=df.columns)
    return metadata

def normalize_counts(df, pseudocount=1):
    df = df.replace(0, pseudocount)
    gm = df.prod(axis=1) ** (1.0 / len(df.columns))
    df = df.div(gm, axis=0)
    size_factors = df.median(axis=0)
    df = df.div(size_factors, axis=1)
    return df

def batch_correction(df, n_components=2):
    pca = PCA(n_components=n_components)
    df_T = df.T
    pca_result = pca.fit_transform(scale(df_T))
    corrected_T = df_T - pca_result.dot(pca.components_)
    corrected = corrected_T.T
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
        result = pd.DataFrame({
            'gene_id': countdata.index,
            'condition': condition,
            'fold_change': fold_changes,
            'p_value': p_values,
            'p_value_corrected': p_values_corrected,
        })
        results.append(result)
    return results

def save_results(results, output_file):
    result_df = pd.concat(results)
    result_df.to_csv(output_file, sep='\t', index=False)

def volcano_plot(results):
    result_df = pd.concat(results)

    plt.figure(figsize=(10,10))

    pvalue_threshold = 0.05
    fc_threshold = 0

    result_df = result_df.rename({"Unnamed: 0" : "gene_id"}, axis='columns')


    plt.xlabel('log2FoldChange') #x label
    plt.ylabel('-log10pvalue') #y label
    plt.title('Differential Expression')

    #plot scatter plot while labeling significant expression with red based on padj<0.05 AND log2fold > 0 or log2fold < 0
    plt.scatter(deseqresults['log2FoldChange'], -np.log10(deseqresults['padj']), c = np.where((deseqresults["padj"] < 0.05)  & ((deseqresults["log2FoldChange"] < 0) | (deseqresults["log2FoldChange"] > 0)), "red", "black"))

    #using threshold of p value 5% and taking top 10
    deseqresults = deseqresults[deseqresults['p_value_corrected'] < 0.05]
    pSort =  deseqresults.sort_values('p_value_corrected')
    top_ten = pSort.iloc[:10]

    for index, row in top_ten.iterrows():
        plt.annotate(row["gene_id"], (row["log2FoldChange"], -np.log10(row["pvalue"])))
        
    plt.savefig("VolcanoPlotDeseq2.png")
    



def main():
    myParser = argparse.ArgumentParser(description='Local alignment program')
    myParser.add_argument('-c', '--countdata', type=str)
    myParser.add_argument('-o', '--output_file', type=str)
    myParser.add_argument('-padj', '--pvalue_adjusted', type=str)
    inputArgs = myParser.parse_args()
    
    #countdata = load_data("~/GSE221626_counts.txt")
    countdata = load_data(inputArgs.countdata)
    outfile = inputArgs.output_file
    pval_adj = inputArgs.pvalue_adjusted
    metadata = define_metadata(countdata)
    countdata = normalize_counts(countdata)
    countdata = batch_correction(countdata)
    results = differential_expression_analysis(countdata, metadata, pval_adj)
    save_results(results, outfile)

if __name__ == "__main__":
    main()