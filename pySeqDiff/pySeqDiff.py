import pandas as pd
import statsmodels.stats.multitest as smm
from sklearn.preprocessing import scale
from sklearn.decomposition import PCA
import argparse
import numpy as np
from scipy.stats import gmean, nbinom, poisson
from scipy.stats import ttest_ind
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

def normalize_counts(df, pseudocount=0.01):
    df = df.replace(0, pseudocount)
    # Find the geometric mean of samples for each gene
    geomeans = gmean(df, axis=1)
    
    df_ratios = df.divide(geomeans, axis = 0)
    sfs = df_ratios.median(axis = 0)
    df = df.divide(sfs, axis = 1)

    return df

def base_means(df):
    return df.mean(axis=1).tolist()

def batch_correction(df, n_components=5, pseudocount=0.01):
    # Identify constant rows
    constant_rows = df.index[df.nunique(axis=1) <= 1]
    
    # Separate constant and variable rows
    df_variable = df.drop(constant_rows)
    df_constant = df.loc[constant_rows]
    
    # Perform PCA and batch correction on variable rows
    df_T = df_variable.T
    pca = PCA(n_components=n_components)
    pca_result = pca.fit_transform(scale(df_T))
    corrected_T = df_T - pca_result.dot(pca.components_)
    corrected_variable = corrected_T.T
    
    corrected_variable = corrected_variable.clip(lower=pseudocount)

    # Concatenate variable and constant rows
    corrected = pd.concat([corrected_variable, df_constant])
    
    return corrected

def calculate_fold_changes(control_data, treatment_data, pseudocount=0.01):
    fc = (treatment_data.mean(axis=1) + pseudocount) / (control_data.mean(axis=1) + pseudocount)
    if fc.min() <= 0:
        print('Invalid values in fold changes:', fc[fc <= 0])
    log2_fc = np.log2(fc)
    return log2_fc

def calculate_p_values(control_data, treatment_data):
    p_values = []
    for gene_id in control_data.index:
        control_mean = control_data.loc[gene_id].mean()
        control_var = control_data.loc[gene_id].var()
        treatment_mean = treatment_data.loc[gene_id].mean()

        control_var_shrunk = 0.25 * control_var + 0.25 * control_data.var().mean()

        # Estimating the parameters of the negative binomial distribution
        if control_var_shrunk > control_mean:
            n = control_mean**2 / (control_var_shrunk - control_mean)
            p = n / (n + control_mean)
            
            # Calculate the two-sided p-value as this study cares for both up and down regulated genes
            if treatment_mean > control_mean:
                p_value = 2 * nbinom.sf(treatment_mean - 1, n, p)
            else:
                p_value = 2 * nbinom.cdf(treatment_mean, n, p)

        else:
            # Calculate the two-sided p-value for the Poisson distribution
            if treatment_mean > control_mean:
                p_value = 2 * poisson.sf(treatment_mean - 1, control_mean)
            else:
                p_value = 2 * poisson.cdf(treatment_mean, control_mean)

        p_values.append(p_value)

    return pd.Series(p_values, index=control_data.index)

def correct_p_values(p_values, padj_type):
    if (padj_type == 'bonferroni'):
        _, p_values_corrected, _, _ = smm.multipletests(p_values, method='bonferroni')
        return p_values_corrected
    else:
        _, p_values_corrected, _, _ = smm.multipletests(p_values, method='fdr_bh')
        return p_values_corrected

def differential_expression_analysis(countdata, metadata, base_means, padj_type, sorted):
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
            'log2fold_change': fold_changes,
            'p_value': p_values,
            'p_value_corrected': p_values_corrected,
            'base_mean': base_means
        })
        if sorted:
            # sort the results by increasing 'p_value_corrected'
            result.sort_values(by='p_value_corrected', inplace=True)
        results.append(result)
    return results

def save_results(results, output_file):
    result_df = pd.concat(results)
    result_df.to_csv(output_file, sep='\t', index=False)

def volcano_plot(results, pval_thresh):
    result_df = pd.concat(results)

    #plt.figure(figsize=(10,10))

    fc_threshold = 0
    plt.ylim(0, 50)
    plt.xlabel('log2FoldChange') #x label
    plt.ylabel('-log10pvalue') #y label
    plt.title('Differential Expression')

    #plot scatter plot while labeling significant expression with red based on padj<0.05 AND log2fold > 0 or log2fold < 0
    plt.scatter(result_df["log2fold_change"], -np.log10(result_df["p_value_corrected"]), c = np.where((result_df["p_value_corrected"] < pval_thresh)  & ((result_df["log2fold_change"] < 0)| (result_df["log2fold_change"] > 0)), "red", "black"))
    #plt.scatter(result_df["log2fold_change"], -np.log10(result_df["p_value_corrected"]))
    #using threshold of p value 5% and taking top 10
    result_df = result_df[result_df['p_value_corrected'] < pval_thresh]
    pSort =  result_df.sort_values('p_value_corrected')
    top_ten = pSort.iloc[:10]

    for index, row in top_ten.iterrows():
        plt.annotate(row["gene_id"], (row["log2fold_change"], -np.log10(row["p_value"])))
        
    plt.savefig("VolcanoPlot.png")
    plt.show()
    
def ma_plot(results, pval_thresh):
    result_df = pd.concat(results)

    plt.figure(figsize=(10,10))

    #color points red if p value is less than 0.1 or choose other threshold value
    plt.ylim(-2, 2)
    plt.xlabel('Mean of Normalized Counts') #x label
    plt.ylabel('log2FoldChange') #y label
    plt.title('MA-Plot')
    plt.savefig("MAPlot.png")
    #plt.scatter(result_df["log2fold_change"], result_df["base_mean"], c = np.where((result_df["p_value_corrected"] < pval_thresh), "red", "black"))
    plt.scatter(result_df["log2fold_change"], result_df["base_mean"])
    plt.savefig("MAPlot.png")
    plt.show()

def main():
    myParser = argparse.ArgumentParser(description='Local alignment program')
    myParser.add_argument('-c', '--countdata', help="Input counts txt file", type=str)
    myParser.add_argument('-o', '--output_file', help="Write output to file", type=str)
    myParser.add_argument('-padj', '--pvalue_adjusted', help="P value adjustment type. " "Default: fdr_bh",type=str)
    myParser.add_argument('-pval_thresh', '--pvalue_threshold', help="P value threshold for MA-Plot" ,type=float)
    myParser.add_argument('-s', '--sorted_output', action='store_true')
    inputArgs = myParser.parse_args()
    
    #countdata = load_data("~/GSE221626_counts.txt")
    countdata = load_data(inputArgs.countdata)
    outfile = inputArgs.output_file
    pval_adj = inputArgs.pvalue_adjusted
    pval_thresh = inputArgs.pvalue_threshold
    metadata = define_metadata(countdata)
    countdata = normalize_counts(countdata)
    countdata = batch_correction(countdata)
    base_means_values = base_means(countdata)
    results = differential_expression_analysis(countdata, metadata, base_means_values, pval_adj, sorted)
    save_results(results, outfile)
    volcano_plot(results, pval_thresh)
    ma_plot(results, pval_thresh)

if __name__ == "__main__":
    main()