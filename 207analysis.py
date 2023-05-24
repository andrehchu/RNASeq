import pandas as pd
from scipy import stats
import statsmodels.stats.multitest as smm
from sklearn.preprocessing import scale
from sklearn.decomposition import PCA
import argparse
import numpy as np

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

def normalize_counts(df, pseudocount=1):
    df = df.replace(0, pseudocount)
    gm = df.prod(axis=1) ** (1.0 / len(df.columns))
    df = df.div(gm, axis=0)
    size_factors = df.median(axis=0)
    df = df.div(size_factors, axis=1)
    if (df < 0).any().any():
        print("Negative values after normalization: ", df[df < 0])
    return df

def batch_correction(df, n_components=2, pseudocount=0.1):
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

def calculate_fold_changes(control_data, treatment_data, pseudocount=0.1):
    fc = (treatment_data.mean(axis=1) + pseudocount) / (control_data.mean(axis=1) + pseudocount)
    if fc.min() <= 0:
        print('Invalid values in fold changes:', fc[fc <= 0])
    log2_fc = np.log2(fc)
    return log2_fc

def calculate_p_values(control_data, treatment_data):
    p_values = []
    for gene_id in control_data.index:
        _, p_value = stats.mannwhitneyu(control_data.loc[gene_id], treatment_data.loc[gene_id], alternative='two-sided')
        p_values.append(p_value)
    return pd.Series(p_values, index=control_data.index)

def correct_p_values(p_values):
    _, p_values_corrected, _, _ = smm.multipletests(p_values, method='fdr_bh')
    return p_values_corrected

def differential_expression_analysis(countdata, metadata, sorted):
    results = []
    for condition in metadata['condition'].unique():
        if condition == "WT":
            continue
        control_data = countdata.loc[:, metadata['condition'] == 'WT']
        treatment_data = countdata.loc[:, metadata['condition'] == condition]
        fold_changes = calculate_fold_changes(control_data, treatment_data)
        p_values = calculate_p_values(control_data, treatment_data)
        p_values_corrected = correct_p_values(p_values.values)
        result = pd.DataFrame({
            'gene_id': countdata.index,
            'condition': condition,
            'log2_fold_change': fold_changes,
            'p_value': p_values,
            'p_value_corrected': p_values_corrected,
        })
        if sorted:
            # sort the results by increasing 'p_value_corrected'
            result.sort_values(by='p_value_corrected', inplace=True)
        results.append(result)
    return results

def save_results(results, output_file):
    result_df = pd.concat(results)
    result_df.to_csv(output_file, sep='\t', index=False)

def main():
    myParser = argparse.ArgumentParser(description='Local alignment program')
    myParser.add_argument('-c', '--countdata', type=str)
    myParser.add_argument('-o', '--output_file', type=str)
    myParser.add_argument('-s', '--sorted_output', type=bool)
    inputArgs = myParser.parse_args()
    sorted = inputArgs.sorted

    countdata = load_data(inputArgs.countdata)
    outfile = inputArgs.output_file
    metadata = define_metadata(countdata)
    countdata = normalize_counts(countdata)
    countdata = batch_correction(countdata)
    results = differential_expression_analysis(countdata, metadata, sorted)
    save_results(results, outfile)

if __name__ == "__main__":
    main()