import pandas as pd
import numpy as np
from scipy import stats
import statsmodels.stats.multitest as smm
from sklearn.preprocessing import scale
from sklearn.decomposition import PCA
import argparse

def load_data(file_path):
    df = pd.read_csv(file_path, sep='\t', index_col=0)
    if df.isnull().sum().sum() > 0:
        raise ValueError("Input data contains missing values.")
    if df.nunique(axis=1).min() == 1:
        raise ValueError("Input data contains rows with constant values.")
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
    pca_result = pca.fit_transform(scale(df.T))
    corrected = df - pca_result.T.dot(pca.components_)
    return corrected

def calculate_fold_changes(control_data, treatment_data):
    return treatment_data.mean(axis=1) / control_data.mean(axis=1)

def calculate_p_values(control_data, treatment_data):
    _, p_values = stats.ttest_ind(control_data.T, treatment_data.T, equal_var=False)
    return p_values

def correct_p_values(p_values):
    _, p_values_corrected, _, _ = smm.multipletests(p_values, method='fdr_bh')
    return p_values_corrected

def differential_expression_analysis(countdata, metadata):
    results = []
    for condition in metadata['condition'].unique():
        if condition == "WT":
            continue
        control_data = countdata.loc[:, metadata['condition'] == 'WT']
        treatment_data = countdata.loc[:, metadata['condition'] == condition]
        fold_changes = calculate_fold_changes(control_data, treatment_data)
        p_values = calculate_p_values(control_data, treatment_data)
        p_values_corrected = correct_p_values(p_values)
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

def main():
    myParser = argparse.ArgumentParser(description='Local alignment program')
    myParser.add_argument('count_data', '--cd', type=str)
    myParser.add_argument('output_file', '--out', type=int)
    
    countdata = load_data("~/GSE221626_counts.txt")
    metadata = define_metadata(countdata)
    countdata = normalize_counts(countdata)
    countdata = batch_correction(countdata)
    results = differential_expression_analysis(countdata, metadata)
    save_results(results, "~/benchmark.txt")

if __name__ == "__main__":
    main()