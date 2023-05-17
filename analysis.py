import pandas as pd
import numpy as np
from scipy import stats
import statsmodels.stats.multitest as smm
import argparse
from sklearn.preprocessing import scale
from sklearn.decomposition import PCA

class DifferentialExpression:
    def __init__(self, data_file, metadata_file):
        self.df = self.load_data(data_file)
        self.metadata = self.load_metadata(metadata_file)
        self.check_data_quality()

    def load_data(self, data_file):
        df = pd.read_csv(data_file, sep='\t', index_col=0)
        return df

    def load_metadata(self, metadata_file):
        metadata = pd.read_csv(metadata_file, sep='\t', index_col=0)
        return metadata

    def check_data_quality(self):
        if self.df.isnull().sum().sum() > 0:
            raise ValueError("Input data contains missing values.")
        if self.df.nunique(axis=1).min() == 1:
            raise ValueError("Input data contains rows with constant values.")

    def normalize_counts(self, df, pseudocount=1):
        df = df.replace(0, pseudocount)
        gm = df.prod(axis=1) ** (1.0 / len(df.columns))
        df = df.div(gm, axis=0)
        size_factors = df.median(axis=0)
        df = df.div(size_factors, axis=1)
        return df

    def batch_correction(self, df, n_components=2):
        pca = PCA(n_components=n_components)
        pca_result = pca.fit_transform(scale(df.T))
        corrected = df - pca_result.T.dot(pca.components_)
        return corrected

    def calculate(self):
        df = self.normalize_counts(self.df)
        df = self.batch_correction(df)
        self.results = []
        for condition in self.metadata['condition'].unique():
            control_df = df.loc[:, self.metadata['condition'] == 'control']
            treatment_df = df.loc[:, self.metadata['condition'] == condition]
            results = self.differential_expression(control_df, treatment_df, condition)
            self.results.append(results)

    def differential_expression(self, control_df, treatment_df, condition):
        fold_changes = treatment_df.mean(axis=1) / control_df.mean(axis=1)
        control_var = control_df.var(axis=1) + 1e-8
        treatment_var = treatment_df.var(axis=1) + 1e-8
        t, p_values = stats.ttest_ind_from_stats(
            control_df.mean(axis=1), np.sqrt(control_var), len(control_df.columns),
            treatment_df.mean(axis=1), np.sqrt(treatment_var), len(treatment_df.columns),
            equal_var=False
        )
        p_values_corrected = self.adjust_p_values(p_values)
        results = pd.DataFrame({
            'gene_id': control_df.index,
            'condition': condition,
            'fold_change': fold_changes,
            'p_value': p_values,
            'p_value_corrected': p_values_corrected,
        })
        return results

    def adjust_p_values(self, p_values):
        reject, pvals_corrected, _, _ = smm.multipletests(p_values, method='fdr_bh')
        return pvals_corrected

    def save_results(self, filename):
        result_df = pd.concat(self.results)
        result_df.to_csv(filename, index=False, sep='\t')

def main():
    parser = argparse.ArgumentParser(description='Calculate differential gene expression.')
    parser.add_argument('data_file', type=str, help='Path to the input data file.')
    parser.add_argument('metadata_file', type=str, help='Path to the input metadata file.')
    parser.add_argument('output_file', type=str, help='Path to the output file.')
    args = parser.parse_args()

    de = DifferentialExpression(args.data_file, args.metadata_file)
    de.calculate()
    de.save_results(args.output_file)

if __name__ == "__main__":
    main()
