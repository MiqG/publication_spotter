import pandas as pd
import numpy as np
import os
import argparse
import process_inputs

class SurvivalAnalysis:
    """
    Parent class that reads inputs.
    """
    def __init__(self,
                 data_file,
                 metadata_file,
                 sample_col,
                 time_col,
                 event_col,
                 output_file=None,
                 subset_col=None,
                 subset_values=None,
                 padj_method='fdr_bh',
                 thresh_std=0,
                 features_as_rows=True):
        
        self.data_file = data_file
        self.metadata_file = metadata_file
        self.sample_col = sample_col
        self.time_col = time_col
        self.event_col = event_col
        self.output_file = output_file
        self.subset_col = subset_col
        self.subset_values = subset_values
        self.padj_method = padj_method
        self.thresh_std = thresh_std
        self.features_as_rows = features_as_rows
        
        
    def prepare_inputs(self):
        """
        Read files, subset metadata according to columns in data, filter data rows by std.
        """
        # read data
        self.data, self.metadata = process_inputs.read_files(self.data_file, self.metadata_file)
        # remove rows with std < 'thresh_std'
        self.data, self.thresh_std = process_inputs.filter_data_by_std(self.data, self.thresh_std)
        # subset metadata
        self.data, self.metadata = process_inputs.prepare_data(self.data, self.metadata, self.sample_col)
        
        if self.subset_col is not None:
            print('Subsetting on',self.subset_col)
            self.data, self.metadata = process_inputs.subset_data(self.data, self.metadata, self.sample_col, self.subset_col, self.subset_values)
    
    
    def run(self, **kws):
        """
        Run full analysis
        """
        print('Preparing inputs...')
        self.prepare_inputs()        
        print('Performing survival analysis...')
        self.do_survival_analysis(**kws)
        print('Preparing outputs...')
        self.prepare_outputs()
        
    def write_table(self, df, file_name, **kws):
        df.to_csv(file_name, **kws)
    
    def write_outputs(self, **kws):
        self.write_table(self.result, self.output_file, **kws)


def parse_args():
    """
    For scripting.
    """
    # required
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_file', type=str, required=True)
    parser.add_argument('--metadata_file', type=str, required=True)
    parser.add_argument('--sample_col', type=str, required=True)
    parser.add_argument('--time_col', type=str, required=True)
    parser.add_argument('--event_col', type=str, required=True)
    parser.add_argument('--output_file',type=str, default=None)
    # optional
    parser.add_argument('--features_as_rows',type=bool,default=True)
    parser.add_argument('--subset_col', type=str, default=None)
    parser.add_argument('--subset_values', type=str, default=None)
    parser.add_argument('--padj_method', type=str, default='fdr_bh')
    parser.add_argument('--thresh_std', type=int, default=0)
    parser.add_argument('--covariates', type=str, default=None)
    parser.add_argument('--strata', type=str, default=None)
    # PCARSF
    parser.add_argument('--pca_kws', type=str, default="{}")
    parser.add_argument('--max_explained_variance', type=float, default=0.9)
    parser.add_argument('--rsf_kws', type=str, default="{}")
    parser.add_argument('--varimp_kws', type=str, default="{}")
    parser.add_argument('--pca_loadings_file',type=str)
    parser.add_argument('--variable_importances_file',type=str)

    args = parser.parse_args()
    return args
