import pandas as pd
import numpy as np
import os
import argparse
import util
import process_inputs


class ABtest:
    """
    Parent class that processes inputs the same way.
    Inputs are common, but outputs are child-dependent.
    """
    def __init__(self,
                 data_file,
                 metadata_file,
                 sample_col,
                 comparison_col,
                 condition_a,
                 condition_b,
                 output_file=None,
                 output_dir='.',
                 subset_col=None,
                 subset_values=None,
                 padj_method='fdr_bh',
                 thresh_std=0):
        
        self.data_file = data_file
        self.metadata_file = metadata_file
        self.sample_col = sample_col
        self.comparison_col = comparison_col
        self.condition_a = condition_a
        self.condition_b = condition_b
        self.conditions = [self.condition_a, self.condition_b]
        self.output_file = output_file
        self.output_dir = output_dir
        self.subset_col = subset_col
        self.subset_values = subset_values
        self.padj_method = padj_method
        self.thresh_std = thresh_std

        
    def get_samples_per_condition(self):
        """
        """
        samples = []
        for condition in self.conditions:
            idx = self.metadata[self.comparison_col] == condition
            samples.append(list(self.metadata.loc[idx, self.sample_col]))
        
        # add as class variable
        self.samples_a = samples[0]
        self.samples_b = samples[1]
        
        
    def prepare_inputs(self):
        """
        Read files, subset metadata according to columns in data, filter data rows by std.
        """
        self.data, self.metadata = process_inputs.read_files(self.data_file, self.metadata_file)
        self.data, self.thresh_std = process_inputs.filter_data_by_std(self.data, self.thresh_std)
        self.data, self.metadata = process_inputs.prepare_data(self.data, self.metadata, self.sample_col)
        
        if self.subset_col is not None:
            print('Subsetting on',self.subset_col)
            self.data, self.metadata = process_inputs.subset_data(self.data, self.metadata, self.sample_col, self.subset_col, self.subset_values)
            
        self.get_samples_per_condition() # adds samples_a and samples_b
        
        
    def compute_conditions_diff(self):
        """
        compute mean and median difference and log2FC between conditions.
        """
        # mean
        mean_diff = self.data[self.samples_a].mean(1) - self.data[self.samples_b].mean(1)
        mean_log2FC = np.log2(self.data[self.samples_a].mean(1) / self.data[self.samples_b].mean(1))

        # median
        median_diff = self.data[self.samples_a].median(1) - self.data[self.samples_b].median(1)
        median_log2FC = np.log2(self.data[self.samples_a].median(1) / self.data[self.samples_b].median(1))

        # prepare result
        result = pd.DataFrame({
            'mean_diff':mean_diff,
            'mean_log2FC':mean_log2FC,
            'median_diff':median_diff,
            'median_log2FC':median_log2FC
        })
        self.conditions_diff = result
        
        
    def run(self, **kws):
        """
        Run full analysis
        """
        print('Preparing inputs...')
        self.prepare_inputs()        
        print('Comparing conditions %s vs %s' % (self.condition_a, self.condition_b))
        self.compare_conditions(**kws)
        print('Preparing outputs...')
        self.prepare_outputs()
        
        
    def create_required_dirs(self, file_name):
        util.create_required_dirs(file_name)
    
    
    def write_table(self, df, file_name, sep='\t', **kws):
        util.write_table(df, file_name, **kws)


def parse_args():
    """
    For scripting.
    """
    # required
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_file', type=str, required=True)
    parser.add_argument('--metadata_file', type=str, required=True)
    parser.add_argument('--sample_col', type=str, required=True)
    parser.add_argument('--comparison_col', type=str, required=True)
    parser.add_argument('--condition_a', type=str, required=True)
    parser.add_argument('--condition_b', type=str, required=True)
    parser.add_argument('--output_file',type=str, default=None)
    parser.add_argument('--output_dir', type=str)
    # optional
    parser.add_argument('--subset_col', type=str, default=None)
    parser.add_argument('--subset_values', type=str, default=None)
    parser.add_argument('--padj_method', type=str, default='fdr_bh')
    parser.add_argument('--thresh_std', type=int, default=0)
    parser.add_argument('--formula', type=str, default='')
    

    args = parser.parse_args()
    return args
