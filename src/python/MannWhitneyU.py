import PairwiseTest
import scipy
import pandas as pd
import numpy as np
import os
from statsmodels.stats.multitest import multipletests
from scipy.stats import kurtosis, skew, median_abs_deviation, iqr

class MannWhitneyU(PairwiseTest.PairwiseTest):
    def __init__(self, **kws):
        super().__init__(**kws)
        self.test_func = scipy.stats.mannwhitneyu
    
    
    def ab_test(self, a, b, **kws_test):
        # do not compute if less than 15 missing
        a = a[~np.isnan(a)]
        b = b[~np.isnan(b)]
        if ((len(a) > 0) & (len(b) > 0)):
            try:
                statistic, pvalue = self.test_func(a, b, **kws_test)
            except:
                statistic, pvalue = np.nan, np.nan
        else:
            statistic, pvalue = np.nan, np.nan
        return pd.Series({'statistic':statistic,'pvalue':pvalue})


    def compute_ab_test(self, **kws_test):
        """
        A/B test across pd.DataFrame.
        padj_method: same as in `statsmodels.stats.multitest.multipletests`
        """
        result = self.data.apply(lambda x: self.ab_test(x[self.samples_a],x[self.samples_b], **kws_test), axis=1)
        result['log10_pvalue'] = - np.log10(result['pvalue'])
        
        # correct for multiple testing
        result['padj'] = np.nan
        is_missing = result['pvalue'].isnull()
        
        if sum(is_missing)!=len(result):
            result.loc[~is_missing,'padj'] = multipletests(result.loc[~is_missing,'pvalue'].values, method=self.padj_method)[1]
            result['log10_padj'] = - np.log10(result['padj'])
    
        # add more info about the test
        result['test_func'] = self.test_func.__name__
        result['padj_metod'] = self.padj_method
        
        self.conditions_test = result
    
    
    def compare_conditions(self, **kws_test):
        self.compute_conditions_diff()
        self.compute_ab_test(**kws_test)
        
    
    def prepare_outputs(self):
        """
        Creates the 'result' instance.
        """
        result = self.conditions_diff.join(self.conditions_test)
        result['comparison_col'] = self.comparison_col
        # summary stats
        result['condition_a'] = self.condition_a
        result['condition_a-median'] = self.data[self.samples_a].median(1)
        result['condition_a-std'] = self.data[self.samples_a].std(1)
        result['condition_a-mad'] = median_abs_deviation(self.data[self.samples_a], axis=1, nan_policy='omit')
        result['condition_a-iqr'] = iqr(self.data[self.samples_a], axis=1, nan_policy='omit')
        result['condition_a-skewness'] = skew(self.data[self.samples_a],axis=1, nan_policy='omit')
        result['condition_a-kurtosis'] = kurtosis(self.data[self.samples_a], axis=1, nan_policy='omit')
        result['condition_a-n_present'] = (~self.data[self.samples_a].isna()).sum(1)
        result['condition_a-n_nan'] = self.data[self.samples_a].isna().sum(1)
        result['condition_a-perc_nan'] = self.data[self.samples_a].isna().sum(1) / len(self.samples_a)
        
        result['condition_b'] = self.condition_b
        result['condition_b-median'] = self.data[self.samples_b].median(1)
        result['condition_b-std'] = self.data[self.samples_b].std(1)
        result['condition_b-mad'] = median_abs_deviation(self.data[self.samples_b], axis=1, nan_policy='omit')
        result['condition_b-iqr'] = iqr(self.data[self.samples_b], axis=1, nan_policy='omit')
        result['condition_b-skewness'] = skew(self.data[self.samples_b],axis=1, nan_policy='omit')
        result['condition_b-kurtosis'] = kurtosis(self.data[self.samples_b], axis=1, nan_policy='omit')
        result['condition_b-n_present'] = (~self.data[self.samples_b].isna()).sum(1)
        result['condition_b-n_nan'] = self.data[self.samples_b].isna().sum(1)
        result['condition_b-perc_nan'] = self.data[self.samples_b].isna().sum(1) / len(self.samples_b)
        
        result['thresh_std'] = self.thresh_std
        
        self.result = result.reset_index()

        
    def write_outputs(self):
        """
        Write outputs as files.
        """
        if self.output_file is None:
            output_filename = ('-'.join(['ab_test',
                                    self.test_func.__name__,
                                    self.condition_a+'_vs_'+self.condition_b
                                   ])+'.tsv').replace(' ','')
            self.output_file = os.path.join(self.output_dir, output_filename)

        self.create_required_dirs(self.output_file)
        self.write_table(self.result, self.output_file, index=False)

        
##### SCRIPT #####
def main():
    args = PairwiseTest.parse_args()
    ab_test = MannWhitneyU(
        data_file = args.data_file,
        metadata_file = args.metadata_file,
        sample_col = args.sample_col,
        comparison_col = args.comparison_col,
        condition_a = args.condition_a,
        condition_b = args.condition_b,
        output_file = args.output_file,
        output_dir = args.output_dir,
        subset_col = args.subset_col,
        subset_values = args.subset_values,
        padj_method = args.padj_method,
        thresh_std = args.thresh_std    
    )
    ab_test.run()
    ab_test.write_outputs()
    
    
if __name__ == '__main__':
    main()
    print('Done!')