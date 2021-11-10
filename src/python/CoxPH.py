import SurvivalAnalysis
import os
import pandas as pd
import numpy as np
from lifelines.fitters.coxph_fitter import CoxPHFitter
import patsy
from statsmodels.stats.multitest import multipletests

SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}

class CoxPH(SurvivalAnalysis.SurvivalAnalysis):
    def __init__(self, covariates, strata, **kws):
        super().__init__(**kws)
        
        # process options
        if covariates is not None: covariates = covariates.split(',')
        if strata is not None: strata = strata.split(',')
        
        # instantiate
        self.covariates = covariates
        self.strata = strata
        
        
    def fit_model(self, x, metadata, sample_col, time_col, event_col, covariates, strata, robust=False):
        """
        Fits our survival model.
        """
        # create input df for fitter
        # combining 'x' from 'data' and 'metadata'
        name = x.name
        
        df = pd.DataFrame({
            name: x,
            time_col: metadata.set_index(sample_col)[time_col],
            event_col: metadata.set_index(sample_col)[event_col]
        })
        
        if strata is not None: 
            df[strata] = metadata[strata]
        if covariates is not None:
            for covariate in covariates:
                df = pd.concat([df, pd.get_dummies(metadata[covariate])], axis=1)
        
        # fit
        model = CoxPHFitter()
        
        try:
            model.fit(df.dropna(), time_col, event_col, strata=strata, robust=robust)
            # prepare output
            # coefficient, std, Zscore, p-value, hazard ratio, concordance
            summary = model.summary.loc[name]
            output = pd.Series({
                'coefficient': summary['coef'],
                'hazard_ratio': summary['exp(coef)'],
                'se': summary['se(coef)'],
                'z_score': summary['z'],
                'pvalue': summary['p'],
                'concordance': model.concordance_index_
            })
        except:
            output = pd.Series({
                'coefficient': np.nan,
                'hazard_ratio': np.nan,
                'se': np.nan,
                'z_score': np.nan,
                'pvalue': np.nan,
                'concordance': np.nan
            })

            
        return output
        
        
    def do_survival_analysis(self):
        """
        Run survival fitter for every row in data.
        """
        result = self.data.apply(
            lambda x: self.fit_model(x,                      
                                     self.metadata, 
                                     self.sample_col,
                                     self.time_col, 
                                     self.event_col, 
                                     self.covariates, 
                                     self.strata), 
            axis=1)
        
        # correct for multiple testing
        result['padj'] = np.nan
        is_missing = result['pvalue'].isnull()
        result.loc[~is_missing,'padj'] = multipletests(result.loc[~is_missing,'pvalue'].values, method=self.padj_method)[1]
        
        # log transform p-values
        result['log10_pvalue'] = - np.log10(result['pvalue'])
        result['log10_padj'] = - np.log10(result['padj'])
        
        # add more info about the model
        result['test_func'] = 'CoxPHFitter'
        result['padj_metod'] = self.padj_method
        
        self.survival_analysis = result
        
        
    def prepare_outputs(self):
        """
        Creates the 'result' instance
        """
        result = self.survival_analysis
        result['thresh_std'] = self.thresh_std
        self.result = result.reset_index()


##### SCRIPT #####
def main():
    args = SurvivalAnalysis.parse_args()
    model = CoxPH(
        data_file = args.data_file,
        metadata_file = args.metadata_file,
        sample_col = args.sample_col,
        time_col = args.time_col,
        event_col = args.event_col,
        output_file = args.output_file,
        subset_col = args.subset_col,
        subset_values = args.subset_values,
        padj_method = args.padj_method,
        thresh_std = args.thresh_std,
        covariates = args.covariates,
        strata = args.strata
    )
    model.run()
    model.write_outputs(**SAVE_PARAMS)


if __name__ == '__main__':
    main()
    print('Done!')
