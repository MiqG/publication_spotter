# Splicing dependency

![](https://fontmeme.com/permalink/210921/2be314f2a86b20b2bc77404529de028a.png)

## Outline
1. [X] Overview of data: 
    - CCLE gene dependencies
    - mRNA levels
    - splicing (maybe split splicing types?)
    - stemness index


2. [X] Fit linear models to predict gene dependency using both splicing and gene expression information.
    - OLS (or MLE with limix??)
        $$ Gene Dependency = \beta_0 + \beta_1*PSI + \beta_2*TPM + \beta_3*PSI*TPM $$
    - use a log-likelihood ratio test to test whether the splicing terms contribute significantly to the prediction.
    
    $H_0: Gene Dependency = \beta_0 + \beta_2*TPM$
    $H_1: Gene Dependency = \beta_0 + \beta_1*PSI + \beta_2*TPM + \beta_3*PSI*TPM$
    
    - summarize them into a single model to get the empirical coefficients
    - sample 500 times the coefficients to be able to make predictions with confidence
    
    - [X] fit with limix using the covariance matrix of dependencies as random effects
        - 100 times on random training/test sets
        - [ ] check that the posterior distributions are normal in general
        - keep the mean, median, quantiles and std of betas
        - keep the mean, median, quantiles and std of pearson correlations
        - which likelihood ratio test counts? --> make likelihood ratio test using robust betas?


3. [X] Model selection and exploration
    - **Model selection**: use different thresholds on LR test's p-value to select the ensemble of linear models that best ranks dependencies and has a positive pearon correlation coefficient. If a more than one model are selected for a gene, pick the one with the worse pearson correlation on the test set.
    - Explore the models:
        - p-value distribution of LR tests w.r.t. threshold selected
        - event std of selected event models?
        - mutation rates of genes with at least one selected model?
        - [ ] mutation rates of exons?
        - expected protein impact of selected event models?
        - in which terms are these genes enriched?
        - [ ] clear tumor suppressor, onco-exons or dual agent?
        - [X] classify interactions
            - heatmap/pca plot with classification, protein impact
            - interactions as if gene expression moderates splicing
            - how many events show significant interaction?
            - classify models into interaction analysis from Aviv's paper --> enrichments
            - classify models into moderation classes --> enrichments
            - protein impact of the interaction

        - [X] Examples
            - sort by LR test p-value
            - top in general
            - top by protein impact category
        - [X] spearman correlation with sample indices
        

4. [X] Cell lines embedded with splicing dependencies
    - use only selected models
    - overlay metadata and sample indices on the embeddings
    - clustering with leiden algorithm
        - cluster size
        - cluster properties using metadata and indices
        - ?? differential analysis across clusters (one vs all)
            - splicing dependencies --> different weak points
            - PSI --> different machineries
            - enrichments


5. [ ] Model validation/exploration
    - [ ] what's going on with events of interest:
        - make list of events know to affect cancer
        - which did we catch? Why not?
    - [X] CRISPR screen (Thomas 2020)
        - preprocess
        - consider selected event models
        - visualize event gene vs. true dependencies
        - include the confidence of our predictions


6. [ ] splicing dependencies in TCGA
    - [X] PT vs STN
        - differential PSI analysis (Mann Whitney U)
        - prioritize targets for every cancer type
            - pick only those for which there is an event model
            - re-calculate FDR (we only look at those)
            - select differentially spliced events
            - sort differentially spliced events by their median splicing dependency predicted in primary tumors (we can only be sure our models work in tumor samples)
        - interesting examples --> negative median splicing dependency and positive delta PSI (more inserted in PT)
            - plot median splicing dependency vs delta PSI with top 10
            - median splicing dependency for all cancers together
            - splicing dependencies for BRCA subtypes
            - ?? embeddings of splicing dependencies for all PT samples together cancers in selected models
            
    - [X] RESPONDER vs NON-RESPONDER
        - preprocess Moiso 2021 adding sampleIDs
        - differential analyses:
            - differential PSI and splicing dependency by cancer and drug
        - plot number of samples per cancer and treatment; move on only with those with a good number of samples
        - volcano plots --> select exons of interest --> gene set enrichment
        - plot deltaPSI vs deltaSplDep
        - any possible causes of different response?
        - [X] are these treatments relevant nowadays? --> talk to Fran
        - [X] build model with selected events
            - for those cancers and treatments with enough samples
            - models: RF, LDA
            - classify patients as responders and non-responders
            - compare this set of events with random sets of events
            - does the splicing profile have the same information?
        - [ ] consider if patients are smokers
        
    - [ ] Survival Analysis (Cox P.H.):
        - consider Patient Free Survival to associate to splicing dependency
        - [Neary 2021](https://www.nature.com/articles/s41598-021-84211-y)
    
    - [ ] response to antiPD1 in different cancers
        - collect data from different studies


7. [ ] causal inference with splicing dependencies
    - which splicing event may trigger cells to stop proliferating when a drug is applied?
    - datasets of cells with treatments (ideally in the splicing machinery)
        - Sophie/Irene
            - /no_backup/jvalcarcel/Sophie/Share/Miquel
            - one with E7820 and E7107
            - several with H3B8800 -/+ SF3B1 mutations including the one from Irene
            - one with spliceostatin
            - one with PRMT inhibitors (I did not yet look very much into it)
        - Interesting datasets
            - combine SSO with drug treatment
    
8. [X] Association of drug sensitivity to splicing dependencies
    - model following [Goncalves 2020](https://www.embopress.org/doi/full/10.15252/msb.20199405)
        - GDCS drug screens; PRISM drug screens??
        - model covatiates
            - binary variate indicating institute of origin of cell line
            - PC1 of the drug response dataset (correlates with cell lines growth rate)
            - growing conditions (adheret, suspension or semi-adherent)
            - gene fitness similarity considered as random effects in the model to account for potential sample structure
    $$ DrugIC50 = \beta_0 + \beta_1*SplDep +  \beta_2*PC1 + \beta_3*GrowCond + \epsilon_{SplDep_cov}$$
            - LR test to see if term $\beta_1*SplDep$ contributes significantly
            - FDR correction per drug tested
        - they used the limix package: https://github.com/EmanuelGoncalves/dtrace
        - for each drug we may get a list of splicing dependencies that are associated to its effects
    - There are drugs whose sensitivity correlates with splicing dependency of their targets but also to other splicing dependencies (more than expected)
    - could somehow splicing dependencies of other genes be mediators of the effect of the drug on cell viability?
    - [ ] check closeness of extreme drug-splicing dependency associations to the drug target using a PPI (HIPPIE, or STRING)
    - [ ] with drugs that have exactly the same targets, check the pattern of associations; maybe the drug has different structural targets (or off target effects). Check correlations of effects among drugs that share their targets.
    - [ ] therapy recommendation:
        - which is the drug with maximum sensitivity in the sample?
        - which is the drug strongly associated with the splicing dependency of a targetable event?
        