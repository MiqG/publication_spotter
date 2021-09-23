# Splicing dependency

![](https://fontmeme.com/permalink/210921/2be314f2a86b20b2bc77404529de028a.png)

## Outline
1. [X] Overview of data: 
    - CCLE gene dependencies
    - mRNA levels
    - splicing (maybe split splicing types?)
    - stemness index


2. [X] Use path analysis principles to decompose gene dependencies into splicing dependency

    $$ Gene Dependency = \beta_0 + \beta_1*PSI + \beta_2*TPM + \beta_3*PSI*TPM $$
    $$ Splicing Dependency = \beta_0 + \beta_1*PSI + \beta_3*PSI*TPM $$

    - Fit linear models using OLS
        - consider using GLS to correct for population structure as in https://github.com/kundajelab/coessentiality/blob/master/cancer_type_dependencies.py
        - then we would have population-unbiased coefficients and we would use them to make marginal predictions from the random effects model.
    - summarize them into a single model to get the empirical coefficients
    - sample 1000 times the coefficients to be able to make predictions with confidence


3. [X] Explore the models
    - p-value distributions
    - Z-score distributions
    - what are we losing? --> events with few observations? small variation?
    - significant events per gene --> possible limitation
    - gene set enrichment --> explore main sets in depth
        - GO BP
        - event protein impact
        - MSigDB: hallmarks and cancer signatures
    - examples of extreme coefficients --> biology?
    - [X] mutation frequency and mutation entropy of genes with significant splicing dependency --> we expect genes with many mutations and high mutation frequency across cell lines to not contain events associated to gene dependency.


4. [X] Cell lines embedded with splicing dependencies
    - overlay metadata and sample indices on the embeddings
    - clustering with leiden algorithm
        - cluster size
        - cluster properties using metadata and indices
        - ??differential splicing dependencies across clusters (one vs all)
            - enrichments


5. [ ] Model validation
    - [ ] what's going on with genes of interest: 
        - cancer genes (get events)
        - known drug targets
        - splicing factors
    - [X] CRISPR screen (Thomas 2020)
        - preprocess
        - re-read to make sure we are selecting exons correctly (those very deviant)
        - include the confidence of our predictions

    
    
6. [X] Explore model interactions of splicing and gene expression
    - interactions as if gene expression moderates splicing
    - how many events show significant interaction?
    - classify models into interaction analysis from Aviv's paper --> enrichments
    - classify models into moderation classes --> enrichments
    - protein impact of the interaction
        

7. [ ] splicing dependencies in TCGA
    - PT vs STN
        - differential PSI
        - differential splicing dependency
        - volcano plots --> select exons of interest --> gene set enrichment
        - prioritize targets for every cancer type --> deltaPSI vs deltaSplDep
    - RESPONDER vs NON-RESPONDER
        - preprocess Moiso 2021 adding sampleIDs
        - differential PSI
        - differential splicing dependency
        - volcano plots --> select exons of interest --> gene set enrichment
        - suggest possible causes of different response
    - BRCA subtypes different splicing dependency
       
       
8. [ ] causal inference with splicing dependencies
    - which splicing event may trigger cells to stop growing when a drug is applied?
    - datasets of cells with treatments (ideally in the splicing machinery)
    

100. [ ] what do exons known to trigger NMD do?
    - overlaps with lists of poison exons
    - dataset with inhibited NMD
    
    