# Splicing dependency

## Outline
1. [X] Overview of data: 
    - CCLE gene dependencies
    - mRNA levels
    - splicing (maybe split splicing types?)
    - stemness index

2. [X] Use path analysis principles to decompose gene dependencies into splicing dependency

3. [X] Explore the models
    - p-value distributions
    - scoring
    - Z-score distributions
    - what are we losing? --> events with few observations
    - significant events per gene
    - gene set enrichment
    - examples of extreme coefficient Z-scores --> biology
    - what's going on with genes of interest: cancer genes, drug targets, splicing factors
    - project cell lines based on splicing dependency predictions (using splicing and interaction term)
    
4. [ ] Explore moderation of splicing and gene expression
    - interaction analysis from Aviv's paper
    - Which events models show moderation? Who is the moderator?
    - NMD:
        - overlaps with lists of poison exons
        - dataset with inhibited NMD
        
5. [ ] Validate models as predictors of splicing dependency
    - our significant models can predict exon CRISPR screen fitness scores
    - our significant models belong to genes with low rate of detrimental mutations

