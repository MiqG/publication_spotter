# Splicing dependency

## Outline
1. Overview of data: 
    - CCLE gene dependencies
    - mRNA levels
    - splicing (maybe split splicing types?)
    - stemness index

2. Use path analysis principles to decompose gene dependencies into splicing dependency

3. Explore the models
    - p-value distributions
    - scoring
    - examples of extreme coefficient Z-scores --> biology
    - what's going on with known examples: cancer genes, drug targets, splicing factors
    - project cell lines based on splicing dependency predictions (using splicing and interaction term)
    
4. Explore moderation of splicing and gene expression
    - interaction analysis from Aviv's paper
    - Which models show moderation? Who is the moderator?
    - NMD:
        - overlaps with lists of poison exons
        - dataset with inhibited NMD
        
5. Validate models as predictors of splicing dependency
    - our significant models can predict exon CRISPR screen fitness scores
    - our significant models belong to genes with low rate of detrimental mutations

