#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# EDA of gene dependencies regressed on event PSI and gene TPMs.
# 
# possible interactions (from DOI: https://doi.org/10.1016/j.cell.2016.11.038)
# 
# Outline
# -------
# - pvalue distributions of linear models? p-values vs n. observations?
# - scoring of linear models? spearman vs pearson?
# - which events/genes/interaction/intercepts have extreme z-scores?
# - to which biological processes are they related?
# - Interesting examples? Cancer fitness, known drug targets, splicing factors
# - project cell lines based on real fitness, fitness predictions vs PSI (impute) and gene expression (make a script)
#  

require(tidyverse)
require(ggrepel)
require(readxl)
require(ggpubr)
require(latex2exp)
require(ggrepel)
require(ComplexHeatmap)
require(cowplot)
require(scattermore)

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

MIN_OBS = 50
THRESH_ZSCORE = 1.96


# Development
# -----------
psi_file = file.path(ROOT,'data','raw','articles','Thomas2020','vast_out','PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz')
genexpr_file = file.path(ROOT,'data','raw','articles','Thomas2020','vast_out','TPM-hg38-2.tab.gz')
dependency_file = file.path(ROOT,'results','model_splicing_dependency','files','splicing_dependencies.tsv.gz')
possible_interactions_file = file.path(ROOT,'support','possible_pairwise_interaction_categories.tsv')
figs_dir = file.path(ROOT,'results','model_splicing_dependency','figures','eda')


##### FUNCTIONS #####
get_interaction_categories = function(dependency, possible_interactions){
    # study interaction between event inclusion and gene expression
    possible_interactions = possible_interactions %>% 
        mutate(combined=paste0(beta_a,beta_b,beta_ab))
    
    ## interactions
    zscores = dependency %>% 
        column_to_rownames("event_gene") %>% 
        dplyr::select(paste0(c('event','gene','interaction'),"_zscore"))
    intcats = (abs(zscores) > THRESH_ZSCORE) * sign(zscores)
    intcats[is.na(intcats)] = 0
    intcats = intcats %>% 
        rownames_to_column("event_gene") %>%
        mutate(combined = paste0(event_zscore,gene_zscore,interaction_zscore)) %>%
        left_join(possible_interactions[,c('combined','category')],by='combined') %>%
        dplyr::select(event_gene,category) %>%
        rename(interaction_category=category) %>% 
        mutate(interaction_category=factor(
            interaction_category, 
            levels = c('buffering','synergistic','dominant','additive','no_regulation')
        ))
    
    return(intcats)
}


plot_interactions = function(){
    
    X = dependency
    
    plts = list()
    
    # across event types, what's the frequency of each type of psi-tpm interaction?
    plts[['interactions-counts_overview']] = X %>% 
        count(event_type, interaction_category) %>%
        ggbarplot(x='interaction_category', y='n', facet.by = 'event_type', 
                  label=TRUE, fill='interaction_category', color=NA, palette='Set2') + 
        guides(fill=FALSE) + 
        yscale('log10', .format=TRUE) + 
        labs(x='Interaction Category', y='Counts') +
        theme_pubr(x.text.angle = 70, border=TRUE)
    
    plts[['interactions-counts_max_zscore']] = X %>% 
        group_by(event_type, interaction_category, max_zscore) %>%
        summarise(n = n()) %>% mutate(freq = 100*(n / sum(n))) %>%
        ggbarplot(x='interaction_category', y='freq', facet.by = 'event_type', 
                  fill='max_zscore', color=NA, palette='lancet') + 
        labs(x='Interaction Category', y='%', fill='Max. Z-score') +
        theme_pubr(x.text.angle = 70, border=TRUE)
    
    # across event types and additive interactions, which term has a higher zscore
    # more frequently?
    effect_sizes = X %>% 
        filter(interaction_category!='no_regulation') %>% 
        replace_na(list(event_zscore = 0, gene_zscore = 0)) %>%
        mutate(
            zscore_diff = event_zscore - gene_zscore,
            zscore_diffsign = ifelse(sign(zscore_diff)>0, 'Splicing', 'Expression') 
        )
    plts[['interactions-effect_size_diffs']] = effect_sizes %>%
        ggviolin(x='interaction_category', y='zscore_diff', 
                 fill = 'zscore_diffsign', color=NA, palette='npg') + 
        geom_boxplot(aes(fill=zscore_diffsign), 
                     position = position_dodge(0.8), 
                     width=0.3, outlier.size = 0.5) +
        facet_wrap(~event_type) +
        labs(x='Interaction Category', 
             y=TeX('$Splicing_{Z-score} - Expression_{Z-score}$'),
             fill='Larger Effect Size') +
        theme_pubr(x.text.angle = 70, border=TRUE)
    
    plts[['interactions-effect_size_counts']] = effect_sizes %>%
        count(event_type, interaction_category, zscore_diffsign) %>%
        ggtexttable()
        
    return(plts)
}


plot_zscores = function(dependency){
    X = dependency 
    
    plts = list()
    plts[['zscores-event_vs_gene']] = X %>% 
        ggplot(aes(x=event_zscore, y=gene_zscore)) +
        geom_hex(bins=200) +
        facet_wrap(event_type) +
        stat_cor(method='spearman') +
        theme_pubr(border=TRUE)

    plts[['zscores-event_vs_interaction']] = X %>% 
        ggplot(aes(x=event_zscore, y=interaction_zscore)) +
        geom_hex(bins=200) +
        facet_wrap(event_type) +
        stat_cor(method='spearman') +
        theme_pubr(border=TRUE)

    plts[['zscores-gene_vs_interaction']] = X %>%
        ggplot(aes(x=gene_zscore, y=interaction_zscore)) +
        geom_hex(bins=200) +
        facet_wrap(event_type) +
        stat_cor(method='spearman') +
        theme_pubr(border=TRUE)
    
    plts[['zscores-event']] = X %>%
        ggviolin(x='event_type', y='event_zscore', fill='event_type', color=NA) +
        geom_boxplot() +
        guides(fill=FALSE)
    
    plts[['zscores-gene']] = X %>%
        group_by(event_type,GENE) %>%
        summarize(median_gene_zscore = median(gene_zscore, na.rm=TRUE)) %>%
        ggviolin(x='event_type', y='median_gene_zscore', 
                 fill='event_type', color=NA) +
        geom_boxplot() +
        guides(fill=FALSE)
    
    plts[['zscores-interaction']] = X %>%
        ggviolin(x='event_type', y='interaction_zscore', fill='event_type', color=NA) +
        geom_boxplot() +
        guides(fill=FALSE)
    
    # through intercepts, we see the baseline essentiality of the gene
    plts[['zscores-intercept']] = X %>%
        group_by(event_type,GENE) %>%
        summarize(median_intercept_zscore = median(intercept_zscore, na.rm=TRUE)) %>%
        ggviolin(x='event_type', y='median_intercept_zscore', 
                 fill='event_type', color=NA) +
        geom_boxplot() +
        guides(fill=FALSE)
        

    
}


make_plots = function(){
    plts = list(

    )
    plts = do.call(c,plts)
    return(plts)
}


save_plt = function(plts, plt_name, extension='.pdf', 
                      directory='', dpi=350, 
                      width = par("din")[1], height = par("din")[2]){
        filename = file.path(directory,paste0(plt_name,extension))
        save_plot(filename, 
                  plts[[plt_name]], 
                  width=width, height=height, dpi=dpi)
}


save_plots = function(plts, figs_dir){
    save_plt(plts, 'event_dependency-volcano', '.png', figs_dir, width=20, height=10)    
}


main = function(){
    args = getParsedArgs()
    dependency_file = args$dependency_file
    figs_dir = args$figs_dir
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    possible_interactions = read_tsv(possible_interactions_file)
    dependency = read_tsv(dependency_file) %>%
        filter(n_obs>MIN_OBS) %>%
        mutate(event_gene = paste0(EVENT,'_',GENE),
               event_type = gsub('Hsa','',gsub("[^a-zA-Z]", "",EVENT)))
    
    # add interaction categories to dependency
    intcats = get_interaction_categories(dependency, possible_interactions)
    dependency = dependency %>% left_join(intcats,by='event_gene')
    
    # add who has the largest effect size?
    effect_size = dependency %>% 
        column_to_rownames('event_gene') %>% 
        dplyr::select(paste0(c('event','gene','interaction'),'_zscore'))
    effect_size[is.na(effect_size)] = 0
    effect_size[['max_zscore']] = c('event','gene','interaction')[apply(abs(effect_size),1,which.max)]
    effect_size = effect_size %>% dplyr::select(max_zscore) %>% rownames_to_column('event_gene')
    dependency = dependency %>% left_join(effect_size, by='event_gene')
    
    # are exons with significant interaction terms associated to NMD?
    
    # plot
    plts = make_plots()

    # save
    save_plots(plts, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}