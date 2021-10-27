#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# EDA of gene dependencies regressed on event PSI and gene TPMs.
# 
# Outline
# -------
#  

require(tidyverse)
require(tidytext)
require(ggpubr)
require(cowplot)
require(ggrepel)
require(ComplexHeatmap)
require(gridExtra)
require(ggplotify)

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

# variables
THRESH_FDR = 0.05
THRESH_PVALUE = 0.05

# Development
# -----------
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# RESULTS_DIR = file.path(ROOT,'results','splicing_dependency_drugs')
# models_file = file.path(RESULTS_DIR,'files','models_drug_response-gdsc-EX.tsv.gz')
# drug_targets_file = file.path(RAW_DIR,'GDSC','screened_compunds_rel_8.2.csv')
# figs_dir = file.path(RESULTS_DIR,'figures','model_drug_screens')

##### FUNCTIONS #####
plot_associations = function(models, drug_targets){
    top_n = 25
    X = models
    
    plts = list()
    
    # what are the distributions of p-values and FDR?
    plts[['associations-lr_pvalues']] = X %>% 
        gghistogram(x='lr_pvalue', fill='lightblue', color=NA, bins=100) + 
        geom_vline(xintercept=median(X[['lr_pvalue']], na.rm=TRUE), 
                   linetype='dashed') +
        labs(x='LR Test p-value')
    
    plts[['associations-lr_fdr']] = X %>% 
        gghistogram(x='lr_padj', fill='darkred', color=NA, bins=100) + 
        geom_vline(xintercept=median(X[['lr_padj']], na.rm=TRUE), 
                   linetype='dashed') +
        labs(x='LR Test FDR')
    
    # how many drugs and events are significantly associated?
    plts[['associations-drug_counts']] = X %>% 
        filter(lr_padj < THRESH_FDR) %>% 
        count(drug_name) %>% 
        arrange(n) %>% 
        ggbarplot(x='drug_name', y='n', fill='#2a9d8f', color=NA) + 
        labs(x='Drug', y='No. Significant Associations')
    
    plts[['associations-top_drug_counts']] = X %>% 
        filter(lr_padj < THRESH_FDR) %>% 
        count(drug_name) %>% 
        arrange(n) %>% 
        tail(top_n) %>%
        ggbarplot(x='drug_name', y='n', fill='#2a9d8f', color=NA) + 
        labs(x='Drug', y='No. Significant Associations', 
             title=sprintf('Top %s', top_n)) +
        coord_flip()
    
    plts[['associations-spldep_counts']] = X %>% 
        filter(lr_padj < THRESH_FDR) %>% 
        count(event_gene) %>% 
        arrange(n) %>% 
        ggbarplot(x='event_gene', y='n', fill='orange', color=NA) + 
        labs(x='Event & Gene', y='No. Significant Associations')
    
    plts[['associations-top_spldep_counts']] = X %>% 
        filter(lr_padj < THRESH_FDR) %>% 
        count(event_gene) %>% 
        arrange(n) %>% 
        tail(top_n) %>%
        ggbarplot(x='event_gene', y='n', fill='orange', color=NA) + 
        labs(x='Event & Gene', y='No. Significant Associations', 
             title=sprintf('Top %s', top_n)) +
        coord_flip()
    
    # is there any drug that is significantly associated to the splicing
    # dependency of an event of a gene that it targets?
    found_targets = X %>% 
        filter(lr_padj < THRESH_FDR) %>% 
        mutate(drug_name=gsub("(.*),.*", "\\1",drug_name)) %>% 
        left_join(drug_targets, by=c('GENE'='TARGET','drug_name'='DRUG_NAME')) %>%
        drop_na(TARGET_PATHWAY) %>%
        mutate(is_target=TRUE)
    drugs_oi = found_targets %>% pull(drug_name) %>% unique()
    events_oi = found_targets %>% pull(EVENT) %>% unique()
    
    plts[['associations-top_found_targets']] = X %>% 
        filter(lr_padj<THRESH_FDR & drug_name%in%drugs_oi) %>% 
        group_by(drug_name) %>% 
        slice_max(order_by=abs(spldep_coefficient), n=15) %>% 
        left_join(distinct(found_targets[,c('drug_name','GENE','is_target')]), 
                  by=c('drug_name','GENE')) %>% 
        left_join(distinct(found_targets[,c('drug_name','TARGET_PATHWAY')]), 
                  by=c('drug_name')) %>% 
        mutate(is_target=replace_na(is_target,FALSE), 
               log10_pvalue=-log10(lr_pvalue), 
               drug_name_clean=sprintf('%s | %s',drug_name,TARGET_PATHWAY), 
               event_gene=reorder_within(event_gene, spldep_coefficient, drug_name)) %>% 
        ggbarplot(x='event_gene', y='spldep_coefficient', fill='is_target', 
                  palette='lancet', color=NA) + 
        facet_wrap(~drug_name_clean, scales='free_y', ncol=2) + 
        scale_x_reordered() + 
        labs(x='Event & Gene', y='Effect Size', fill='Is Drug Target', 
             title='Top 15 Significant Effect Sizes') + 
        coord_flip() +
        theme_pubr(border=TRUE)
    
    return(plts)
}


plot_examples = function(models, drug_targets){
    plts = list()
    # SRSF7 is important in LUAD treated with carboxi + paclitaxel
    #models %>% 
    
    # overlaps between targets of MASITINIB, PAZOPANIB, PONATINIB
    drugs_oi = c('MASITINIB','PAZOPANIB','PONATINIB')
    m = drug_targets %>% 
        filter(DRUG_NAME %in% drugs_oi) %>% 
        group_by(DRUG_NAME) %>%
        with(.,split(TARGET,DRUG_NAME)) %>% 
        list_to_matrix() %>% 
        make_comb_mat()
    plts[['examples-common_targets-upset']] = UpSet(
        m, comb_order = order(comb_size(m))
    ) %>%
        draw() %>%
        grid.grabExpr() %>%
        as.ggplot()
    
    # associations with HsaEX0034998_KRAS
    event_oi = 'HsaEX0034998_KRAS'
    ptls[['examples-kras']] = models %>% 
        filter(event_gene==event_oi & lr_padj<THRESH_FDR) %>%
        slice_max(abs(spldep_coefficient), n=15) %>%
        arrange(spldep_coefficient) %>%
        ggbarplot(x='drug_name', y='spldep_coefficient', 
                  fill='orange', color=NA) + 
        labs(x='Drug Name', y='Effect Size', title=event_oi) +
        coord_flip()
    
    drugs_oi = c('TOZASERTIB','BI-2536')
    ptls[['examples-kras_top_drugs']] = models %>% 
        filter(drug_name%in%drugs_oi & lr_padj<THRESH_FDR) %>%
        group_by(drug_name) %>%
        slice_max(abs(spldep_coefficient), n=15) %>%
        mutate(event_gene=reorder_within(event_gene, spldep_coefficient, drug_name)) %>%
        ggbarplot(x='event_gene', y='spldep_coefficient', 
                  fill='#2a9d8f', color=NA)+ 
        facet_wrap(~drug_name, scales='free_y', ncol=2) +
        scale_x_reordered() + 
        labs(x='Event & Gene', y='Effect Size', title='Top 15') +
        coord_flip()
    
    # associations with events in TCGA
    events_oi = c(
        'HsaEX0034998_KRAS',
        'HsaEX1036699_SFPQ',
        'HsaEX0050350_PRPF3',
        'HsaEX0060707_SNRNP70',
        'HsaEX0072698_ZFR',
        'HsaEX0044199_NUF2',
        'HsaEX0060960_SON',
        'HsaEX0037668_MAP4K4'
    )
    models %>% 
        filter(event_gene%in%events_oi & lr_padj<THRESH_FDR) %>%
        group_by(event_gene) %>%
        slice_max(abs(spldep_coefficient), n=15) %>%
        mutate(drug_name=reorder_within(drug_name, spldep_coefficient, event_gene)) %>%
        ggbarplot(x='drug_name', y='spldep_coefficient',
                  fill='orange', color=NA) + 
        facet_wrap(~event_gene, scales='free_y', ncol=2) +
        scale_x_reordered() + 
        labs(x='Drug Name', y='Effect Size', title='Top 15') +
        coord_flip()
    
    return(plts)
}


make_plots = function(models, drug_targets){
    plts = list(
        plot_associations(models, drug_targets)
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
                  base_width=width, base_height=height, dpi=dpi)
}


save_plots = function(plts, figs_dir){
    save_plt(plts, 'associations-lr_pvalues', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'associations-lr_fdr', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'associations-drug_counts', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'associations-top_drug_counts', '.pdf', figs_dir, width=8, height=5)
    save_plt(plts, 'associations-spldep_counts', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'associations-top_spldep_counts', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'associations-top_found_targets', '.pdf', figs_dir, width=10, height=10)

}


make_figdata = function(results_enrich){
    figdata = list(
        'gsea-hallmarks'= results_enrich[['hallmarks']],
    )
    return(figdata)
}


save_figdata = function(figdata, dir){
    lapply(names(figdata), function(x){
        filename = file.path(dir,'figdata',paste0(x,'.xlsx'))
        dir.create(dirname(filename), recursive=TRUE)
        write_xlsx(figdata[[x]], filename)
    })
}


main = function(){
    args = getParsedArgs()
    models_file = args$models_file
    drug_targets_file = args$drug_targets_file
    figs_dir = args$figs_dir
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    models = read_tsv(models_file) %>% 
        mutate(event_gene = paste0(EVENT,'_',GENE),
               event_type = gsub('Hsa','',gsub("[^a-zA-Z]", "",EVENT)))
    drug_targets = read_csv(drug_targets_file) %>% 
        mutate(DRUG_NAME=toupper(DRUG_NAME)) %>%
        dplyr::select(DRUG_NAME,TARGET,TARGET_PATHWAY) %>%
        separate_rows(TARGET) %>%
        distinct()
    
    # make plots
    plts = make_plots(models, drug_targets)

    # make figdata
    # figdata = make_figdata(results_enrich)
    
    # save
    save_plots(plts, figs_dir)
    # save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}