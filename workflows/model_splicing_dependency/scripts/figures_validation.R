#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# EDA on event fitness from 
#   - Gonatopoulos-Pournatzis (2020) (DOI: https://doi.org/10.1038/s41587-020-0437-z) in RPE1 cells (retina-derived cell line)
# 
# Outline
# -------
# - prediction of fitness effect of exon CRISPR screen Thomas 2020: if we run the model with less data, do we get the same predictions? (How robust are the predictions?)
# - cancer associated genes tend to rank at the top
# - tested targets
# - HsaEX0044216 (NUMB) in certain cancer cell lines


require(tidyverse)
require(ggrepel)
require(readxl)
require(ggpubr)
require(latex2exp)
require(ggrepel)
require(ComplexHeatmap)
require(limma)

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

MIN_OBS = 50

# variables
PSI_COLS = c('EVENT','CL_HeLa')

# Development
# -----------
psi_file = file.path(ROOT,'data','raw','articles','Thomas2020','vast_out','PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz')

genexpr_file = file.path(ROOT,'data','raw','articles','Thomas2020','vast_out','TPM-hg38-2.tab.gz')

event_annotation_file = file.path(ROOT,'data','raw','VastDB','EVENT_INFO-hg38_noseqs.tsv')

thomas_exon_crispr_screen_file = file.path(ROOT,'data','raw','articles','Thomas2020','crispr_screen.xlsx')
thomas_event_mapping_file = file.path(ROOT,'data','raw','articles','Thomas2020','event_mapping_vastdb.tsv')

dependency_file = file.path(ROOT,'results','model_splicing_dependency','files','splicing_dependencies.tsv.gz')

#event_properties_ccle_file = file.path(ROOT,'results','exon_associations_ignoring_confounders','CCLE','files','exon_sample_properties_correlation-spearman-pt.tsv')

#event_properties_tcga_file = file.path(ROOT,'results','pancancer_exon_properties','files','tcga-analysis_summaries','merged.tsv.gz')

cancer_related_as_file = file.path(ROOT,'support','cancer_related_as_isoforms.tsv')

figs_dir = file.path(ROOT,'results','model_splicing_dependency','figures','validation')


##### FUNCTIONS #####
prep_thomas = function(X){
    cols = grep('\\|', colnames(X), value = TRUE)
    conditions = do.call(rbind,strsplit(cols, '\\|'))
    stats = unique(conditions[,1])
    #
    #colnames(conditions) = c('cell_line','comparison','replicate')
    cond = conditions
    screen = lapply(conditions, function(cond){
        screen = X %>% dplyr::select(contains(cond))
        colnames(screen) = stats
        cond = unlist(strsplit(cond, '\\.'))
        screen = screen %>% 
            mutate(cell_line=cond[1], 
                   comparison=cond[2], 
                   replicate=cond[3])
        cols_oi = c('event','EVENT','target_type','geneName')
        screen = cbind(X[,cols_oi],screen)
        return(screen)
    })
    screen = do.call(rbind,screen)
    return(screen)
}


load_profiles = function(psi_file, genexpr_file, event_annotation_file){
    psi_ccle = read_tsv(psi_file) %>% 
        dplyr::rename(HeLa=SRR7946515_1, PC9=SRR7946516_1) %>% 
        pivot_longer(cols=c('HeLa','PC9'), names_to='cell_line', 
                     values_to='event_psi', values_drop_na=TRUE)
        
    genexpr_ccle = read_tsv(genexpr_file) %>% 
        dplyr::select(-one_of('ID')) %>%
        dplyr::rename(GENE=NAME, HeLa=SRR7946515_1, PC9=SRR7946516_1) %>%
            pivot_longer(cols=c('HeLa','PC9'), names_to='cell_line', 
                         values_to='gene_tpm', values_drop_na=TRUE) %>%
            mutate(gene_tpm=log2(gene_tpm+1))
    
    annotation = read_tsv(event_annotation_file)
    
    profiles = psi_ccle %>% 
        left_join(annotation, by='EVENT') %>% 
        left_join(genexpr_ccle, by=c('GENE','cell_line'))
    
    return(profiles)
}


load_crispr_screens = function(
    thomas_exon_crispr_screen_file, 
    thomas_event_mapping_file,
    profiles
    ){
    
    # load
    thomas_crispr_screen = read_excel(thomas_exon_crispr_screen_file, skip=1, sheet=5)
    thomas_mapping = read_tsv(thomas_event_mapping_file)
    
    # prepare Thomas
    thomas_crispr_screen = thomas_crispr_screen %>% 
        left_join(thomas_mapping, by=c('event','target_type')) %>%
        filter(!is.na(EVENT))
    thomas_crispr_screen = prep_thomas(thomas_crispr_screen)
    thomas_crispr_screen = thomas_crispr_screen %>% mutate(
        fitness_score = log2(fc_norm_mean),
        study = 'Thomas (2020)',
        pvalue = pval,
        hit = fdr < 0.05
    )
    
    # add info from reads
    thomas_crispr_screen = thomas_crispr_screen %>% 
        left_join(profiles, by=c('EVENT','cell_line'))
    
    # prepare main dataframe
    cols_oi = c('EVENT','GENE','fitness_score','pvalue','fdr','cell_line',
                'comparison','replicate','hit','target_type','study',
                'event_psi','gene_tpm')
    crispr_screen = thomas_crispr_screen[,cols_oi] %>% 
      filter(!is.na(fitness_score)) %>% 
      mutate(reference=sprintf('%s - %s',cell_line,comparison)) %>%
      distinct(EVENT,fitness_score,study,hit,target_type, .keep_all=TRUE)
    
    return(crispr_screen)
}


plot_dependency = function(dependency, cancer_related_as){
    X = dependency
    genes_oi = unique(cancer_related_as[['gene_symbol']])
    
    plts = list()
#     plts[['dependency-volcano']] = X %>% 
#         ggscatter(x='event_zscore', y='log10_fdr', size=1, facet.by='event_type', alpha=0.8) + 
# #         geom_text_repel(aes(label=event_gene), 
# #                         X %>% filter(abs(event_zscore)>5)) +
#         ggtitle(sprintf('Association PSI with dependency score | n=%s',nrow(X))) + 
#         labs(x='Z-score', y=TeX('$-log_{10}(p-value)$'))
    
    
#     med_pvalue = X %>% group_by(event_type) %>% summarise(med=median(event_pvalue,na.rm=TRUE))
#     plts[['dependency-event_pvalue']] = X %>% 
#         gghistogram(x='event_pvalue', facet.by='event_type',fill='grey') + 
#         geom_vline(aes(xintercept=med), med_pvalue, linetype='dashed') + 
#         geom_text(aes(x=med-0.07, y=6000, label=round(med,2)), med_pvalue)
        
    
    plts[['cancer_related_genes']] = X %>% 
        group_by(GENE) %>% 
        summarize(coef = event_coefficient[which.max(abs(event_coefficient))]) %>% 
        ungroup() %>%
        mutate(cancer_related = factor(GENE %in% genes_oi, levels=c(TRUE,FALSE)), 
               rank = rank(-abs(coef))) %>%
        ggviolin(y='rank',x='cancer_related',
                 fill='cancer_related',color=NA,
                 palette='uchicago', trim=TRUE) +
        geom_boxplot(width=0.1) +
        ggtitle(sprintf('n. cancer genes = %s',length(genes_oi))) + 
        stat_compare_means(label.y = 100) +
        scale_y_reverse() +
        labs(x='Is Cancer Related', y='Rank') +
        guides(fill=FALSE)
    
    
    return(plts)
}


plot_crispr_validation = function(dependency, crispr_screen){
    X = crispr_screen %>%
        left_join(dependency, by='EVENT') %>% 
        mutate(is_significant = factor(event_pvalue < 0.01, c(TRUE,FALSE))) %>%
        filter(!is.na(fitness_prediction)) %>%
        filter(comparison %in% c('0d-vs-14d','0d-vs-8d'))
    
    plts = list() 
    
    # is the distribution of dependency scores the same?
    plts[['']] = 
    dependency %>% filter(event_pvalue < 0.01) %>% 
        mutate(fitness_score=fitness_prediction, score_type='predicted') %>%
        dplyr::select(fitness_score, score_type, cell_line) %>% 
        bind_rows(
            crispr_screen %>% 
                mutate(score_type='crispr_screen') %>%
                dplyr::select(fitness_score, score_type, cell_line)) %>%
        drop_na(cell_line) %>%
        ggboxplot(x='cell_line', y='fitness_score', fill='score_type') + stat_compare_means(aes(group=score_type), method='wilcox.test')
        ggdensity(x='fitness_score', fill='score_type', facet.by = c('cell_line')) %>%
        
    
    # how spearman correlation of fitness with significance threshold
    threshs = c(0.1,0.05,0.025,0.01,0.005,0.0025,0.001)
    groups = X %>% pull(reference) %>% unique()
    res = lapply(threshs, function(thresh){
        res = lapply(groups, function(group){
            tmp = X %>% filter(event_pvalue<thresh & reference==group)
            res = cor.test(tmp$fitness_prediction, tmp$fitness_score, method='spearman')
            df = data.frame(
                correlation=res$estimate,
                pvalue=-log10(res$p.value),
                group=group,
                thresh=thresh,
                n=nrow(tmp)
            )
            return(df)
        })
        res = do.call(rbind,res)
    })
    res = do.call(rbind, res)
    
    plts[['crispr-correlations_vs_thresholds_event_pvalue']] = res %>% 
        mutate(thresh=as.factor(thresh)) %>% 
        ggscatter(x='thresh', y='correlation', 
                  size = 'n', shape = 'group', color='pvalue') +
        gradient_color(c("red", "yellow","blue")) +
        theme(aspect.ratio = 1) +
        labs(x='Significance Threshold', y='Spearman Corr.', 
             shape='Experiment', size='n', color='-log10(p-value)')
    
    
    # validate predictions with event screens
    plts[['crispr-fitness_true_vs_pred']] = X %>% 
        ggscatter(x='fitness_prediction', y='fitness_score', 
                  alpha=0.5, color='is_significant', 
                  facet.by=c('comparison','cell_line'), 
                  add = 'reg.line', conf.int=TRUE, 
                  add.params = list('linetype'='dashed'), palette='uchicago') + 
        stat_cor(aes(color=is_significant), label.y=c(-3,-3.5), method='spearman') + 
#         geom_text_repel(aes(x=fitness_prediction, y=fitness_score, label=event_gene), 
#                        X %>% filter(as.logical(is_significant))) +
        labs(x='Pred. Fitness Score', y='Fitness Score', color='Is Significant') +
        guides(fill=FALSE)
    
    return(plts)
}


make_plots = function(dependency, cancer_related_as, crispr_screen){
    plts = list(
        plot_dependency(dependency, cancer_related_as),
        plot_crispr_validation(dependency, crispr_screen)
    )
    plts = do.call(c,plts)
    return(plts)
}


save_plot = function(plt, plt_name, extension='.pdf', 
                      directory='', dpi=300, 
                      width = par("din")[1], height = par("din")[2]){
        filename = file.path(directory,paste0(plt_name,extension))
        ggsave(filename, 
               plt, 
               width=width, height=height, dpi=dpi)
}


save_plots = function(plts, figs_dir){
    # save manually
    save_plot(plts[['cancer_related_genes']], 'cancer_related_genes', '.pdf', figs_dir, 350, width=8, height=8)
    save_plot(plts[['crispr-correlations_vs_thresholds_event_pvalue']]+theme(aspect.ratio = 1), 'crispr-correlations_vs_thresholds_event_pvalue', '.pdf', figs_dir, 350, width=12, height=6)
    save_plot(plts[['crispr-fitness_true_vs_pred']], 'crispr-fitness_true_vs_pred', '.png', figs_dir, 350, width=15, height=15)    
}


main = function(){
    args = getParsedArgs()
    thomas_crispr_screen = args$thomas_crispr_screen
    thomas_event_mapping_file = args$thomas_event_mapping_file
    dependency_file = args$dependency_file
    cancer_related_as = args$cancer_related_as_file
    figs_dir = args$figs_dir
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    dependency = read_tsv(dependency_file) %>%
        filter(n_obs>MIN_OBS) %>%
        mutate(
            event_zscore=event_coefficient/event_stderr,
            event_gene = paste0(EVENT,'_',GENE),
            log10_pvalue = (-1)*log10(event_pvalue),
            event_type = ifelse(grepl('EX',EVENT),'EX','ALT')
        ) %>%
        group_by(event_type) %>%
        mutate(event_fdr = p.adjust(event_pvalue, method='fdr'),
               log10_fdr = -log10(event_fdr)) %>%
        ungroup()
        
    cancer_related_as = read_tsv(cancer_related_as_file)
    profiles = load_profiles(psi_file, genexpr_file, event_annotation_file)
    crispr_screen = load_crispr_screens(
        thomas_exon_crispr_screen_file, 
        thomas_event_mapping_file,
        profiles
    )
    
    dependency = dependency %>%
        left_join(profiles, by=c('EVENT','GENE')) %>%
        mutate(norm_psi = (event_psi - event_mean) / event_std,
               norm_genexpr = (gene_tpm - gene_mean) / gene_std,
               fitness_prediction = event_coefficient * norm_psi + intercept_coefficient)
    
    # plot
    plts = make_plots(dependency, cancer_related_as, crispr_screen)

    # save
    save_plots(plts, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}