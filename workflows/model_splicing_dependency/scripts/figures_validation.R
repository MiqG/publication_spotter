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
# 1. Load 
#     - exon CRISPR screen fitness table
#         - Gonatopoulos-Pournatzis (2020)
#         - Thomas (2020)
#     - event properties
# 2. How does every property correlate with fitness across cancers?
# 3. Get the detailed correlations.

require(tidyverse)
require(ggrepel)
require(readxl)
require(ggpubr)
require(latex2exp)
require(ggrepel)
require(pheatmap)

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

MIN_OBS = 100


# Development
# -----------
psi_ccle_file = file.path(ROOT,'data','raw','VastDB','cell_lines-Hs2136-VastDB-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab')

genexpr_ccle_file = file.path(ROOT,'data','raw','DepMap','achilles_ccle','CCLE_expression_transposed.tsv.gz')

event_annotation_file = file.path(ROOT,'data','references','vastdb_events_annotation.tsv')

gonatopoulos_exon_crispr_screen_file = file.path(ROOT,'data','raw','articles','Gonatopoulos-Pournatzis2020','exon_targeting_library_scores.xlsx')

thomas_exon_crispr_screen_file = file.path(ROOT,'data','raw','articles','Thomas2020','crispr_screen.xlsx')
thomas_event_mapping_file = file.path(ROOT,'data','raw','articles','Thomas2020','event_mapping_vastdb.tsv')

event_dependency_file = file.path(ROOT,'results','exon_associations_ignoring_confounders','CCLE','files','exon_association_with_rnai_fitness.tsv.gz')

event_properties_ccle_file = file.path(ROOT,'results','exon_associations_ignoring_confounders','CCLE','files','exon_sample_properties_correlation-spearman-pt.tsv')

event_properties_tcga_file = file.path(ROOT,'results','pancancer_exon_properties','files','tcga-analysis_summaries','merged.tsv.gz')

cancer_related_as_file = file.path(ROOT,'support','cancer_related_as_isoforms.tsv')

figs_dir = file.path(ROOT,'results','pancancer_exon_properties','figures','event_dependency')


##### FUNCTIONS #####
prep_thomas = function(X){
    cols = grep('\\|', colnames(X), value = TRUE)
    conditions = do.call(rbind,strsplit(cols, '\\|'))
    stats = unique(conditions[,1])
    conditions = unique(conditions[,2])
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


load_crispr_screens = function(
    gonatopoulos_exon_crispr_screen_file, 
    thomas_exon_crispr_screen_file, 
    thomas_event_mapping_file,
    psi_ccle_file,
    genexpr_ccle_file,
    event_annotation_file
    )
{
    # load
    psi_ccle = read_tsv(psi_ccle_file)
    genexpr_ccle = read_tsv(genexpr_ccle_file) %>% 
        dplyr::rename(GENE=X1) %>%
        mutate(GENE=gsub(' .*', '', GENE))
    annotation = read_tsv(event_annotation_file)
    gonatopoulos_crispr_screen = read_excel(gonatopoulos_exon_crispr_screen_file, skip=1)
    thomas_crispr_screen = read_excel(thomas_exon_crispr_screen_file, skip=1, sheet=5)
    thomas_mapping = read_tsv(thomas_event_mapping_file)
    
    # prepare Gonatopoulos-Pournatzis
    gonatopoulos_crispr_screen = gonatopoulos_crispr_screen %>% mutate(
        fitness_score = as.numeric(ess.RPE1.log2FC),
        pvalue = NA,
        fdr = NA,
        EVENT = Event,
        study = 'Gonatopoulos-Pournatzis (2020)',
        target_type = type.RPE1,
        cell_line = 'RPE1',
        comparison = 'before-vs-post',
        replicate = NA,
        event_psi = as.numeric(PSI.RPE1),
        gene_tpm = NA
    )
    
    # prepare Thomas
    thomas_crispr_screen = merge(
        thomas_crispr_screen, thomas_mapping, by.x='event', by.y='original'
    ) %>%  filter(!is.na(EVENT))
    thomas_crispr_screen = prep_thomas(thomas_crispr_screen)
    thomas_crispr_screen = thomas_crispr_screen %>% mutate(
        fitness_score = log2(fc_norm_mean),
        study = 'Thomas (2020)',
        pvalue = pval,
        hit = fdr < 0.05
    ) %>% 
    filter(target_type == 'conservedPE') %>% 
    left_join(psi_ccle[,c('EVENT','CL_HeLa')] %>% 
              mutate(cell_line='HeLa'), by=c('EVENT','cell_line')) %>% 
    dplyr::rename(event_psi = CL_HeLa) %>%
    left_join(annotation[,c('EVENT','GENE')], by='EVENT') %>% 
    left_join(genexpr_ccle[,c('GENE','ACH-001086')] %>% 
              mutate(cell_line='HeLa'), by=c('GENE','cell_line')) %>% 
    dplyr::rename(gene_tpm = `ACH-001086`)
    
    # prepare main dataframe
    cols_oi = c('EVENT','fitness_score','pvalue','fdr','cell_line',
                'comparison','replicate','hit','target_type','study',
                'event_psi','gene_tpm')
    crispr_screen = rbind(gonatopoulos_crispr_screen[,cols_oi],
                          thomas_crispr_screen[,cols_oi]) %>% 
      filter(!is.na(fitness_score)) %>% 
      mutate(reference=sprintf('%s - %s',cell_line,comparison)) %>%
      distinct(EVENT,fitness_score,study,hit,target_type, .keep_all=TRUE)
    
    
    
    
    return(crispr_screen)
}


plot_event_dependency = function(event_dependency, cancer_related_as){
    X = event_dependency
    genes_oi = unique(cancer_related_as[['gene_symbol']])
    
    plts = list()
    plts[['event_dependency-volcano']] = X %>% 
        ggscatter(x='event_zscore', y='log10_fdr', size=1, facet.by='event_type', alpha=0.8) + 
        geom_text_repel(aes(label=event_gene), 
                        X %>% filter(abs(event_zscore)>5)) +
        ggtitle(sprintf('Association PSI with dependency score | n=%s',nrow(X))) + 
        labs(x='Z-score', y=TeX('$-log_{10}(p-value)$'))
    
    
    med_pvalue = X %>% group_by(event_type) %>% summarise(med=median(event_pvalue,na.rm=TRUE))
    plts[['event_dependency-event_pvalue']] = X %>% 
        gghistogram(x='event_pvalue', facet.by='event_type',fill='grey') + 
        geom_vline(aes(xintercept=med), med_pvalue, linetype='dashed') + 
        geom_text(aes(x=med-0.07, y=6000, label=round(med,2)), med_pvalue)
        
    
    plts[['event_dependency-cancer_related_genes']] = X %>% 
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


plot_crispr_validation = function(event_dependency, crispr_screen){
    X = merge(event_dependency, crispr_screen, by='EVENT') %>% 
        mutate(norm_psi = (event_psi - event_mean) / event_std,
               norm_genexpr = (gene_tpm - gene_mean) / gene_std,
               fitness_prediction_ = event_coefficient * norm_psi + gene_coefficient * norm_genexpr + intercept_coefficient,
               fitness_prediction = event_coefficient * norm_psi + intercept_coefficient,
               is_significant = factor(event_pvalue < 0.001, c(TRUE,FALSE))) %>%
        filter(!is.na(fitness_prediction) & cell_line=='HeLa') 
    
    plts = list() 
    # how spearman correlation of fitness with significance threshold
    threshs = c(0.1,0.05,0.01,0.005,0.001)
    groups = X %>% pull(reference) %>% unique()
    res = lapply(threshs, function(thresh){
        res = lapply(groups, function(group){
            tmp = X %>% filter(event_pvalue<thresh & reference==group)
            res = cor.test(tmp$fitness_prediction, tmp$fitness_score, method='spearman')
            df = data.frame(
                correlation=res$estimate,
                pvalue=res$p.value,
                group=group,
                thresh=thresh,
                n=nrow(tmp)
            )
            return(df)
        })
        res = do.call(rbind,res)
    })
    res = do.call(rbind, res)
    
    plts[['event_dependency-validation_correlations']] = res %>% 
        mutate(thresh=as.factor(thresh)) %>% 
        ggscatter(x='thresh', y='correlation', 
                  size = 'n', shape = 'group', color='pvalue') +
        gradient_color(c("red", "white", "blue")) +
        theme(aspect.ratio = 1) +
        labs(x='Significance Threshold', y='Spearman Corr.', 
             shape='Experiment', size='n', color='p-value')
    
    
    # validate predictions with event screens
    plts[['event_dependency-validation_scatter']] = X %>% 
        ggscatter(x='fitness_prediction', y='fitness_score', 
                  alpha=0.5, facet.by='reference', color='is_significant', 
                  add = 'reg.line', conf.int=TRUE, 
                  add.params = list('linetype'='dashed'), palette='uchicago') + 
        stat_cor(aes(color=is_significant), label.y = c(-1.6,-1.5), method='spearman') + 
        geom_text_repel(aes(x=fitness_prediction, y=fitness_score, label=event_gene), 
                       X %>% filter(as.logical(is_significant))) +
        labs(x='Pred. Fitness Score', y='Fitness Score', color='Is Significant') +
        guides(fill=FALSE)
    
    return(plts)
}


scale2 = function(x){ 
    x = as.numeric(scale(x))
    x = (x - min(x, na.rm=TRUE)) / (max(x, na.rm=TRUE) - min(x, na.rm=TRUE)) 
    return(x)
}


plot_diffsplicing_dependency = function(event_dependency, event_properties_tcga){
    palette = get_palette('Paired',length(unique(event_properties_tcga[['cancer_type']])))
    n_cancer_types = event_properties_tcga %>% 
        filter(is_diff_spliced) %>% 
        distinct(cancer_type) %>% 
        nrow()
    event_freq = event_properties_tcga %>% 
        filter(is_diff_spliced) %>% 
        distinct(EVENT,cancer_type) %>% 
        drop_na() %>% 
        group_by(EVENT) %>%
        summarize(event_freq=n()) %>% 
        ungroup() %>%
        mutate(rank_diff_spliced = rank(-event_freq))
    events_oi = event_freq %>% 
        filter(event_freq>round(0.5*n_cancer_types)) %>% 
        pull(EVENT)
    
    df = merge(event_dependency, event_freq, by=c('EVENT')) %>% 
        group_by(EVENT) %>%
        mutate(rank_combined = mean(c(rank,rank_diff_spliced))) %>% 
        arrange(rank_combined) %>% 
        ungroup()
    
    plts = list()
    plts[['diffsplicing_vs_dependency-recurrent_violin']] = df %>%
        ggviolin(x='event_freq', y='event_zscore') + 
        geom_point(data=df %>% filter(abs(event_zscore)>3)) + 
        geom_text_repel(aes(label=event_gene),
                        df %>% filter(abs(event_zscore)>4.5)) +
        labs(x='Diff. Splicing Freq.', y='Z-score')
    
    plts[['diffsplicing_vs_dependency-recurrent_rank']] = df %>%
        ggviolin(x='event_freq', y='rank_combined', trim=TRUE, add='median') +
        scale_y_reverse() +
        labs(x='Diff. Splicing Freq.', y='Combined Rank')
    
    # scatter
    X = event_properties_tcga %>% 
        filter(is_diff_spliced) %>% 
        inner_join(event_dependency, by=c('EVENT','GENE'))
    for (cancer in unique(X[['cancer_type']])){
        tmp = X %>% filter(cancer_type==cancer &
                           abs(`das-median_diff`) > 10 & 
                           abs(event_zscore) > 1.96
                          ) %>% 
            mutate(
                rank_fitness = rank(-abs(event_zscore)),
                rank_das = rank(-abs(`das-median_diff`)),
                rank_combined = rank(-(rank_fitness + rank_das))
            ) %>% arrange(rank_combined)
        
        plts[[sprintf('diffsplicing_vs_dependency-scatter-%s',cancer)]] = tmp %>%
            ggscatter(x='event_zscore', y='das-median_diff', size=1, alpha=0.5) +
            geom_text_repel(
                aes(label=event_gene),
                tmp %>% filter(`das-median_diff` > 10 & event_zscore < 1.96 |
                               `das-median_diff` < 10 & event_zscore > 1.96)) +
            ggtitle(cancer) +
            labs(x='Z-score', y=TeX('$\\Delta PSI$'))
        
        mat = tmp %>% column_to_rownames('event_gene') %>% dplyr::select(all_of(PROPERTIES_OI))
        plts[[sprintf('diffsplicing_vs_dependency-heatmap-%s',cancer)]] = pheatmap(t(mat), scale = 'row', cutree_cols = 2, silent=TRUE, main=cancer)
    }
       
    
    # top events with same direction across cancers
    event_freq = event_properties_tcga %>%
                    filter(is_diff_spliced) %>% 
                    inner_join(event_dependency, by=c('EVENT','GENE')) %>%
                    filter(grepl('HsaEX',event_gene)) %>%
                    mutate(direction = ifelse((`das-median_diff`<0), 'Spliced Out', 'Spliced In')) %>%
                    distinct(cancer_type, event_gene, GENE, direction, event_zscore) %>% 
                    count(event_gene, GENE, event_zscore, direction) %>%
                    mutate(freq=n/n_cancer_types) %>%
                    arrange(desc(freq)) %>%
                    filter(freq>0.5)
    events_oi = event_freq %>% group_by(direction) %>% slice_max(order_by=freq, n=10) %>% arrange(direction) %>% ungroup() %>% pull(event_gene)
    
    X = event_properties_tcga %>% 
        inner_join(event_dependency, by=c('EVENT','GENE')) %>% 
        filter(event_gene %in% events_oi)
    labels = X %>% distinct(event_gene,event_zscore) %>% 
              mutate(event_zscore=round(event_zscore,1),
                     color=ifelse(abs(event_zscore)>1.96,'red','black'))
    
    plts[['diffsplicing_vs_dependency-top_das']] = X %>% 
        ggstripchart(x='event_gene',y='das-median_diff',color='cancer_type',palette=palette,order=events_oi) +
        theme_pubr(x.text.angle = 70) + 
        geom_text(aes(x=event_gene, y=80, label=event_zscore),
                  labels,
                  size=3.5, color=labels$color) +
        geom_hline(yintercept = c(-5,5), linetype='dashed') +
        labs(x='Event & Gene', y=TeX('$\\Delta PSI$'),color='Cancer Type')
    
    # extreme event dependency
    X = event_properties_tcga %>% 
        inner_join(event_dependency, by=c('EVENT','GENE')) %>% 
        filter(abs(event_zscore)>1.96 & event_gene %in% event_freq[['event_gene']]) %>% filter(!is.na(`das-median_diff`))
    labels = X %>% distinct(event_gene,event_zscore) %>% 
              mutate(event_zscore=round(event_zscore,1),
                     color=ifelse(abs(event_zscore)>1.96,'red','black'))
    
    plts[['diffsplicing_vs_dependency-top_das_extreme']] = X %>% 
        ggstripchart(x='event_gene',y='das-median_diff',
                     color='cancer_type',palette=palette) +
        theme_pubr(x.text.angle = 70) + 
        geom_text(aes(x=event_gene, y=80, label=event_zscore),
                  labels,
                  size=3.5, color=labels$color) +
        geom_hline(yintercept = c(-5,5), linetype='dashed') +
        labs(x='Event & Gene', y=TeX('$\\Delta PSI$'),color='Cancer Type')
    
    return(plts)
}


make_plots = function(event_dependency, cancer_related_as, crispr_screen, event_properties_tcga){
    plts = list(
        plot_event_dependency(event_dependency, cancer_related_as),
        plot_crispr_validation(event_dependency, crispr_screen),
        plot_diffsplicing_dependency(event_dependency, event_properties_tcga)
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
    # save heatmaps
    heatmaps = grep('heatmap',names(plts),value = TRUE)
    lapply(heatmaps, function(plt_name){
        save_plot(plts[[plt_name]], plt_name, '.pdf', figs_dir, 300, width=15, height=10)
    })
        
    # save scatters
    scatters = grep('scatter',names(plts),value = TRUE)
    lapply(scatters, function(plt_name){
        save_plot(plts[[plt_name]], plt_name, '.png', figs_dir, 300, width=15, height=15)
    })
    
    # save manually
    save_plot(plts[['event_dependency-volcano']], 'event_dependency-volcano', '.png', figs_dir, 300, width=20, height=10)
    save_plot(plts[['event_dependency-event_pvalue']], 'event_dependency-event_pvalue', '.pdf', figs_dir, 300, width=12, height=6)
    save_plot(plts[['event_dependency-cancer_related_genes']], 'event_dependency-cancer_related_genes', '.pdf', figs_dir, 300, width=6, height=6)
    save_plot(plts[['event_dependency-validation_correlations']], 'event_dependency-validation_correlations', '.pdf', figs_dir, 300, width=12, height=6)
    save_plot(plts[['event_dependency-validation_scatter']], 'event_dependency-validation_scatter', '.png', figs_dir, 300, width=18, height=10)
    save_plot(plts[['diffsplicing_vs_dependency-recurrent_violin']], 'diffsplicing_vs_dependency-recurrent_violin', '.pdf', figs_dir, 300, width=6, height=6)
    save_plot(plts[['diffsplicing_vs_dependency-recurrent_rank']], 'diffsplicing_vs_dependency-recurrent_rank', '.pdf', figs_dir, 300, width=6, height=6)
    save_plot(plts[['diffsplicing_vs_dependency-top_das']], 'diffsplicing_vs_dependency-top_das', '.pdf', figs_dir, 300, width=10, height=6)
    save_plot(plts[['diffsplicing_vs_dependency-top_das_extreme']], 'diffsplicing_vs_dependency-top_das_extreme', '.pdf', figs_dir, 300, width=10, height=6)
    
}


main = function(){
    args = getParsedArgs()
    gonatopoulos_exon_crispr_screen_file = args$gonatopoulos_exon_crispr_screen_file
    thomas_crispr_screen = args$thomas_crispr_screen
    thomas_event_mapping_file = args$thomas_event_mapping_file
    event_dependency_file = args$event_dependency_file
    cancer_related_as = args$cancer_related_as_file
    event_properties_tcga = args$event_properties_tcga_file
    figs_dir = args$figs_dir
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    event_dependency = read_tsv(event_dependency_file) %>%
        filter(n_obs>MIN_OBS) %>%
        mutate(
            event_zscore=event_coefficient/event_stderr,
            event_gene = paste0(EVENT,'_',GENE),
            log10_pvalue = (-1)*log10(event_pvalue),
            event_type = ifelse(grepl('EX',EVENT),'EX','ALT')
        ) %>%
        group_by(event_type) %>%
        mutate(event_fdr = p.adjust(event_pvalue, method='fdr'),
               log10_fdr = -log10(event_fdr))
        
    cancer_related_as = read_tsv(cancer_related_as_file)
    crispr_screen = load_crispr_screens(
        gonatopoulos_exon_crispr_screen_file, 
        thomas_exon_crispr_screen_file, 
        thomas_event_mapping_file,
        psi_ccle_file,
        genexpr_ccle_file,
        event_annotation_file
    )
    #event_properties_tcga = read_tsv(event_properties_tcga_file) %>% 
                    filter(!is.na(EVENT)) %>% 
                    distinct(EVENT, GENE, cancer_type, .keep_all=TRUE) %>%
                    filter((!is.na(`das-median_diff`) | !is.na(`surv-coefficient`))) %>% 
                    distinct() 
    
    # plot
    plts = make_plots(event_dependency, cancer_related_as, crispr_screen, event_properties_tcga)

    # save
    save_plots(plts, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}