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
# - pvalue distributions of linear models? p-values vs n. observations?
# - scoring of linear models? spearman vs pearson? 
# - How many events don't have enough observations?
# - to which biological processes are significant events related?
# - which events/genes/interaction/intercepts have extreme z-scores?
# - Interesting examples? Cancer fitness, known drug targets, splicing factors
# - project cell lines based on real fitness, fitness predictions vs PSI (impute) and gene expression (make a script)
#  

require(tidyverse)
require(writexl)
require(ggpubr)
require(cowplot)
require(scattermore)
require(ggrepel)
require(clusterProfiler)
require(ComplexHeatmap)
require(ggplotify)
require(grid)
require(umap)
require(extrafont)
require(gtools)
require(tidytext)

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

# variables
THRESH_ZSCORE = 1.96
THRESH_PVALUE = 0.05
THRESH_LR_PVALUE = 0.005
THRESH_CORR = 0.2
SIZE_CTL = 100

# Development
# -----------
# PREP_DIR = file.path(ROOT,'data','prep')
# RESULTS_DIR = file.path(ROOT,'results','model_splicing_dependency')
# models_file = file.path(RESULTS_DIR,'files','models_gene_dependency-EX','model_summaries.tsv.gz')
# ccle_stats_file = file.path(PREP_DIR,'stats','CCLE.tsv.gz')
# rnai_file = file.path(PREP_DIR,'demeter2','CCLE.tsv.gz')
# spldep_file = file.path(RESULTS_DIR,'files','splicing_dependency-EX','mean.tsv.gz')
# msigdb_dir = file.path(ROOT,'data','raw','MSigDB','msigdb_v7.4','msigdb_v7.4_files_to_download_locally','msigdb_v7.4_GMTs')
# protein_impact_file = file.path(ROOT,'data','raw','VastDB','PROT_IMPACT-hg38-v3.tab.gz')
# gene_mut_freq_file = file.path(ROOT,'data','prep','gene_mutation_freq','CCLE.tsv.gz')
# event_mut_freq_file = file.path(ROOT,'data','prep','event_mutation_freq','CCLE-EX.tsv.gz')
# figs_dir = file.path(RESULTS_DIR,'figures','model_selection')
# possible_interactions_file = file.path(ROOT,'support','possible_pairwise_interaction_categories.tsv')
# cancer_events_file = file.path(ROOT,'support','cancer_events.tsv')
# indices_file = file.path(RESULTS_DIR,'files','correlation_spldep_indices-EX.tsv.gz')

##### FUNCTIONS #####
get_sets = function(df, set_names, set_values){
    sets = df[,c(set_values,set_names)] %>%
        distinct() %>%
        with(., split(get(set_values),get(set_names)))
    sets = sapply(sets, unique, simplufy=FALSE)
    return(sets)
}


run_enrichment = function(genes, events, universe, ontologies){
    enrichments = list()
    enrichments[['hallmarks']] = enricher(genes, TERM2GENE=ontologies[['hallmarks']], universe=universe[['genes']])
    enrichments[['oncogenic_signatures']] = enricher(genes, TERM2GENE=ontologies[['oncogenic_signatures']], universe=universe[['genes']])
    enrichments[['GO_BP']] = enricher(genes, TERM2GENE=ontologies[['GO_BP']], universe=universe[['genes']])
    enrichments[['protein_impact']] = enricher(events, TERM2GENE=ontologies[['protein_impact']], universe=universe[['events']])
    
    return(enrichments)
}


run_enrichments = function(genes_lists, events_lists, universe, ontologies){
    enrichments = sapply(names(genes_lists), function(cond){
            genes = genes_lists[[cond]]
            events = events_lists[[cond]]
            run_enrichment(genes, events, universe, ontologies)
        }, simplify=FALSE)
    return(enrichments)
}


get_enrichment_result = function(enrich_list, thresh=0.05){
    ## groups are extracted from names
    ontos = names(enrich_list[[1]])
    groups = names(enrich_list)
    results = sapply(
        ontos, function(onto){
        result = lapply(groups, function(group){
            result = enrich_list[[group]][[onto]]
            if(!is.null(result)){
                result = result@result
                result$Cluster = group
            }
            return(result)
        })
        result[sapply(result, is.null)] = NULL
        result = do.call(rbind,result)
        ## filter by p.adjusted
        result = result %>% filter(p.adjust<thresh)
    }, simplify=FALSE)
    
    return(results)
}


plot_model_selection = function(models, rnai, spldep, gene_mut_freq, event_mut_freq, ontologies, cancer_events){
    # prep data
    rnai = rnai %>% filter(index%in%models[['GENE']]) %>% column_to_rownames('index')
    spldep = spldep %>% filter(index%in%models[['EVENT']]) %>% column_to_rownames('index')
    rnai_stats = data.frame(
            rnai_med = apply(rnai,1,median, na.rm=TRUE),
            rnai_std = apply(rnai,1,sd, na.rm=TRUE),
            n_missing = rowSums(is.na(rnai))
        ) %>% rownames_to_column('GENE')
        
    plts = list()
    # as a positive control, we got ~300 exons whose modulation affects cell
    # proliferation. Of course, we need them to belong to genes that vary to
    # be likely to have an effect across cancers. So we took the top 100.
    ctl_pos_genes = cancer_events %>% pull(GENE) %>% unique()
    plts[['model_selection-deps_sorted_vs_std_ctl_pos']] = rnai_stats %>% 
        arrange(-rnai_std) %>% 
        mutate(index=row_number(), 
               is_ctl=GENE %in% ctl_pos_genes) %>%
        ggscatter(x="index", y="rnai_std", color="is_ctl", 
                  size="is_ctl", palette=c('grey','darkgreen')) +
        scale_size_discrete(range=c(0.5,1)) +
        labs(x='Index', y='Gene Std. Demeter2', 
             color='In Positive Control', size='In Positive Control', 
             title=sprintf("Total Control Genes = %s", length(ctl_pos_genes))) +
        theme(aspect.ratio=1)
    
    plts[['model_selection-deps_med_vs_std_ctl_pos']] = rnai_stats %>% 
        mutate(is_ctl=GENE %in% ctl_pos_genes) %>%
        ggplot(aes(x=rnai_med, y=rnai_std, color=is_ctl)) +
        geom_scattermore(pixels = c(1000,1000), pointsize=8) +
        theme_pubr() +
        color_palette(palette = c('grey','darkgreen')) +
        labs(x='Gene Median Demeter2', y='Gene Std. Demeter2', 
             title=sprintf("Total Control Genes = %s", length(ctl_pos_genes)),
             color='In Positive Control') +
        theme(aspect.ratio=1)
    
    ctl_pos = rnai_stats %>% 
        filter(GENE %in% ctl_pos_genes) %>% 
        slice_max(order_by=rnai_std, n=SIZE_CTL) %>% 
        left_join(cancer_events, by="GENE") %>%
        pull(EVENT) %>%
        unique()
    
    plts[['model_selection-lr_pvalue_ctl_pos']] = models %>%
        mutate(is_ctl_pos = EVENT %in% ctl_pos) %>%
        ggviolin(x="is_ctl_pos", y="lr_pvalue", trim = TRUE, 
                 fill='is_ctl_pos', color=NA, palette = c('grey','darkgreen')) + 
        geom_boxplot(width=0.1) +
        stat_compare_means(method="wilcox.test") +
        guides(fill="none") + 
        labs(x='Is Positive Control', y='LR Test p-value')
    
    # as a negative control, we selected 100 uniform gene dependencies
    ctl_neg = rnai_stats %>% 
        filter(!(GENE %in% ctl_pos_genes)) %>% # do not include positive controls
        slice_min(order_by = rnai_std, n=SIZE_CTL) %>% # selected top 100
        pull(GENE)
    
    plts[['model_selection-deps_sorted_vs_std_ctl_neg']] = rnai_stats %>% 
        arrange(-rnai_std) %>% 
        mutate(index=row_number(), 
               is_ctl=GENE %in% ctl_neg) %>%
        ggscatter(x="index", y="rnai_std", color="is_ctl", 
                  size="is_ctl", palette=c('grey','brown')) +
        scale_size_discrete(range=c(0.5,1)) +
        labs(x='Index', y='Gene Std. Demeter2', 
             title=sprintf("Total Control Genes = %s", length(ctl_neg)),
             color="In Negative Control", size="In Negative Control") +
        theme(aspect.ratio=1)
    
    plts[['model_selection-deps_med_vs_std_ctl_neg']] = rnai_stats %>% 
        mutate(is_ctl=GENE %in% ctl_neg) %>%
        ggplot(aes(x=rnai_med, y=rnai_std, color=is_ctl)) +
        geom_scattermore(pixels = c(1000,1000), pointsize=8) +
        theme_pubr() +
        color_palette(palette = c('grey','brown')) +
        labs(x='Gene Median Demeter2', y='Gene Std. Demeter2', 
             title=sprintf("Total Control Genes = %s", length(ctl_neg)),
             color='In Negative Control') +
        theme(aspect.ratio=1)
        
    # we selected the ensemble of models by their ability to rank dependencies
    # using likelihood ratio tests' p-values as thresholds
    plts[['model_selection-lr_pvalue']] = models %>%
        gghistogram(x='lr_pvalue', bins=100, fill='grey', color=NA) +
        geom_vline(xintercept=median(models[['lr_pvalue']], na.rm=TRUE),
                   linetype='dashed') +
        labs(x='LR Test p-value', y='Count')
    
    plts[['model_selection-pearson_corr']] = models %>%
        gghistogram(x='pearson_correlation_mean', bins=100, fill='darkred', color=NA) +
        geom_vline(xintercept=median(models[['pearson_correlation_mean']], na.rm=TRUE),
                   linetype='dashed') +
        labs(x='Pearson Correlation Mean (Test Set)', y='Count')
    
    # find the threshold for the p-value
    threshs = c(
        1e-4,2.5e-4,5e-4,
        1e-3,2.5e-3,5e-3,
        1e-2,2.5e-2,5e-2,
        1e-1,2.5e-1,5e-1,
        1
    )
    eval_pvalue = lapply(threshs, function(thresh){
        models_possible = models %>%
            filter(pearson_correlation_mean>0)
        ctl_neg_corrected = ctl_neg[ctl_neg %in% (models_possible %>% pull(GENE))]
        ctl_pos_corrected = ctl_pos[ctl_pos %in% (models_possible %>% pull(EVENT))]
        
        models_filtered = models_possible %>% 
            filter(lr_pvalue<thresh) 
        
        # from each gene, pick de event-level model that generalizes the worst
        models_selected = models_filtered %>%
            group_by(GENE) %>%
            slice_min(order_by=pearson_correlation_mean, n=1)
        
        events = models_selected %>% pull(EVENT)
        genes = models_selected %>% pull(GENE)
        
        df = data.frame(
                thresh = thresh,
                total_events = nrow(models_filtered),
                total_genes = length(genes),
                total_ctl_neg = sum(genes %in% ctl_neg_corrected),
                total_ctl_pos = sum(models_filtered %>% pull(EVENT) %in% ctl_pos_corrected)
            )
        df[['tpr']] = df[['total_ctl_pos']] / length(ctl_pos_corrected)
        df[['fpr']] = df[['total_ctl_neg']] / length(ctl_neg_corrected)
        df[['tpr_vs_fpr']] = df[['tpr']] - df[['fpr']]
        return(df)
    })
    eval_pvalue = do.call(rbind,eval_pvalue)
    
    plts[['model_selection-pvalue_vs_n_genes']] = eval_pvalue %>%
        ggbarplot(x='thresh', y='total_genes', label=TRUE, lab.size=1, 
                  fill='darkblue', color=NA) +
        labs(x='Thresholds LR Test p-value', y='No. Genes Selected') +
        theme_pubr(x.text.angle=70)
    
    plts[['model_selection-pvalue_vs_n_events']] = eval_pvalue %>%
        ggbarplot(x='thresh', y='total_events', label=TRUE, lab.size=1, 
                  fill='lightblue', color=NA) +
        labs(x='Thresholds LR Test p-value', y='No. Events Selected') +
        theme_pubr(x.text.angle=70)
    
    plts[['model_selection-pvalue_vs_ctl_neg']] = eval_pvalue %>%
        ggbarplot(x='thresh', y='total_ctl_neg', label=TRUE, lab.size=1, 
                  fill='brown', color=NA) +
        labs(x='Thresholds LR Test p-value', y='No. Neg. Ctl. Genes Selected') +
        theme_pubr(x.text.angle=70)
    
    plts[['model_selection-pvalue_vs_ctl_pos']] = eval_pvalue %>%
        ggbarplot(x='thresh', y='total_ctl_pos', label=TRUE, lab.size=1, 
                  fill='darkgreen', color=NA) +
        labs(x='Thresholds LR Test p-value', y='No. Positive Control Exons Selected') +
        theme_pubr(x.text.angle=70)
    
    # tpr vs fpr vs median correlation
    plts[['model_selection-roc_curve']] = eval_pvalue %>% 
        ggline(x='fpr', y='tpr', linetype='dashed', numeric.x.axis=TRUE, size=0.25,
               label="thresh", repel=FALSE, font.label = c(6, "plain")) + 
        geom_abline(intercept=0, slope=1, linetype='dashed') + 
        labs(x='FPR', y='TPR') + 
        xlim(0,1) + 
        ylim(0,1)
    
    plts[['model_selection-tpr_vs_fpr']] = eval_pvalue %>% 
        ggbarplot(x='thresh', y='tpr_vs_fpr', fill="lightblue", color=NA) +
        labs(x='Thresholds LR Test p-value', y='TPR - FPR') +
        geom_text(aes(y=0.10, label=total_events), size=0.85, color="#994F00") +
        geom_text(aes(y=0.093, label=total_genes), size=0.85, color="#006CD1")
    
    # find the threshold for the mean pearson correlation
    threshs = seq(0,0.4,0.05)
    samples_oi = intersect(colnames(rnai),colnames(spldep))
    eval_corr = lapply(threshs, function(thresh){
        models_possible = models %>%
            filter(lr_pvalue < THRESH_LR_PVALUE)
        ctl_neg_corrected = ctl_neg[ctl_neg %in% (models_possible %>% pull(GENE))]
        ctl_pos_corrected = ctl_pos[ctl_pos %in% (models_possible %>% pull(EVENT))]
        
        models_filtered = models_possible %>% 
            filter(pearson_correlation_mean > thresh) 
        
        # from each gene, pick de event-level model that generalizes the worst
        models_selected = models_filtered %>%
            group_by(GENE) %>%
            slice_min(order_by=pearson_correlation_mean, n=1)
        
        events = models_selected %>% pull(EVENT)
        genes = models_selected %>% pull(GENE)
        
        idx = sample(1:length(samples_oi), size=250)
        tmp = lapply(samples_oi[idx], function(sample_oi){
            true_dep = rnai[genes,sample_oi]
            pred_dep = spldep[events,sample_oi]
            test = cor.test(true_dep, pred_dep, 
                       method='spearman', use='pairwise.complete.obs')
            df = data.frame(
                id = sample_oi,
                corr = test[['estimate']][['rho']],
                pvalue = test[['p.value']],
                thresh = thresh,
                total_events = nrow(models_filtered),
                total_genes = length(genes),
                total_ctl_neg = sum(genes %in% ctl_neg_corrected),
                total_ctl_pos = sum(models_filtered %>% pull(EVENT) %in% ctl_pos_corrected)
            )
            df[['tpr']] = df[['total_ctl_pos']] / length(ctl_pos_corrected)
            df[['fpr']] = df[['total_ctl_neg']] / length(ctl_neg_corrected)
            return(df)
        })
        tmp = do.call(rbind,tmp)
        return(tmp)
    })
    eval_corr = do.call(rbind, eval_corr)
    eval_corr = eval_corr %>% mutate(thresh_fct = as.factor(thresh))
    
    plts[['model_selection-pearson_corr_vs_spearman']] = eval_corr %>%
        ggviolin(x='thresh_fct', y='corr', fill='orange', color=NA, trim=TRUE) + 
        geom_boxplot(fill=NA, outlier.size=0.1, width=0.2) +
        labs(x='Thresholds Single-Model Avg. Pearson Correlation', y='Spearman Correlation') +
        theme_pubr(x.text.angle=70) +
        geom_text(
            aes(x=thresh_fct, y=1.04, label=total_events), 
            eval_corr %>% 
                distinct(thresh_fct,total_events), 
            size=1, color="#994F00") +
        geom_text(
            aes(x=thresh_fct, y=1.01, label=total_genes), 
            eval_corr %>% 
                distinct(thresh_fct,total_genes), 
            size=1, color="#006CD1")
    
    # with this set of models, we expect to see certain properties
    models = models %>%
        drop_na(is_selected)
    
    ## selected models come from variant events
    plts[['model_selection-selected_vs_event_std']] = models %>%
        ggviolin(x='is_selected', y='event_std', color=NA, trim=TRUE,
                 fill='is_selected', palette='npg') +
        geom_boxplot(fill=NA, outlier.size = 0.5) +
        stat_compare_means(method='wilcox.test') +
        guides(color='none', fill='none') +
        labs(x='Selected Model', y='Event Std.')
    
    ## selected models are in genes less prone to have deleterious mutations
    X = models %>%
        left_join(gene_mut_freq, by='GENE') %>%
        group_by(Variant_Classification,GENE,mut_freq_per_kb) %>%
        summarize(is_selected = any(is_selected)) %>%
        drop_na() %>%
        ungroup()
        
    plts[['model_selection-mutation_gene_count']] = X %>% 
        count(Variant_Classification, is_selected) %>%
        ggbarplot(x='Variant_Classification', y='n', 
                  label=TRUE, lab.size=2, palette='npg',
                  fill='is_selected', color=NA, position=position_dodge(0.9)) + 
        yscale('log10', .format=TRUE) + 
        labs(x='Mutation Effect', y='No. Genes', 
             fill='Selected Model') +
        theme_pubr(x.text.angle=70)
    
    plts[['model_selection-mutation_gene_frequency']] = X %>% 
        ggplot(aes(x=Variant_Classification, y=mut_freq_per_kb, 
                   group=interaction(Variant_Classification,is_selected))) +
        geom_boxplot(aes(fill=is_selected), outlier.size=0.1, 
                     position=position_dodge(0.7)) +
        stat_compare_means(aes(group=is_selected), method='wilcox.test', 
                           label='p.signif', size=2) +
        yscale('log10', .format=TRUE) + 
        fill_palette('npg') +
        labs(x='Mutation Effect', y='log10(Mut. Freq. per Kb)', 
             fill='Selected Model') +
        theme_pubr(x.text.angle=70)
    
    
    # how often do selected exons get hit when the gene is mutated?
    X = models %>%
        left_join(event_mut_freq, by=c('EVENT','GENE')) %>%
        # keep only genes with a selected exon
        group_by(GENE) %>%
        filter(any(is_selected)) %>%
        ungroup() %>%
        left_join(ontologies[['protein_impact']], by='EVENT') %>%
        drop_na(Variant_Classification, term_clean, is_selected)

    notsel_mut_freq = X %>%
        # average mutation frequency per kb of not selected events
        filter(!is_selected) %>%
        group_by(Variant_Classification, GENE) %>%
        summarize(notsel_mut_freq = mean(event_mut_freq_per_kb, na.rm=TRUE))
    
    X = X %>% 
        left_join(notsel_mut_freq, by=c("Variant_Classification","GENE")) %>%
        # Fold change difference 
        mutate(fc_mut_freq = log2(event_mut_freq_per_kb / notsel_mut_freq))
    
    ## do the same with a null distribution of 1000 random exons
    null = models %>%
        left_join(ontologies[['protein_impact']], by='EVENT') %>%
        left_join(event_mut_freq, by=c('EVENT','GENE')) %>%
        # shuffle selected exons by protein impact and mutation variant
        group_by(Variant_Classification,GENE) %>%
        mutate(is_selected = sample(is_selected)) %>%
        ungroup() %>%
        # keep only genes with a selected exon
        group_by(GENE) %>%
        filter(any(is_selected)) %>%
        ungroup() %>%
        drop_na(Variant_Classification, term_clean, is_selected)
    
    notsel_mut_freq = null %>%
        # average mutation frequency per kb of not selected events
        filter(!is_selected) %>%
        group_by(Variant_Classification, GENE) %>%
        summarize(notsel_mut_freq = mean(event_mut_freq_per_kb, na.rm=TRUE))
    
    null = null %>% 
        left_join(notsel_mut_freq, by=c("Variant_Classification","GENE")) %>%
        # Fold change difference 
        mutate(fc_mut_freq = log2(event_mut_freq_per_kb / notsel_mut_freq))
    
    X = X %>% 
        mutate(dataset = "Real") %>% 
        bind_rows(null %>% mutate(dataset = "Random"))
    
    
    plts[['model_selection-mutation_event_count']] = X %>% 
        count(dataset, Variant_Classification, is_selected) %>%
        ggbarplot(x='Variant_Classification', y='n', 
                  label=TRUE, palette='npg', lab.size=2,
                  fill='is_selected', color=NA, position=position_dodge(0.9)) + 
        yscale('log10', .format=TRUE) + 
        labs(x='Mutation Effect', y='No. Events', fill='Selected Model') +
        theme_pubr(x.text.angle=70) +
        facet_wrap(~dataset, ncol=1)
    
    
    plts[['model_selection-mutation_event_frequency']] = X %>% 
        filter(is_selected) %>%
        ggplot(aes(x=Variant_Classification, y=fc_mut_freq, 
                   group=interaction(Variant_Classification,dataset))) +
        geom_boxplot(aes(fill=dataset), outlier.size=0.1, 
                     position=position_dodge(0.7)) +
        fill_palette("npg") + 
        geom_hline(yintercept=0, linetype='dashed') +
        stat_compare_means(aes(group=dataset), method='wilcox.test', 
                           label='p.signif', size=2) +
        labs(x='Mutation Effect', y='log2(FC Mut. Freq. per Kb)', 
             fill='Selected Model') +
        theme_pubr(x.text.angle=70)
    
    
    plts[['model_selection-mutation_event_frequency-by_protein_impact']] = X %>% 
        drop_na(fc_mut_freq) %>%
        filter(is_selected) %>%
        ggplot(aes(x=Variant_Classification, y=fc_mut_freq, 
                   group=interaction(Variant_Classification,dataset))) +
        geom_boxplot(aes(fill=dataset), outlier.size=0.1, 
                     position=position_dodge(0.7)) +
        fill_palette("npg") + 
        geom_hline(yintercept=0, linetype='dashed') +
        stat_compare_means(aes(group=dataset), method='wilcox.test', 
                           label='p.signif', size=2) +
        labs(x='Mutation Effect', y='log2(FC Mut. Freq. per Kb)', 
             fill='Selected Model') +
        theme_pubr(x.text.angle=70) +
        facet_wrap(~term_clean)
    
    # what are the exons/genes selected?
    ## protein impact
    prot_imp = models %>% 
        left_join(ontologies[['protein_impact']], by='EVENT')
    
    plts[['model_selection-protein_impact-counts']] = prot_imp %>%
        filter(is_selected) %>%
        group_by(is_selected,term) %>% 
        summarize(n=n()) %>% 
        arrange(n) %>%
        ggbarplot(x='term', y='n', label=TRUE, lab.size=2,
                  fill='is_selected', color=NA, palette='npg') + 
        theme_pubr(x.text.angle = 45, legend='right') +
        guides(color='none') + 
        labs(x='Protein Impact', y='Count')
    
    plts[['model_selection-protein_impact-freqs']] = prot_imp %>%
        group_by(term,is_selected) %>%
        summarize(n=n()) %>%
        drop_na() %>%
        mutate(freq=n/sum(n)) %>%
        filter(is_selected) %>%
        ggbarplot(x='term', y='freq', fill='is_selected', 
                  color=NA, palette='npg') + 
        theme_pubr(x.text.angle = 45, legend='right') +
        guides(color='none') + 
        labs(x='Protein Impact', y='Proportion')
    
#      prot_imp %>%
#         group_by(is_selected,term) %>%
#         summarize(n=n()) %>%
#         drop_na() %>%
#         mutate(freq=n/sum(n)) %>%
#         ggbarplot(x='is_selected', y='freq', fill='term', color=NA,
#                   palette=get_palette('BrBG', k=length(unique(prot_imp[['term']])))) + 
#         theme_pubr(x.text.angle = 45, legend='right') +
#         guides(color='none') + 
#         labs(x='Selected Model', y='Proportion', fill='Protein Impact')
    
    
    plts[['model_selection-protein_impact_clean-counts']] = prot_imp %>%
        filter(is_selected) %>%
        group_by(is_selected,term_clean) %>%
        summarize(n=n()) %>%
        arrange(n) %>%
        drop_na() %>%
        ggbarplot(x='term_clean', y='n', label=TRUE, lab.size=2,
                  fill='is_selected', color=NA, palette='npg') + 
        theme_pubr(x.text.angle = 45, legend='right') +
        guides(color='none') + 
        labs(x='Protein Impact', y='Count')
    
    plts[['model_selection-protein_impact_clean-freqs']] = prot_imp %>%
        group_by(term_clean,is_selected) %>%
        summarize(n=n()) %>%
        mutate(freq=n/sum(n)) %>%
        filter(is_selected) %>%
        ggbarplot(x='term_clean', y='freq', fill='is_selected', color=NA,
                  palette='npg') + 
        theme_pubr(x.text.angle = 45, legend='right') +
        guides(color='none') + 
        labs(x='Selected Model', y='Proportion')
    
    ## GSEA
    genes_oi = models %>% 
            filter(is_selected) %>% 
            pull(GENE) %>% 
            unique()
    events_oi = models %>% 
            filter(is_selected) %>% 
            pull(EVENT)
    universe = list(
        'genes' = models %>% pull(GENE) %>% unique(),
        'events' = models %>% pull(EVENT) %>% unique()
    )
    enrichment = run_enrichment(genes_oi, events_oi, universe, ontologies)
    enrichment[sapply(enrichment, nrow)<1] = NULL
    plts_enrichment = lapply(names(enrichment), function(e_name){
        res = enrichment[[e_name]]
        plts = list()
        plts[['dotplot']] = dotplot(res)
        plts[['cnetplot']] = cnetplot(res, cex_label_category=0.5, cex_label_gene=0.5)
        names(plts) = sprintf('%s-%s',e_name,names(plts))
        return(plts)
    })
    plts_enrichment = do.call(c,plts_enrichment)
    names(plts_enrichment) = sprintf('model_selection-enrichment-%s',
                                     names(plts_enrichment))
    plts = c(plts, plts_enrichment)
    plts[['model_selection-enrichment-GO_BP-dotplot']] = plts[['model_selection-enrichment-GO_BP-dotplot']] + 
        theme_pubr()
    
    # are selected exons associated to transcriptomic indices?
    X = models %>%
        left_join(indices, by=c("EVENT"="index")) %>% 
        mutate(abs_corr = abs(correlation),
               corr_sign = ifelse(sign(correlation)>0,"Positive","Negative")) %>%
        drop_na(correlation, is_selected)
    
    plts[["model_selection-indices-violin"]] = X %>% 
        ggplot(aes(x=index_name, y=correlation, 
                   group=interaction(index_name,is_selected))) + 
        geom_violin(aes(fill=is_selected), color=NA) + 
        geom_boxplot(fill=NA, outlier.size=0.1, 
                     width=0.2, position=position_dodge(0.9)) +
        stat_compare_means(aes(group=is_selected), method='wilcox.test', 
                           label='p.signif', size=2) + 
        fill_palette("npg") + 
        geom_hline(yintercept=0, linetype="dashed") + 
        #facet_wrap(~corr_sign) + 
        theme_pubr() + 
        labs(x="Transcriptomic Index", y="Spearman Correlation")
    
    plts[["model_selection-indices-top_pos"]] = X %>% 
        filter(is_selected) %>%
        group_by(index_name) %>%
        slice_max(correlation, n=15) %>%
        ungroup() %>%
        mutate(event_gene = as.factor(event_gene),
               name = reorder_within(event_gene, correlation, index_name)) %>%
        ggbarplot(x="name", y="correlation", fill="#4DBBD5FF", color=NA) +
        facet_wrap(~index_name, scales='free') +
        scale_x_reordered() +
        labs(x="Event & Gene", y="Spearman Correlation") +
        coord_flip()
    
    plts[["model_selection-indices-top_neg"]] = X %>% 
        filter(is_selected) %>%
        group_by(index_name) %>%
        slice_max(-correlation, n=15) %>%
        ungroup() %>%
        mutate(event_gene = as.factor(event_gene),
               name = reorder_within(event_gene, correlation, index_name)) %>%
        ggbarplot(x="name", y="correlation", fill="#4DBBD5FF", color=NA) +
        facet_wrap(~index_name, scales='free') +
        scale_x_reordered() +
        labs(x="Event & Gene", y="Spearman Correlation") +
        coord_flip()
    
    return(plts)
}


get_interaction_categories = function(models, possible_interactions){
    # study interaction between event inclusion and gene expression
    possible_interactions = possible_interactions %>% 
        mutate(combined=paste0(beta_a,beta_b,beta_ab))
    
    ## interactions
    zscores = models %>% 
        column_to_rownames("event_gene") %>% 
        dplyr::select(paste0(c('event','gene','interaction'),"_zscore"))
    intcats = (abs(zscores) > THRESH_ZSCORE) * sign(zscores) # p-value<0.05
    intcats[is.na(intcats)] = 0
    intcats = intcats %>% 
        rownames_to_column("event_gene") %>%
        mutate(combined = paste0(event_zscore,gene_zscore,interaction_zscore)) %>%
        left_join(possible_interactions,by='combined') %>%
        dplyr::select(event_gene,category,combined) %>%
        rename(interaction_category=category, interaction_subcategory=combined)
    
    return(intcats)
}


plot_interactions = function(models, ontologies){
    
    X = models %>% 
        mutate(interaction_subcategory_clean = gsub('-1','1',interaction_subcategory)) %>%
        left_join(ontologies[['protein_impact']], by='EVENT')
        
    plts = list()
    
    # across event types, what's the frequency of each type of psi-tpm interaction?
    plts[['interactions-counts']] = X %>% 
        filter(is_selected) %>%
        count(event_type, interaction_category) %>%
        ggbarplot(x='interaction_category', y='n', facet.by = 'event_type', 
                  label=TRUE, fill='interaction_category', color=NA, palette='Set1') + 
        guides(fill='none') + 
        yscale('log10', .format=TRUE) + 
        labs(x='Interaction Category', y='Counts') +
        theme_pubr(x.text.angle = 45, border=TRUE)
    
    plts[['interactions-freqs_subcategories']] = X %>%
        filter(is_selected) %>%
        group_by(interaction_category, interaction_subcategory_clean) %>%
        summarize(n=n()) %>%
        mutate(freq=n/sum(n)) %>%
        ggbarplot(x='interaction_category', y='freq', fill='interaction_subcategory_clean', 
                  color=NA, palette='Set2') +
        theme_pubr(x.text.angle = 45, legend = 'right') +
        labs(x='Interaction Category', y='Proportion', fill='Interaction\nSubcategory')
    
    # Is some interaction category associated to some protein impact?
    plts[['interactions-protein_impact-freqs']] = X %>% 
        filter(is_selected) %>%
        group_by(term,interaction_category) %>% 
        summarize(n=n()) %>% 
        drop_na() %>%
        mutate(freq= n/sum(n)) %>% 
        ggbarplot(x='term', y='freq', fill='interaction_category', 
                  color=NA, palette='Set1') +
        labs(x='Interaction Category', y='Proportion', fill='Interaction\nCategory') +
        theme_pubr(x.text.angle = 45, legend = 'right')
    
    plts[['interactions-protein_impact_clean-freqs']] = X %>% 
        filter(is_selected) %>%
        group_by(term_clean,interaction_category) %>% 
        summarize(n=n()) %>% 
        drop_na() %>%
        mutate(freq= n/sum(n)) %>% 
        ggbarplot(x='term_clean', y='freq', fill='interaction_category', 
                  color=NA, palette='Set1') +
        labs(x='Interaction Category', y='Proportion', fill='Interaction\nCategory') +
        theme_pubr(x.text.angle = 45, legend = 'right')
    
    # PCA of selected model coefficients
    pca = X %>%
        filter(is_selected) %>%
        column_to_rownames('event_gene') %>%
        dplyr::select(matches('coefficient')) %>%
        mutate_all(scale) %>%
        t() %>%
        prcomp()
    pcs = pca[['rotation']] %>% 
        as.data.frame() %>% 
        rownames_to_column('event_gene') %>%
        left_join(X, by='event_gene')
    plts[['interactions-category-pca']] = pcs %>% 
        ggscatter(x='PC1', y='PC2', color='interaction_category', palette='Set1') 
    
    umaps = pca[['rotation']] %>% umap()
    umaps = umaps[['layout']] %>% 
        as.data.frame() %>% 
        rownames_to_column('event_gene') %>%
        left_join(X, by='event_gene')
    plts[['interactions-category-umap']] = umaps %>% 
        ggscatter(x='V1', y='V2', color='interaction_category', palette='Set1') + 
        labs(x='UMAP1', y='UMAP2')
    
    # is any of the interaction classes enriched in some BP term? (TODO)
    genes_oi = get_sets(X %>% filter(is_selected), 'interaction_category', 'GENE')
    events_oi = get_sets(X %>% filter(is_selected), 'interaction_category', 'EVENT')
    universe = list(
        'genes' = X %>% pull(GENE) %>% unique(),
        'events' = X %>% pull(EVENT) %>% unique()
    )
    ## how repetitive are genes among interaction categories?
    m = genes_oi %>% list_to_matrix() %>% make_comb_mat()
    plts[['interactions-category-upset']] = UpSet(m, comb_order = order(comb_size(m)))
    plts[['interactions-category-upset']] = as.ggplot(grid.grabExpr(draw(plts[['interactions-category-upset']])))

    ## GSEA
    enrichments = run_enrichments(genes_oi, events_oi, universe, ontologies)
    results = get_enrichment_result(enrichments)
    results = results[sapply(results,nrow)>0]
    plts_enrichment = sapply(names(results), function(onto){
        result = results[[onto]]
        res = new("compareClusterResult", compareClusterResult = result)
        plt = dotplot(res) + labs(title=onto, x='Interaction Category')
        return(plt)
    }, simplify=FALSE)
    names(plts_enrichment) = sprintf('interactions-enrichment-%s',
                                     names(plts_enrichment))
    plts = c(plts,plts_enrichment)
    
    return(plts)
}


plot_examples = function(models, protein_impact){
    X = models %>%
        filter(is_selected) %>%
        arrange(pearson_correlation_mean) %>%
        mutate(index=row_number()) %>%
        left_join(protein_impact, by='EVENT')
    
    plts = list()
    
    cols_oi = c('EVENT','GENE','ENSEMBL','interaction_category',
                'interaction_subcategory','term','event_type','lr_pvalue')
    
    plts[['examples-top_overall-scatter']] = X %>%
        ggplot(aes(x=index, y=pearson_correlation_mean)) +
        geom_scattermore(pixels=c(1000,1000), pointsize=8) +
        geom_text_repel(
            aes(label=event_gene),
            X %>%
            slice_max(order_by=pearson_correlation_mean, n=15),
            max.overlaps = 50,
            size=1,
            segment.size=0.1
        ) +
        theme_pubr() +
        labs(x='Index', y='Avg. Pearson Correlation (Test Set)', 
             title=sprintf('No. Selected Models=%s',nrow(X)))
    
     plts[['examples-top_overall-table']] = X %>%
        dplyr::select(cols_oi) %>%
        slice_min(order_by=lr_pvalue, n=100) %>%
        ggtexttable()
    
    plts[['examples-top_by_protein_impact-scatter']] = X %>%
        ggscatter(x='index', y='lr_pvalue', facet.by = 'term_clean') +
        yscale('log10', .format=TRUE) +
        geom_text_repel(
            aes(label=event_gene),
            X %>%
            group_by(term_clean) %>%
            slice_min(order_by=lr_pvalue, n=5),
            max.overlaps = 50
        ) +
        labs(x='Index', y='log10(LR Test p-value)', 
             title=sprintf('No. Selected Models=%s',nrow(X)))
    
     plts[['examples-top_by_protein_impact-table']] = X %>%
        group_by(term_clean) %>%
        slice_min(order_by=lr_pvalue, n=5) %>%
        dplyr::select(cols_oi) %>%
        ggtexttable()
    
    plts[['examples-sign_events_gene-bar']] = X %>%
        ungroup() %>%
        count(GENE) %>%
        count(n) %>%
        ggbarplot(x='n', y='nn', numeric.x.axis=TRUE, 
                  fill='orange', color=NA) +
        labs(x='Events per Gene', y='Count') + 
        geom_text_repel(
            aes(x=n,y=y,label=GENE),
            X %>%
            ungroup() %>%
            count(GENE) %>%
            group_by(n) %>%
            slice_head(n=5) %>%
            mutate(y=100)
        )
    
    return(plts)
}

plot_tumorigenesis = function(models, spldep){
    # there are exons with different behaviors dependencing 
    # on the splicing dependency:
    # - oncoexons: SplDep<0
    # - tumor suppressor exons: SplDep>0
    # - double agent exons: SplDep ambiguous
    
    X = spldep %>% 
        filter(index%in%(models %>% filter(is_selected) %>% pull(EVENT))) %>% 
        column_to_rownames('index')
    spldep_stats = data.frame(
            med = apply(X,1,median, na.rm=TRUE),
            std = apply(X,1,sd, na.rm=TRUE),
            min = apply(X,1,min, na.rm=TRUE),
            max = apply(X,1,max, na.rm=TRUE),
            q95 = apply(X,1,quantile, na.rm=TRUE, probs=0.95),
            q05 = apply(X,1,quantile, na.rm=TRUE, probs=0.05),
            range = abs(apply(X,1,max, na.rm=TRUE)) - abs(apply(X,1,min, na.rm=TRUE)),
            n_missing = rowSums(is.na(X))
        ) %>% rownames_to_column('EVENT') %>%
        left_join(models[c('EVENT','event_gene')], by='EVENT')
    
    plts = list()

    plts[['tumorigenesis-scatter']] = spldep_stats %>% 
        ggplot(aes(x=q05, y=q95)) +
        geom_scattermore(pointsize=8, pixels=c(1000,1000)) +
        theme_pubr() +
        geom_hline(yintercept=0, linetype='dashed') + 
        geom_vline(xintercept=0, linetype='dashed') +
        labs(x='0.05 Quantile', y='0.95 Quantile') +
        geom_text_repel(aes(label=event_gene), size=1, segment.size=0.1,
                        spldep_stats %>% slice_max(order_by = q05*q95, n=5)) +
        geom_text_repel(aes(label=event_gene), size=1, segment.size=0.1,
                        spldep_stats %>% filter(q95>1)) +
        geom_text_repel(aes(label=event_gene), size=1, segment.size=0.1,
                        spldep_stats %>% filter(q95>0 & q05<(-0.9)))
    
    return(plts)
}

make_plots = function(models, rnai, spldep, gene_mut_freq, 
                      event_mut_freq, ontologies, cancer_events){
    plts = list(
        plot_model_selection(models, rnai, spldep, gene_mut_freq, 
                             event_mut_freq, ontologies, cancer_events),
        plot_interactions(models, ontologies),
        plot_examples(models, protein_impact=ontologies[['protein_impact']]),
        plot_tumorigenesis(models, spldep)
    )
    plts = do.call(c,plts)
    return(plts)
}


save_plt = function(plts, plt_name, extension='.pdf', 
                    directory='', dpi=350, format=TRUE,
                    width = par("din")[1], height = par("din")[2]){
    plt = plts[[plt_name]]
    if (format){
        plt = ggpar(plt, font.title=8, font.subtitle=8, font.caption=8, 
                    font.x=8, font.y=8, font.legend=8,
                    font.tickslab=6, font.family='Arial')    
    }
    filename = file.path(directory,paste0(plt_name,extension))
    save_plot(filename, plt, base_width=width, base_height=height, dpi=dpi, units='cm')
}


save_plots = function(plts, figs_dir){
    
    # model selection
    save_plt(plts, 'model_selection-deps_sorted_vs_std_ctl_neg', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'model_selection-deps_med_vs_std_ctl_neg', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'model_selection-deps_sorted_vs_std_ctl_pos', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'model_selection-deps_med_vs_std_ctl_pos', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'model_selection-lr_pvalue', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'model_selection-lr_pvalue_ctl_pos', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'model_selection-pearson_corr', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'model_selection-pvalue_vs_n_genes', '.pdf', figs_dir, width=5, height=4)
    save_plt(plts, 'model_selection-pvalue_vs_ctl_neg', '.pdf', figs_dir, width=5, height=4)
    save_plt(plts, 'model_selection-pvalue_vs_ctl_pos', '.pdf', figs_dir, width=5, height=4)
    save_plt(plts, 'model_selection-roc_curve', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'model_selection-tpr_vs_fpr', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'model_selection-pvalue_vs_n_events', '.pdf', figs_dir, width=5, height=4)
    save_plt(plts, 'model_selection-pearson_corr_vs_spearman', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'model_selection-selected_vs_event_std', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'model_selection-mutation_gene_count', '.pdf', figs_dir, width=8, height=8)
    save_plt(plts, 'model_selection-mutation_event_count', '.pdf', figs_dir, width=8, height=10)
    save_plt(plts, 'model_selection-mutation_gene_frequency', '.pdf', figs_dir, width=8, height=8)
    save_plt(plts, 'model_selection-mutation_event_frequency', '.pdf', figs_dir, width=8, height=8)
    save_plt(plts, 'model_selection-mutation_event_frequency-by_protein_impact', '.pdf', figs_dir, width=14, height=9)
    save_plt(plts, 'model_selection-protein_impact-counts', '.pdf', figs_dir, width=8, height=8)
    save_plt(plts, 'model_selection-protein_impact-freqs', '.pdf', figs_dir, width=8, height=8)
    save_plt(plts, 'model_selection-protein_impact_clean-counts', '.pdf', figs_dir, width=6, height=6)
    save_plt(plts, 'model_selection-protein_impact_clean-freqs', '.pdf', figs_dir, width=6, height=6)
    save_plt(plts, 'model_selection-enrichment-GO_BP-dotplot', '.pdf', figs_dir, width=12, height=8)
    save_plt(plts, 'model_selection-enrichment-GO_BP-cnetplot', '.pdf', figs_dir, width=20, height=20)
    # save_plt(plts, 'model_selection-enrichment-GO_BP_vs_prot_imp', '.pdf', figs_dir, width=9, height=8)
    save_plt(plts, 'model_selection-indices-violin', '.pdf', figs_dir, width=6, height=6)
    save_plt(plts, 'model_selection-indices-top_pos', '.pdf', figs_dir, width=10, height=6)
    save_plt(plts, 'model_selection-indices-top_neg', '.pdf', figs_dir, width=10, height=6)
    
    
    # interactions
    save_plt(plts, 'interactions-counts', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'interactions-freqs_subcategories', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'interactions-protein_impact-freqs', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'interactions-protein_impact_clean-freqs', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'interactions-category-pca', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'interactions-category-umap', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'interactions-category-upset', '.pdf', figs_dir, width=5, height=5, format=FALSE)
    # save_plt(plts, 'interactions-enrichment-oncogenic_signatures', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'interactions-enrichment-GO_BP', '.pdf', figs_dir, width=10, height=8)
    
    # examples
    save_plt(plts, 'examples-top_overall-scatter', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'examples-top_overall-table', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'examples-top_by_protein_impact-scatter', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'examples-top_by_protein_impact-table', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'examples-sign_events_gene-bar', '.pdf', figs_dir, width=5, height=5)
    
    # tumorigenesis
    save_plt(plts, 'tumorigenesis-scatter', '.pdf', figs_dir, width=5, height=5)
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
    ccle_stats_file = args$ccle_stats_file
    msigdb_dir = args$msigdb_dir
    protein_impact_file = args$protein_impact_file
    gene_mut_freq_file = args$gene_mut_freq_file
    event_mut_freq_file = args$event_mut_freq_file
    spldep_file = args$spldep_file
    rnai_file = args$rnai_file
    possible_interactions_file = args$possible_interactions_file
    cancer_events_file = args$cancer_events_file
    figs_dir = args$figs_dir
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    models = read_tsv(models_file) %>% 
        mutate(event_gene = paste0(EVENT,'_',GENE),
               event_type = gsub('Hsa','',gsub("[^a-zA-Z]", "",EVENT)),
               event_zscore = event_coefficient_mean / event_coefficient_std,
               gene_zscore = gene_coefficient_mean / gene_coefficient_std,
               interaction_zscore = interaction_coefficient_mean / interaction_coefficient_std,
               intercept_zscore = intercept_coefficient_mean / intercept_coefficient_std)
    ccle_stats = read_tsv(ccle_stats_file)
    indices = read_tsv(indices_file) %>%
        filter(index_name %in% c("stemness","mitotic_index"))
    gene_mut_freq = read_tsv(gene_mut_freq_file) %>% 
        dplyr::rename(GENE=Hugo_Symbol) %>% 
        mutate(log_rel_entropy=log2(max_rel_entropy))
    event_mut_freq = read_tsv(event_mut_freq_file)
    ontologies = list(
        "hallmarks" = read.gmt(file.path(msigdb_dir,'h.all.v7.4.symbols.gmt')),
        "oncogenic_signatures" = read.gmt(file.path(msigdb_dir,'c6.all.v7.4.symbols.gmt')),
        "GO_BP" = read.gmt(file.path(msigdb_dir,'c5.go.bp.v7.4.symbols.gmt')),
        "protein_impact" = read_tsv(protein_impact_file) %>%
                            dplyr::rename(EVENT=EventID, term=ONTO) %>%
                            dplyr::select(term,EVENT) %>%
                            mutate(term_clean=gsub(' \\(.*','',term),
                                   term_clean=gsub('ORF disruption upon sequence inclusion','ORF disruption',term_clean),
                                   term_clean=gsub('ORF disruption upon sequence exclusion','ORF disruption',term_clean))
    )
    spldep = read_tsv(spldep_file)
    rnai = read_tsv(rnai_file)
    possible_interactions = read_tsv(possible_interactions_file)
    cancer_events = read_tsv(cancer_events_file)
    
    # add interaction categories, summary stats and correlation with indices to models
    models = models %>% 
        left_join(get_interaction_categories(models, possible_interactions), 
                  by='event_gene') %>%
        mutate(is_selected = pearson_correlation_mean>THRESH_CORR & lr_pvalue<THRESH_LR_PVALUE) %>%
        dplyr::select(-c(event_mean,event_std,gene_mean,gene_std)) %>% 
        left_join(ccle_stats, by=c("EVENT","ENSEMBL","GENE"))
    
    plts = make_plots(models, rnai, spldep, gene_mut_freq, event_mut_freq, ontologies, cancer_events)

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