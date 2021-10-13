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
require(ComplexHeatmap)
require(clusterProfiler)
require(umap)
require(gridExtra)
require(ggplotify)

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

# variables
THRESH_ZSCORE = 1.96
THRESH_PVALUE = 0.05
THRESH_LR_PVALUE = 0.001

# Development
# -----------
# PREP_DIR = file.path(ROOT,'data','prep')
# RESULTS_DIR = file.path(ROOT,'results','model_splicing_dependency')
# models_file = file.path(RESULTS_DIR,'files','models_gene_dependency-EX.tsv.gz')
# rnai_file = file.path(PREP_DIR,'demeter2','CCLE.tsv.gz')
# spldep_file = file.path(RESULTS_DIR,'files','splicing_dependency_mean-EX.tsv.gz')
# msigdb_dir = file.path(ROOT,'data','raw','MSigDB','msigdb_v7.4','msigdb_v7.4_files_to_download_locally','msigdb_v7.4_GMTs')
# protein_impact_file = file.path(ROOT,'data','raw','VastDB','PROT_IMPACT-hg38-v3.tab.gz')
# mut_freq_file = file.path(ROOT,'data','prep','mutation_freq','CCLE.tsv.gz')
# figs_dir = file.path(RESULTS_DIR,'figures','model_selection')
# possible_interactions_file = file.path(ROOT,'support','possible_pairwise_interaction_categories.tsv')

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


plot_model_selection = function(models, rnai, spldep, mut_freq, ontologies){
    # prep data
    rnai = rnai %>% filter(index%in%models[['GENE']]) %>% column_to_rownames('index')
    spldep = spldep %>% filter(index%in%models[['EVENT']]) %>% column_to_rownames('index')
    
    plts = list()
    
    # as a negative control, we selected 100 uniform gene dependencies
    rnai_stats = data.frame(
            rnai_med = apply(rnai,1,median, na.rm=TRUE),
            rnai_std = apply(rnai,1,sd, na.rm=TRUE),
            n_missing = rowSums(is.na(rnai))
        ) %>% rownames_to_column('GENE')
    
    uniform_genes = rnai_stats %>% 
        slice_min(order_by = rnai_std, n=100) %>% 
        pull(GENE)
    
    plts[['model_selection-deps_sorted_vs_std']] = rnai_stats %>% 
        arrange(-rnai_std) %>% 
        mutate(index=row_number(), 
               is_uniform=GENE %in% uniform_genes) %>%
        ggplot(aes(x=index, y=rnai_std, color=is_uniform)) +
        geom_scattermore(pixels=c(1000,1000), pointsize=4) +
        theme_pubr() +
        ggpubr::color_palette(palette=c('black','orange')) +
        labs(x='Index', y='Gene Std. Demeter2', color='Uniform Dependency')
    
    plts[['model_selection-deps_med_vs_std']] = rnai_stats %>% 
        mutate(is_uniform=GENE %in% uniform_genes) %>%
        ggplot(aes(x=rnai_med, y=rnai_std, color=is_uniform)) +
        geom_scattermore(pixels=c(1000,1000), pointsize=4) +
        theme_pubr() +
        ggpubr::color_palette(palette=c('black','orange')) +
        labs(x='Gene Median Demeter2', y='Gene Std. Demeter2', 
             color='Uniform Dependency')
    
    # we selected the ensemble of models by their ability to rank dependencies
    # using likelihood ratio tests' p-values as thresholds
    plts[['model_selection-lr_pvalue']] = models %>%
        gghistogram(x='lr_pvalue', bins=100, fill='grey', color=NA) +
        geom_vline(xintercept=median(models[['lr_pvalue']], na.rm=TRUE),
                   linetype='dashed') +
        labs(x='LR Test p-value', y='Count')
    
    plts[['model_selection-pearson_corr']] = models %>%
        gghistogram(x='pearson_correlation', bins=100, fill='darkred', color=NA) +
        geom_vline(xintercept=median(models[['pearson_correlation']], na.rm=TRUE),
                   linetype='dashed') +
        labs(x='Pearson Correlation (Test Set)', y='Count')
    
    threshs = c(
        1e-4,5e-4,
        1e-3,5e-3,
        1e-2,5e-2,
        1e-1,5e-1
    )
    samples_oi = intersect(colnames(rnai),colnames(spldep))
    corrs = lapply(threshs, function(thresh){
        models_filtered = models %>% 
            filter(pearson_correlation>0 & lr_pvalue<thresh) 
        
        # from each gene, pick de event-level model that generalizes the worst
        models_selected = models_filtered %>%
            group_by(GENE) %>%
            slice_min(order_by=pearson_correlation, n=1)
        
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
                total_uniforms = sum(genes %in% uniform_genes)
            )
            return(df)
        })
        tmp = do.call(rbind,tmp)
        return(tmp)
    })
    corrs = do.call(rbind,corrs)
    
    plts[['model_selection-pvalue_vs_spearman']] = corrs %>% 
        ggviolin(x='thresh', y='corr', fill='orange', color=NA) + 
        geom_boxplot(fill=NA) +
        labs(x='Thresholds LR Test p-value', y='Spearman Correlation',
             title='Sample Size = 250') +
        theme_pubr(x.text.angle=70)
    
    plts[['model_selection-pvalue_vs_n_genes']] = corrs %>%
        dplyr::select(-one_of(c('id','corr','pvalue'))) %>%
        distinct() %>%
        ggbarplot(x='thresh', y='total_genes', label=TRUE) +
        labs(x='Thresholds LR Test p-value', y='No. Genes Selected') +
        theme_pubr(x.text.angle=70)
    
    plts[['model_selection-pvalue_vs_n_uniforms']] = corrs %>%
        dplyr::select(-one_of(c('id','corr','pvalue'))) %>%
        distinct() %>%
        ggbarplot(x='thresh', y='total_uniforms', label=TRUE) +
        labs(x='Thresholds LR Test p-value', y='No. Uniform Dep. Genes Selected') +
        theme_pubr(x.text.angle=70)
    
    plts[['model_selection-pvalue_vs_n_events']] = corrs %>%
        dplyr::select(-one_of(c('id','corr','pvalue'))) %>%
        distinct() %>%
        ggbarplot(x='thresh', y='total_events', label=TRUE) +
        labs(x='Thresholds LR Test p-value', y='No. Events Selected') +
        theme_pubr(x.text.angle=70)
    
    # with this set of models, we expect to see certain properties
    THRESH_LR_PVALUE = 0.001
    models = models %>%
        drop_na(is_selected)
    
    ## selected models come from variant events
    plts[['model_selection-selected_vs_event_std']] = models %>%
        ggviolin(x='is_selected', y='event_std', color=NA, trim=TRUE,
                 fill='is_selected', palette='lancet') +
        geom_boxplot(fill=NA) +
        stat_compare_means(method='wilcox.test') +
        guides(color='none', fill='none') +
        labs(x='Selected Model', y='Event Std.')
    
    ## selected models are in genes less prone to have deleterious mutations
    X = models %>%
        left_join(mut_freq, by='GENE') %>%
        group_by(Variant_Classification,GENE,mut_freq_per_kb) %>%
        summarize(is_selected = any(is_selected)) %>%
        drop_na() %>%
        ungroup()
    
    plts[['model_selection-mutation_gene_count']] = X %>% 
        count(Variant_Classification, is_selected) %>%
        ggbarplot(x='Variant_Classification', y='n', label=TRUE, palette='lancet',
                  fill='is_selected', color=NA, position=position_dodge(0.9)) + 
        yscale('log10', .format=TRUE) + 
        labs(x='Mutation Effect', y='No. Genes', 
             fill='Selected Model') +
        theme_pubr(x.text.angle=70)
    
    plts[['model_selection-mutation_frequency']] = X %>% 
        ggplot(aes(x=Variant_Classification, y=mut_freq_per_kb, 
                   group=interaction(Variant_Classification,is_selected))) +
        geom_violin(aes(fill=is_selected), color=FALSE) +
        geom_boxplot(width=0.2, outlier.size=0.5, 
                     position=position_dodge(0.9), fill=NA) +
        stat_compare_means(aes(group=is_selected), 
                           method='wilcox.test', label='p.signif') +
        yscale('log10', .format=TRUE) + 
        fill_palette('lancet') +
        labs(x='Mutation Effect', y='log10(Mut. Freq. per Kb)', 
             fill='Selected Model') +
        theme_pubr(x.text.angle=70)
        
    
    # what are the exons/genes selected?
    ## protein impact
    prot_imp = models %>% 
        left_join(ontologies[['protein_impact']], by='EVENT')
    
    plts[['model_selection-protein_impact-counts']] = prot_imp %>%
        filter(is_selected) %>%
        group_by(is_selected,term) %>% 
        summarize(n=n()) %>% 
        arrange(n) %>%
        ggbarplot(x='term', y='n', label=TRUE,
                  fill='is_selected', color=NA, palette='lancet') + 
        theme_pubr(x.text.angle = 45, legend='right') +
        guides(color='none') + 
        labs(x='Protein Impact', y='Count')
    
    plts[['model_selection-protein_impact-freqs']] = prot_imp %>%
        group_by(is_selected,term) %>%
        summarize(n=n()) %>%
        drop_na() %>%
        mutate(freq=n/sum(n)) %>%
        ggbarplot(x='is_selected', y='freq', fill='term', color=NA,
                  palette=get_palette('BrBG', k=length(unique(prot_imp[['term']])))) + 
        theme_pubr(x.text.angle = 45, legend='right') +
        guides(color='none') + 
        labs(x='Selected Model', y='Proportion', fill='Protein Impact')
    
    plts[['model_selection-protein_impact_clean-counts']] = prot_imp %>%
        filter(is_selected) %>%
        group_by(is_selected,term_clean) %>%
        summarize(n=n()) %>%
        arrange(n) %>%
        ggbarplot(x='term_clean', y='n', label=TRUE,
                  fill='is_selected', color=NA, palette='lancet') + 
        theme_pubr(x.text.angle = 45, legend='right') +
        guides(color='none') + 
        labs(x='Protein Impact', y='Count')
    
    plts[['model_selection-protein_impact_clean-freqs']] = prot_imp %>%
        group_by(is_selected,term_clean) %>%
        summarize(n=n()) %>%
        mutate(freq=n/sum(n)) %>%
        ggbarplot(x='is_selected', y='freq', fill='term_clean', color=NA,
                  palette='BrBG') + 
        theme_pubr(x.text.angle = 45, legend='right') +
        guides(color='none') + 
        labs(x='Selected Model', y='Proportion', fill='Protein Impact')
    
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
    plts_enrichment = sapply(enrichment, function(res){
        plt = dotplot(res)
        return(plt)
    }, simplify=FALSE)
    names(plts_enrichment) = sprintf('model_selection-enrichment-%s',
                                     names(plts_enrichment))
    plts = c(plts,plts_enrichment)
    
    
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


plot_interactions = function(models, protein_impact){
    
    X = models %>% 
        mutate(interaction_subcategory_clean = gsub('-1','1',interaction_subcategory)) %>%
        left_join(protein_impact, by='EVENT')
        
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
        group_by(interaction_category,term) %>% 
        summarize(n=n()) %>% 
        drop_na() %>%
        mutate(freq= n/sum(n)) %>% 
        ggbarplot(x='interaction_category', y='freq', fill='term', 
                  color=NA, palette=get_palette('BrBG',length(unique(X[['term']])))) +
        labs(x='Interaction Category', y='Proportion', fill='Protein Impact') +
        theme_pubr(x.text.angle = 45, legend = 'right')
    
    plts[['interactions-protein_impact_clean-freqs']] = X %>% 
        filter(is_selected) %>%
        group_by(interaction_category,term_clean) %>% 
        summarize(n=n()) %>% 
        drop_na() %>%
        mutate(freq= n/sum(n)) %>% 
        ggbarplot(x='interaction_category', y='freq', fill='term_clean', 
                  color=NA, palette='BrBG') +
        labs(x='Interaction Category', y='Proportion', fill='Protein Impact') +
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
        arrange(lr_pvalue) %>%
        mutate(index=row_number()) %>%
        left_join(protein_impact, by='EVENT')
    
    plts = list()
    
    cols_oi = c('EVENT','GENE','ENSEMBL','interaction_category',
                'interaction_subcategory','term','event_type','lr_pvalue')
    
    plts[['examples-top_overall-scatter']] = X %>%
        ggscatter(x='index', y='lr_pvalue') +
        yscale('log10', .format=TRUE) +
        geom_text_repel(
            aes(label=event_gene),
            X %>%
            slice_min(order_by=lr_pvalue, n=25),
            max.overlaps = 50
        ) +
        labs(x='Index', y='log10(LR Test p-value)', 
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


make_plots = function(models, rnai, spldep, mut_freq, ontologies){
    plts = list(
        plot_model_selection(models, rnai, spldep, mut_freq, ontologies),
        plot_interactions(models, protein_impact=ontologies[['protein_impact']]),
        plot_examples(models, protein_impact=ontologies[['protein_impact']])
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
    
    # model selection
    save_plt(plts, 'model_selection-deps_sorted_vs_std', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'model_selection-deps_med_vs_std', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'model_selection-lr_pvalue', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'model_selection-pearson_corr', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'model_selection-pvalue_vs_spearman', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'model_selection-pvalue_vs_n_genes', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'model_selection-pvalue_vs_n_uniforms', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'model_selection-pvalue_vs_n_events', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'model_selection-selected_vs_event_std', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'model_selection-mutation_frequency', '.pdf', figs_dir, width=8, height=8)
    save_plt(plts, 'model_selection-protein_impact-counts', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'model_selection-protein_impact-freqs', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'model_selection-protein_impact_clean-counts', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'model_selection-protein_impact_clean-freqs', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'model_selection-enrichment-GO_BP', '.pdf', figs_dir, width=9, height=8)
    
    # interactions
    save_plt(plts, 'interactions-counts', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'interactions-freqs_subcategories', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'interactions-protein_impact-freqs', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'interactions-protein_impact_clean-freqs', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'interactions-category-pca', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'interactions-category-umap', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'interactions-category-upset', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'interactions-enrichment-oncogenic_signatures', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'interactions-enrichment-GO_BP', '.pdf', figs_dir, width=5, height=5)
    
    # examples
    save_plt(plts, 'examples-top_overall-scatter', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'examples-top_overall-table', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'examples-top_by_protein_impact-scatter', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'examples-top_by_protein_impact-table', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'examples-sign_events_gene-bar', '.pdf', figs_dir, width=5, height=5)
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
    msigdb_dir = args$msigdb_dir
    protein_impact_file = args$protein_impact_file
    mut_freq_file = args$mut_freq_file
    figs_dir = args$figs_dir
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    models = read_tsv(models_file) %>% 
        mutate(event_gene = paste0(EVENT,'_',GENE),
               event_type = gsub('Hsa','',gsub("[^a-zA-Z]", "",EVENT)))
    mut_freq = read_tsv(mut_freq_file) %>% 
        dplyr::rename(GENE=Hugo_Symbol) %>% 
        mutate(log_rel_entropy=log2(max_rel_entropy))
    ontologies = list(
        "hallmarks" = read.gmt(file.path(msigdb_dir,'h.all.v7.4.symbols.gmt')),
        "oncogenic_signatures" = read.gmt(file.path(msigdb_dir,'c6.all.v7.4.symbols.gmt')),
        "GO_BP" = read.gmt(file.path(msigdb_dir,'c5.go.bp.v7.4.symbols.gmt')),
        "protein_impact" = read_tsv(protein_impact_file) %>%
                            dplyr::rename(EVENT=EventID, term=ONTO) %>%
                            dplyr::select(term,EVENT) %>%
                            mutate(term_clean=gsub(' \\(.*','',term))
    )
    spldep = read_tsv(spldep_file)
    rnai = read_tsv(rnai_file)
    possible_interactions = read_tsv(possible_interactions_file)
    
    # add interaction categories to models
    models = models %>% 
        left_join(get_interaction_categories(models, possible_interactions), 
                  by='event_gene') %>%
        mutate(is_selected = pearson_correlation>0 & lr_pvalue<THRESH_LR_PVALUE)
    
    plts = make_plots(models, rnai, spldep, mut_freq, ontologies)

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