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
require(enrichplot)

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

# variables
MIN_OBS = 50
THRESH_ZSCORE = 1.96
THRESH_PVALUE = 0.05

PALETTE_VARS = 'uchicago'


# Development
# -----------
# PREP_DIR = file.path(ROOT,'data','prep')
# RESULTS_DIR = file.path(ROOT,'results','model_splicing_dependency')
# models_file = file.path(RESULTS_DIR,'files','models_gene_dependency-EX.tsv.gz')
# msigdb_dir = file.path(ROOT,'data','raw','MSigDB','msigdb_v7.4','msigdb_v7.4_files_to_download_locally','msigdb_v7.4_GMTs')
# protein_impact_file = file.path(ROOT,'data','raw','VastDB','PROT_IMPACT-hg38-v3.tab.gz')
# mut_freq_file = file.path(ROOT,'data','prep','mutation_freq','CCLE.tsv.gz')
# figs_dir = file.path(RESULTS_DIR,'figures','eda_models')


##### FUNCTIONS #####
get_genes_lists = function(models){
    df = models 
    cols_oi = c('event_zscore','gene_zscore',
                'interaction_zscore','intercept_zscore')
    
    genes_lists= df[,c('GENE',cols_oi)] %>% 
        pivot_longer(cols_oi, names_to='model_variable', values_to='zscore') %>%
        mutate(significant = abs(zscore)>THRESH_ZSCORE) %>%
        filter(significant) %>%
        distinct(model_variable,GENE) %>%
        with(., split(GENE, model_variable))
        
    return(genes_lists)
}


get_events_lists = function(models){
    df = models 
    cols_oi = c('event_zscore','gene_zscore',
                'interaction_zscore','intercept_zscore')
    
    events_lists= df[,c('EVENT',cols_oi)] %>% 
        pivot_longer(cols_oi, names_to='model_variable', values_to='zscore') %>%
        mutate(significant = abs(zscore)>THRESH_ZSCORE) %>%
        filter(significant) %>%
        distinct(model_variable,EVENT) %>%
        with(., split(EVENT, model_variable))
        
    return(events_lists)
}


get_universe = function(models){
    df = models
    universe = df[,c('EVENT','GENE')] %>% apply(., 2, unique)
    names(universe) = c('events','genes')
    return(universe)
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


plot_qc = function(models){
    df = models
    
    plts = list()
    
    # p-value distributions
    cols_oi = c('event_pvalue','gene_pvalue',
                'interaction_pvalue','intercept_pvalue')
    
    X = df[,c('EVENT','n_obs','event_std',cols_oi)] %>% 
        pivot_longer(cols_oi, names_to='model_variable', values_to='pvalue')
    
    plts[['qc-pvalues']] = X %>%
        gghistogram(x='pvalue', bins=100, color=NA, 
                    fill='model_variable', palette=PALETTE_VARS) +
        geom_vline(
            data = X %>% group_by(model_variable) %>% 
                summarize(med = median(pvalue, na.rm=TRUE)), 
            mapping = aes(xintercept=med, group=model_variable),
            linetype='dashed'
        ) +
        facet_wrap(~model_variable, scales = 'free') + 
        guides(fill=FALSE) +
        labs(x='p-value', y='Count')
    
    # unseen due to few observations
    plts[['qc-n_obs_vs_pvalues']] = X %>%
        mutate(log10_pvalue = -log10(pvalue)) %>%
        ggplot(aes(x=n_obs, y=log10_pvalue, color=model_variable)) +
        geom_scattermore(pointsize = 2, pixels = c(1000,1000)) +
        facet_wrap(~model_variable, scales = 'free') +
        theme_pubr(border = TRUE) +
        guides(color=FALSE) +
        labs(x='No. Observations', y='-log10(p-value)') +
        geom_hline(yintercept = -log10(0.05), linetype='dashed') +
        geom_vline(xintercept = 50, linetype='dashed') +
        stat_cor(method='spearman') +
        ggpubr::color_palette(palette=PALETTE_VARS)
    
    plts[['qc-n_obs_vs_coefs']] = df %>%
        pivot_longer(matches('coefficient'), 
                     names_to='model_variable', 
                     values_to = 'coefficient') %>%
        ggplot(aes(x=n_obs, y=coefficient, color=model_variable)) +
        geom_scattermore(pointsize = 2, pixels = c(1000,1000)) +
        facet_wrap(~model_variable, scales = 'free') +
        theme_pubr(border = TRUE) +
        guides(color=FALSE) +
        labs(x='No. Observations', y='-log10(p-value)') +
        geom_hline(yintercept = -log10(0.05), linetype='dashed') +
        geom_vline(xintercept = 50, linetype='dashed') +
        stat_cor(method='spearman') +
        ggpubr::color_palette(palette=PALETTE_VARS)
    
    plts[['qc-n_obs_vs_event_std']] = df %>%
        ggplot(aes(x=n_obs, y=event_std)) +
        geom_scattermore(pointsize = 2, pixels = c(1000,1000), color='orange') +
        theme_pubr(border = TRUE) +
        labs(x='No. Observations', y='Event PSI Std.') +
        geom_vline(xintercept = 50, linetype='dashed') +
        stat_cor(method='spearman') +
        geom_smooth(method='lm', linetype='dashed', color='black')
    
    # events per gene vs significant events per gene
    X = df %>% 
        mutate(is_significant = event_pvalue<THRESH_PVALUE |
                                interaction_pvalue<THRESH_PVALUE) %>% 
        group_by(GENE,is_significant) %>% 
        summarize(n = n()) %>% 
        group_by(GENE) %>% 
        mutate(total = sum(n)) %>% 
        filter(is_significant)
    
    plts[['qc-sign_events_gene-bar']] = X %>%
        dplyr::rename(significant=n) %>%
        pivot_longer(c(significant,total)) %>%
        ungroup() %>%
        count(name,value) %>%
        ggbarplot(x='value', y='n', numeric.x.axis=TRUE, 
                  fill='name', color=NA, palette='jco') +
        labs(x='Events per Gene', y='Count', fill='')
    
    plts[['qc-sign_events_gene-scatter']] = X %>%
        ggscatter(
            x='total', y='n', margin.params = list(fill = "lightgray"),
            add='reg.line', add.params=list(linetype='dashed'), conf.int=TRUE
        ) +
        stat_cor(method='spearman') +
        labs(x='Total Events per Gene', y='Significant Events per Gene') +
        geom_text_repel(
            aes(label=GENE),
            X %>% ungroup() %>% slice_max(order_by = abs(n - total), n=20)
        )
    
    return(plts)
}


plot_coefs = function(models){
    df = models %>% filter(n_obs>=MIN_OBS)
    cols_oi = c('event_zscore','gene_zscore',
                'interaction_zscore','intercept_zscore')
    
    plts = list()
    
    # Z-score distributions
    X = df[,c('event_gene',cols_oi)] %>% 
        pivot_longer(cols_oi, names_to='model_variable', values_to='zscore')
    
    plts[['zscores-violins']] = X %>%
        ggviolin(x='model_variable', y='zscore', color=NA, 
                 fill='model_variable', palette=PALETTE_VARS) +
        geom_hline(yintercept=0, linetype='dashed') +
        facet_wrap(~model_variable, scales = 'free') + 
        guides(fill=FALSE) +
        labs(x='Model Variable', y='Count') +
        theme_pubr(border = TRUE) +
        geom_text_repel(
            data = X %>% 
                group_by(model_variable) %>% 
                slice_max(order_by = abs(zscore), n = 5),
            aes(model_variable, zscore, label=event_gene)
        )
    
    # coefficient volcano plots
    X = df %>% 
        dplyr::select(matches('event_gene|coefficient|pvalue')) %>% 
        pivot_longer(cols = -event_gene) %>% 
        separate(name, c("variable","measurement")) %>% 
        pivot_wider(
            id_cols = c('event_gene','variable'), 
            names_from = 'measurement', 
            values_from = 'value') %>%
        mutate(log10_pvalue = -log10(pvalue))
    
    plts[['coefs-volcanos']] = X %>%
        ggplot(aes(x=coefficient, y=log10_pvalue)) +
            geom_scattermore(pixels=c(1000,1000), pointsize = 2, alpha=0.8) +
            facet_wrap(~variable, scales='free') +
            theme_pubr(border=TRUE) +
            geom_hline(yintercept = -log10(0.05), linetype='dashed') +
            geom_text_repel(
                data = X %>% 
                group_by(variable) %>% 
                filter(pvalue<THRESH_PVALUE & !is.infinite(log10_pvalue)) %>% 
                slice_max(order_by = abs(coefficient), n=5),
            aes(x=coefficient, y=log10_pvalue, label=event_gene))
    
    plts[['coefs-volcanos_zoom']] = plts[['coefs-volcanos']] + xlim(-10,10)
    
    # violins of coefficients? barplot of top extreme s coefficients?
    plts[['coefs-violins']] = X %>%
        filter(pvalue < THRESH_PVALUE) %>%
        ggviolin(x='variable', y='coefficient', color=NA, 
                 fill='variable', palette=PALETTE_VARS) +
        geom_hline(yintercept=0, linetype='dashed') +
        facet_wrap(~variable, scales = 'free') + 
        guides(fill=FALSE) +
        labs(x='Model Variable', y='Count') +
        theme_pubr(border = TRUE) +
        geom_text_repel(
            data = X %>% 
                group_by(variable) %>% 
                slice_max(order_by = abs(coefficient), n = 5),
            aes(label=event_gene)
        )
    
    plts[['coefs-barplots_top']] = X %>% 
        filter(pvalue < THRESH_PVALUE) %>%
        group_by(variable) %>% 
        slice_max(order_by = abs(coefficient), n = 15) %>% 
        ggbarplot(x='event_gene', y='coefficient', 
                  fill='variable', color=NA, palette=PALETTE_VARS) + 
        facet_wrap(~variable, scales='free') + 
        guides(fill=FALSE) +
        labs(x='Event & Gene', y='Coefficient') + 
        theme_pubr(x.text.angle = 70, border=TRUE) 
    
    return(plts)
}


plot_mutations = function(models, mut_freq){
    X = models %>% filter(n_obs<MIN_OBS) %>% left_join(mut_freq, by='GENE')
    
    plts = list()
    plts[['mutations-entropy_scatter']] = mut_freq %>% 
        ggplot(aes(x=length, y=entropy)) + 
        geom_scattermore(pointsize = 2, pixels = c(1000,1000)) + 
        theme_pubr() + 
        labs(x='Gene Length', y='Entropy')
    
    plts[['mutations-max_rel_entropy-scatter']] = mut_freq %>% 
        ggplot(aes(x=length, y=max_rel_entropy)) + 
        geom_scattermore(pointsize = 2, pixels = c(1000,1000)) + 
        theme_pubr() + 
        labs(x='Gene Length', y='Max. Rel. Entropy')
    
    plts[['mutations-min_rel_entropy-scatter']] = mut_freq %>% 
        ggplot(aes(x=length, y=min_rel_entropy)) + 
        geom_scattermore(pointsize = 2, pixels = c(1000,1000)) + 
        theme_pubr() + 
        labs(x='Gene Length', y='Min. Rel. Entropy')
    
    plts[['mutations-rel_entropy-mut_effect']] = mut_freq %>%
        ggboxplot(x='Variant_Classification', y='log_rel_entropy') + 
        theme_pubr(x.text.angle = 45) + 
        labs(x='Mutation Effect', y='log2(Rel. Entropy)')
    
    plts[['mutations-rel_mut_freq-mut_effect']] = mut_freq %>%
        ggboxplot(x='Variant_Classification', y='mut_freq_per_kb') + 
        theme_pubr(x.text.angle = 45) + 
        yscale('log10') +
        labs(x='Mutation Effect', y='log10(Mut. Freq. per Kb)')
    
    # do significant models have low mutation rate?
    plts[['mutations-rel_mut_freq_vs_sign_models']] = X %>% 
        mutate(is_significant = event_pvalue<THRESH_PVALUE | 
                                interaction_pvalue<THRESH_PVALUE) %>%
        group_by(GENE, Variant_Classification, mut_freq_per_kb) %>%
        summarize(is_significant = any(is_significant, na.rm=TRUE)) %>%
        drop_na() %>%
        ggplot(aes(x=Variant_Classification, y=mut_freq_per_kb, 
                   group=interaction(Variant_Classification,is_significant))) +
        geom_violin(aes(fill=is_significant), color=FALSE) +
        geom_boxplot(width=0.3, outlier.size=0.5, position=position_dodge(0.9)) +
        stat_compare_means(aes(group=is_significant), method='wilcox.test', label='p.signif') +
        yscale('log10') + 
        fill_palette('lancet') +
        labs(x='Mutation Effect', y='log10(Mut. Freq. per Kb)', fill='Signif. Model') +
        theme_pubr(x.text.angle=70)
    
    return(plts)
}


plot_protein_impact = function(models, protein_impact){
    X = models %>% 
        filter(n_obs<MIN_OBS) %>% 
        mutate(is_significant = event_pvalue<THRESH_PVALUE | 
                                interaction_pvalue<THRESH_PVALUE) %>% 
        left_join(protein_impact, by='EVENT') %>%
        group_by(term, is_significant) %>%
        summarize(n=n()) %>%
        drop_na() %>%
        mutate(freq=n/sum(n))
    
    plts = list()
    plts[['protein_impact-counts']] = X %>%
        ggbarplot(x='term', y='n', fill='is_significant', palette='lancet',
                  position=position_dodge(0.8), color=NA, label=TRUE) +
        labs(x='Protein Impact', y='Count', fill='Signif. Model') +
        theme_pubr(x.text.angle=45)
    
    plts[['protein_impact-ratios']] = X %>%
        ggbarplot(x='term', y='freq', fill='is_significant',
                  palette='lancet', color=NA) +
        labs(x='Protein Impact', y='Proportion', fill='Signif. Model') +
        theme_pubr(x.text.angle=45)
    
    return(plts)
}


plot_enrichments = function(result, 
                            pattern='',
                            palette='lancet', 
                            legend_label='Model Variables'){
    res = new("compareClusterResult", compareClusterResult = result)
    
    # prepare palette
    n = length(unique(result$Cluster))
    palette = get_palette(palette, n)
    
    # plot
    plts = list(
        'enrichment-dotplot' = dotplot(res),
        'enrichment-cnetplot' = cnetplot(res) + 
            scale_fill_manual(values=palette) + labs(fill=legend_label)
    )
    names(plts) = paste0(pattern,'-',names(plts))
    
    return(plts)
}


make_plots = function(models, results_enrich, mut_freq, protein_impact){
    plts = list(
        plot_qc(models),
        plot_coefs(models),
        plot_enrichments(results_enrich[['hallmarks']], 'hallmarks'),
        plot_enrichments(results_enrich[['oncogenic_signatures']], 'oncogenic_signatures'),
        plot_enrichments(results_enrich[['GO_BP']], 'GO_BP'),
        plot_mutations(models, mut_freq),
        plot_protein_impact(models, protein_impact)
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
    # QC
    save_plt(plts, 'qc-pvalues', '.pdf', figs_dir, width=10, height=10)
    save_plt(plts, 'qc-n_obs_vs_pvalues', '.pdf', figs_dir, width=10, height=10)
    save_plt(plts, 'qc-n_obs_vs_coefs', '.pdf', figs_dir, width=10, height=10)
    save_plt(plts, 'qc-n_obs_vs_event_std', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'qc-sign_events_gene-bar', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'qc-sign_events_gene-scatter', '.pdf', figs_dir, width=5, height=5)
    
    # coefficients
    save_plt(plts, 'zscores-violins', '.pdf', figs_dir, width=10, height=10)
    save_plt(plts, 'coefs-volcanos', '.pdf', figs_dir, width=10, height=10)
    save_plt(plts, 'coefs-volcanos_zoom', '.pdf', figs_dir, width=10, height=10)
    save_plt(plts, 'coefs-violins', '.pdf', figs_dir, width=10, height=10)
    save_plt(plts, 'coefs-barplots_top', '.pdf', figs_dir, width=10, height=10)
    
    # enrichments
    ## hallmarks
    save_plt(plts, 'hallmarks-enrichment-dotplot', '.pdf', figs_dir, width=15, height=10)
    save_plt(plts, 'hallmarks-enrichment-cnetplot', '.png', figs_dir, width=20, height=20)
    ## oncogenic signatures
    save_plt(plts, 'oncogenic_signatures-enrichment-dotplot', '.pdf', figs_dir, width=15, height=10)
    save_plt(plts, 'oncogenic_signatures-enrichment-cnetplot', '.png', figs_dir, width=20, height=20)
    ## GO BP
    save_plt(plts, 'GO_BP-enrichment-dotplot', '.pdf', figs_dir, width=15, height=10)
    save_plt(plts, 'GO_BP-enrichment-cnetplot', '.png', figs_dir, width=20, height=20)
    
    # mutations
    save_plt(plts, 'mutations-entropy_scatter', '.pdf', figs_dir, width=6, height=6)
    save_plt(plts, 'mutations-max_rel_entropy-scatter', '.pdf', figs_dir, width=6, height=6)
    save_plt(plts, 'mutations-min_rel_entropy-scatter', '.pdf', figs_dir, width=6, height=6)
    save_plt(plts, 'mutations-rel_entropy-mut_effect', '.pdf', figs_dir, width=8, height=6)
    save_plt(plts, 'mutations-rel_mut_freq-mut_effect', '.pdf', figs_dir, width=8, height=6)
    save_plt(plts, 'mutations-rel_mut_freq_vs_sign_models', '.pdf', figs_dir, width=8, height=6)
    
    # protein impact
    save_plt(plts, 'protein_impact-counts', '.pdf', figs_dir, width=8, height=6)
    save_plt(plts, 'protein_impact-ratios', '.pdf', figs_dir, width=8, height=6)
    
}


make_figdata = function(results_enrich){
    figdata = list(
        'gsea-hallmarks'= results_enrich[['hallmarks']],
        'gsea-oncogenic_signatures' = results_enrich[['oncogenic_signatures']],
        'gsea-GO_BP' = results_enrich[['GO_BP']]
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
                            dplyr::select(term,EVENT)
    )
    
    # run enrichments
    genes_lists = get_genes_lists(models)
    events_lists = get_events_lists(models)
    universe = get_universe(models)
    enrichments = run_enrichments(genes_lists, events_lists, universe, ontologies)
    results_enrich = get_enrichment_result(enrichments)
    
    # plot
    plts = make_plots(models, results_enrich, mut_freq, 
                      protein_impact=ontologies[['protein_impact']])

    # make figdata
    figdata = make_figdata(results_enrich)
    
    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}