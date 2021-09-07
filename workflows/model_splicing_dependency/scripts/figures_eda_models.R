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
require(pheatmap)
require(scattermore)
require(ggrepel)

require(clusterProfiler)
require(org.Hs.eg.db)

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

# variables
MIN_OBS = 50
THRESH_ZSCORE = 1.96
ORGDB = org.Hs.eg.db

METADATA_OI = c(
    'sex',
    'source',
    'culture_type',
    'culture_medium',
    'sample_collection_site',
    'primary_or_metastasis',
    'age',
    'cancer_type',
    'mitotic_index',
    'stemness'
)

METADATA_OI_CONT = c(
    'age',
    'mitotic_index',
    'stemness'
)

METADATA_PALETTES = list(
    'sex' = 'lancet',
    'source' = 'default',
    'culture_type' = 'Set2',
    'culture_medium' = 'Dark2',
    'sample_collection_site' = 'simpsons',
    'primary_or_metastasis' = 'npg',
    'age' = 'Greens',
    'cancer_type' = 'Paired',
    'mitotic_index' = 'Reds',
    'stemness' = 'Blues'
)

# Development
# -----------
# PREP_DIR = file.path(ROOT,'data','prep')
# RESULTS_DIR = file.path(ROOT,'results','model_splicing_dependency')
# psi_file = file.path(PREP_DIR,'event_psi','CCLE-EX.tsv.gz')
# genexpr_file = file.path(PREP_DIR,'genexpr_tpm','CCLE.tsv.gz')
# models_file = file.path(RESULTS_DIR,'files','models_gene_dependency-EX.tsv.gz')
# rnai_file = file.path(PREP_DIR,'demeter2','CCLE.tsv.gz')
# embedded_dependency_file = file.path(RESULTS_DIR,'files','embedded_splicing_dependency-EX.tsv.gz')
# metadata_file = file.path(PREP_DIR,'metadata','CCLE.tsv')
# indices_file = file.path(PREP_DIR,'transcriptome_indices','CCLE.tsv.gz')
# genes_splicing_factors_file = file.path(ROOT,'support','splicing_factors-symbols.txt')
# genes_cancer_file = file.path()
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


run_enrichment = function(genes, universe){
    enrichment = enrichGO(genes, keyType = "SYMBOL", universe=universe,
                          ont = 'BP', OrgDb = ORGDB)
    return(enrichment)
}


run_enrichments = function(genes_lists, universe){
    enrichments = sapply(genes_lists, function(genes){
            run_enrichment(genes, universe)
        }, simplify=FALSE)
    return(enrichments)
}


plot_qc = function(models){
    df = models
    
    plts = list()
    
    # p-value distributions
    cols_oi = c('event_pvalue','gene_pvalue',
                'interaction_pvalue','intercept_pvalue')
    
    X = df[,cols_oi] %>% 
        pivot_longer(everything(), names_to='model_variable', values_to='pvalue')
    
    plts[['qc-pvalues']] = X %>%
        gghistogram(x='pvalue', bins=100, color=NA, 
                    fill='model_variable', palette='uchicago') +
        geom_vline(
            data = X %>% group_by(model_variable) %>% 
                summarize(med = median(pvalue, na.rm=TRUE)), 
            mapping = aes(xintercept=med, group=model_variable),
            linetype='dashed'
        ) +
        facet_wrap(~model_variable, scales = 'free') + 
        guides(fill=FALSE) +
        labs(x='p-value', y='Count')
    
    # Scoring
    X = df[,c('pearson_coefficient','spearman_coefficient')] %>% 
        pivot_longer(everything(), names_to='model_variable', values_to='correlation')
    
    plts[['qc-correlations-distrs']] = X %>%
        gghistogram(x='correlation', bins=100, color=NA, 
                    fill='model_variable', palette='uchicago') +
        geom_vline(
            data = X %>% group_by(model_variable) %>% 
                summarize(med = median(correlation, na.rm=TRUE)), 
            mapping = aes(xintercept=med, group=model_variable),
            linetype='dashed'
        ) +
        facet_wrap(~model_variable, scales = 'free') + 
        guides(fill=FALSE) +
        labs(x='Correlation Coef.', y='Count')
    
    plts[['qc-correlations-scatter']] = df %>%
        ggplot(aes(x=spearman_coefficient, y=pearson_coefficient)) +
        geom_scattermore() +
        stat_cor(method='pearson') + 
        geom_abline(intercept = 0, slope = 1, color='orange', linetype='dashed') +
        theme_pubr() +
        labs(x='Spearman Coef.', y='Pearson Coef.')
    
    return(plts)
}


plot_zscores = function(models){
    df = models 
    cols_oi = c('event_zscore','gene_zscore',
                'interaction_zscore','intercept_zscore')
    
    plts = list()
    
    # Z-score distributions
    X = df[,c('event_gene',cols_oi)] %>% 
        pivot_longer(cols_oi, names_to='model_variable', values_to='zscore')
    
    plts[['zscores-violins']] = X %>%
        ggviolin(x='model_variable', y='zscore', color=NA, 
                 fill='model_variable', palette='uchicago') +
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
    
    # significant events per gene
    X = df[,c('GENE',cols_oi)] %>% 
        pivot_longer(cols_oi, names_to='model_variable', values_to='zscore') %>%
        mutate(significant = abs(zscore)>THRESH_ZSCORE) %>%
        drop_na() %>%
        count(GENE,model_variable,significant)
    
    plts[['zscores-sign_events_per_gene-counts']] = X %>%
        filter(significant) %>%
        gghistogram(x='n', color=NA, facet.by = 'model_variable',
                    fill='model_variable', palette='uchicago') +
        guides(fill=FALSE) +
        labs(x='No. Significant Events per Gene', y='Count')
    
    
    plts[['zscores-sign_events_per_gene-percs']] = X %>%
        group_by(GENE, model_variable) %>%
        summarize(perc = n/sum(n)) %>%
        gghistogram(x='perc', color=NA, facet.by = 'model_variable',
                    fill='model_variable', palette='uchicago') +
        guides(fill=FALSE) +
        labs(x='% Significant Events per Gene', y='Count')
    
    
    return(plts)
}


get_enrichment_result = function(enrich_list, thresh=0.05){
    ## groups are extracted from names
    groups = names(enrich_list)
    result = lapply(groups, function(group){
        result = enrich_list[[group]]
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
    return(result)
}


plot_enrichments = function(enrichments, 
                            palette='uchicago', 
                            legend_label='Model Variables'){
    result = get_enrichment_result(enrichments)
    res = new("compareClusterResult", compareClusterResult = result)
    
    # prepare palette
    n = length(unique(result$Cluster))
    palette = get_palette(palette, n)
    
    # plot
    plts = list(
        'enrichment-BP-dotplot' = dotplot(res) + labs(x=legend_label),
        'enrichment-BP-cnetplot' = cnetplot(res) + 
            scale_fill_manual(values=palette) + labs(fill=legend_label)
    )
    
    return(plts)
}


plot_extreme_examples = function(models, psi){
    # interesting events
}


plot_genes_oi = function(genes_oi, key, 
                         models, psi, genexpr, rnai){
    genes = genes_oi[[key]]
    events = models %>% filter(GENE %in% genes) %>% pull(EVENT)
    
    plts = list()
    
    ## model zscores
    cols_oi = c('event_zscore','gene_zscore',
                'interaction_zscore','intercept_zscore')
    X = models %>%
        filter(GENE %in% genes) %>%
        dplyr::select(event_gene, cols_oi) %>%
        pivot_longer(cols_oi, names_to='model_variable', values_to='zscore')
    
    plts[['zscores-violins']] =  X %>%
        ggviolin(x='model_variable', y='zscore', color=NA, 
                 fill='model_variable', palette='uchicago') +
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
    
    ## genexpr
    plts[['genexpr-heatmap']] = genexpr %>%
        filter(index %in% genes) %>%
        mutate_at(vars(-index), scale) %>%
        column_to_rownames('index') %>%
        pheatmap(show_colnames = FALSE, silent=TRUE)
    
    ## splicing
    X = psi %>%
        filter(EVENT %in% events)
    X[is.na(X)] = 0
    plts[['psi-heatmap']] = X %>%
        column_to_rownames('EVENT') %>%
        pheatmap(show_colnames = FALSE, silent=TRUE)
    
    ## gene dependency
    X = rnai %>%
        filter(index %in% genes)
    X[is.na(X)] = 0
    plts[['rnai-heatmap']] = X %>%
        filter(index %in% genes) %>%
        column_to_rownames('index') %>%
        pheatmap(show_colnames = FALSE, color = rev(get_palette('Reds',15)), silent=TRUE)
    
    names(plts) = paste0(key,'-',names(plts))
    return(plts)
}


plot_splicing_dependency = function(df, metadata, pattern, figtitle){
    X = df %>% left_join(metadata, by='index')
    
    plts = list()
    plts[['pca']] = X %>% 
        ggplot(aes_string(x='PC0', y='PC1')) + 
        geom_scattermore(pixels=c(1000,1000), pointsize = 2, alpha=0.8) +
        theme_pubr() +
        labs(title = figtitle)
    
    plts[['umap']] = X %>% 
        ggplot(aes_string(x='UMAP0', y='UMAP1')) + 
        geom_scattermore(pixels=c(1000,1000), pointsize = 2, alpha=0.8) +
        theme_pubr() +
        labs(title = figtitle)
    
    plts[['umap_clusters']] = X %>%
        ggplot(aes_string(x='UMAP0', y='UMAP1', color='leiden_labels')) + 
        geom_scattermore(pixels=c(1000,1000), pointsize = 2, alpha=0.8) +
        theme_pubr() +
        labs(title = figtitle)
    plts[['umap_clusters']] = set_palette(plts[['umap_clusters']], palette = 'Paired')
    
    
    # with metadata information    
    tmp = sapply(METADATA_OI, function(col_oi){
        # make plot
        plt = X %>%
            ggplot(aes_string(x='PC0', y='PC1', color=col_oi)) + 
            geom_scattermore(pixels = c(1000,1000), pointsize = 2, alpha=0.8) +
            theme_pubr() +
            labs(title = figtitle)
        # set palette
        if (col_oi %in% METADATA_OI_CONT){
            palette = METADATA_PALETTES[[col_oi]]
            plt = plt + gradient_color(palette)
        }else{
            palette = get_palette(METADATA_PALETTES[[col_oi]], 
                              length(unique(metadata[[col_oi]])))
            plt = set_palette(plt, palette=palette)
        }
        return(plt)
    }, simplify=FALSE)
    names(tmp) = paste0('pca-',names(tmp))
    plts = c(plts,tmp)
    
    tmp = sapply(METADATA_OI, function(col_oi){
        # make plot
        plt = X %>%
            ggplot(aes_string(x='UMAP0', y='UMAP1', color=col_oi)) + 
            geom_scattermore(pixels = c(1000,1000), pointsize = 2, alpha=0.8) +
            theme_pubr() +
            labs(title = figtitle)
        # set palette
        if (col_oi %in% METADATA_OI_CONT){
            palette = METADATA_PALETTES[[col_oi]]
            plt = plt + gradient_color(palette)
        }else{
            palette = get_palette(METADATA_PALETTES[[col_oi]], 
                              length(unique(metadata[[col_oi]])))
            plt = set_palette(plt, palette=palette)
        }
        return(plt)
    }, simplify=FALSE)
    names(tmp) = paste0('umap-',names(tmp))
    plts = c(plts,tmp)
    
    # prepend pattern name
    names(plts) = paste0(pattern,'-',names(plts))
    
    return(plts)
}


make_plots = function(models, psi, genexpr, rnai, 
                      enrichments, genes_oi, 
                      embedded_dependency, metadata){
    plts = list(
        plot_qc(models),
        plot_zscores(models),
        plot_enrichments(enrichments),
        #plot_extreme_examples(models, psi, genexpr, rnai),
        plot_genes_oi(genes_oi, 'splicing_factors', models, psi, genexpr, rnai),
        # plot_genes_oi(genes_oi, 'cancer_genes', models, psi, genexpr, rnai),
        plot_splicing_dependency(embedded_dependency, metadata, 'splicing_dependency', 'Splicing Dependency')
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
    save_plt(plts, 'qc-correlations-distrs', '.pdf', figs_dir, width=10, height=5)
    save_plt(plts, 'qc-correlations-scatter', '.pdf', figs_dir, width=5, height=5)
    
    # zscores
    save_plt(plts, 'zscores-violins', '.pdf', figs_dir, width=10, height=10)
    save_plt(plts, 'zscores-sign_events_per_gene-counts', '.pdf', figs_dir, width=10, height=10)
    save_plt(plts, 'zscores-sign_events_per_gene-percs', '.pdf', figs_dir, width=10, height=10)
    
    # enrichment
    save_plt(plts, 'enrichment-BP-dotplot', '.pdf', figs_dir, width=15, height=10)
    save_plt(plts, 'enrichment-BP-cnetplot', '.png', figs_dir, width=20, height=20)
    
    # genes of interest
    ## splicing_factors
    save_plt(plts, 'splicing_factors-zscores-violins', '.pdf', figs_dir, width=10, height=10)
    save_plt(plts, 'splicing_factors-psi-heatmap', '.png', figs_dir, width=10, height=15)
    save_plt(plts, 'splicing_factors-rnai-heatmap', '.png', figs_dir, width=10, height=15)
    save_plt(plts, 'splicing_factors-genexpr-heatmap', '.png', figs_dir, width=10, height=15)
    
    # splicing dependency
    lapply(grep('splicing_dependency',names(plts), value = TRUE), function(plt_oi){
        save_plt(plts, plt_oi, '.pdf', figs_dir, width=10, height=10)
    })
    
}


make_figdata = function(enrichments){
    figdata = list(
        'gsea_models'= sapply(enrichments, as.data.frame, simplify=FALSE)
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
    psi_file = args$psi_file
    genexpr_file = args$genexpr_file
    rnai_file = args$rnai_file
    embedded_dependency_file = args$embedded_dependency_file
    genes_splicing_factors_file = args$genes_splicing_factors_file
    metadata_file = args$metadata_file
    indices_file = args$indices_file
    figs_dir = args$figs_dir
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    models = read_tsv(models_file) %>% 
        filter(n_obs>=MIN_OBS) %>%
        mutate(event_gene = paste0(EVENT,'_',GENE),
               event_type = gsub('Hsa','',gsub("[^a-zA-Z]", "",EVENT)))
    psi = read_tsv(psi_file)
    genexpr = read_tsv(genexpr_file)
    rnai = read_tsv(rnai_file)
    embedded_dependency = read_tsv(embedded_dependency_file) %>% 
        mutate(leiden_labels=as.factor(leiden_labels))
    indices = read_tsv(indices_file)
    metadata = read_tsv(metadata_file) %>%
        dplyr::rename(index=DepMap_ID) %>% 
        left_join(indices, by='index')
    genes_splicing_factors = readLines(genes_splicing_factors_file)
    
    # prepare genes oi
    genes_oi = list(
        'splicing_factors' = genes_splicing_factors
    )
    
    # run enrichments
    genes_lists = get_genes_lists(models)
    universe = unique(models[['GENE']])
    enrichments = run_enrichments(genes_lists, universe)
    
    # plot
    plts = make_plots(models, psi, genexpr, rnai, 
                      enrichments, genes_oi, 
                      embedded_dependency, metadata)

    # make figdata
    figdata = make_figdata(enrichments)
    
    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}