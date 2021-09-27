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