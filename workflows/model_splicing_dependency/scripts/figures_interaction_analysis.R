# - Interaction analysis:
#   - frequencies of type of interaction
#   - which variable has the more weight?
#   - are events that interact with gene expression more prone to be related to NMD? 
#      - https://www.nature.com/articles/s41588-019-0555-z?proof=t#Sec10
#      - https://www.sciencedirect.com/science/article/pii/S1097276520307267
#   - models with splicing factors? overlap with known posion exons
#   - moderator analysis:
#      - from those models with significant mRNA levels and PSI 
#        interaction, who acts as a moderator? 
#        Is it complete or partial moderation?

require(tidyverse)
require(ggrepel)
require(ggpubr)
require(cowplot)
require(scattermore)
require(pheatmap)
require(clusterProfiler)

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

MIN_OBS = 50
THRESH_ZSCORE = 1.96
THRESH_PVALUE = 0.05

# Development
# -----------
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# RESULTS_DIR = file.path(ROOT,'results','model_splicing_dependency')
# models_file = file.path(RESULTS_DIR,'files','models_gene_dependency-EX.tsv.gz')
# possible_interactions_file = file.path(ROOT,'support','possible_pairwise_interaction_categories.tsv')
# protein_impact_file = file.path(RAW_DIR,'VastDB','PROT_IMPACT-hg38-v3.tab.gz')
# msigdb_dir = file.path(RAW_DIR,'MSigDB','msigdb_v7.4','msigdb_v7.4_files_to_download_locally','msigdb_v7.4_GMTs')
# figs_dir = file.path(ROOT,'results','model_splicing_dependency','figures','interaction_analysis')


##### FUNCTIONS #####
get_genes_lists = function(models){
    df = models 
    genes_lists = df[,c('GENE','interaction_category')] %>% 
        filter(interaction_category!='no_regulation') %>%
        mutate(interaction_category=factor(interaction_category)) %>% 
        distinct() %>%
        with(., split(GENE, interaction_category))
    return(genes_lists)
}


get_events_lists = function(models){
    df = models 
    events_lists = df[,c('EVENT','interaction_category')] %>% 
        filter(interaction_category!='no_regulation') %>%
        mutate(interaction_category=factor(interaction_category)) %>% 
        distinct() %>%
        with(., split(EVENT, interaction_category))
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


get_interaction_categories = function(models, possible_interactions){
    # study interaction between event inclusion and gene expression
    possible_interactions = possible_interactions %>% 
        mutate(combined=paste0(beta_a,beta_b,beta_ab))
    
    ## interactions
    zscores = models %>% 
        column_to_rownames("event_gene") %>% 
        dplyr::select(paste0(c('event','gene','interaction'),"_zscore"))
    intcats = (abs(zscores) > THRESH_ZSCORE) * sign(zscores)
    intcats[is.na(intcats)] = 0
    intcats = intcats %>% 
        rownames_to_column("event_gene") %>%
        mutate(combined = paste0(event_zscore,gene_zscore,interaction_zscore)) %>%
        left_join(possible_interactions,by='combined') %>%
        dplyr::select(event_gene,category,combined) %>%
        rename(interaction_category=category, interaction_subcategory=combined)
    
    return(intcats)
}


plot_interactions = function(models){
    
    X = models
    
    plts = list()
    
    # across event types, what's the frequency of each type of psi-tpm interaction?
    plts[['interactions-counts_overview']] = X %>% 
        count(event_type, interaction_category) %>%
        ggbarplot(x='interaction_category', y='n', facet.by = 'event_type', 
                  label=TRUE, fill='interaction_category', color=NA, palette='Set1') + 
        guides(fill=FALSE) + 
        yscale('log10', .format=TRUE) + 
        labs(x='Interaction Category', y='Counts') +
        theme_pubr(x.text.angle = 45, border=TRUE)
    
    plts[['interactions-counts_subcategories']] = X %>%
        mutate(interaction_subcategory = gsub('-1','1',interaction_subcategory)) %>%
        group_by(interaction_category, interaction_subcategory) %>%
        summarize(n=n()) %>%
        mutate(freq=n/sum(n)) %>%
        ggbarplot(x='interaction_category', y='freq', fill='interaction_subcategory', 
                  color=NA, palette='Set2') +
        theme_pubr(x.text.angle = 45, legend = 'right') +
        labs(x='Interaction Category', y='Proportion', fill='Interaction\nSubcategory')
    
    # Is some interaction category associated to some protein impact?
    plts[['interactions-rel_protein_impact']] = X %>% 
        group_by(interaction_category,term) %>% 
        summarize(n=n()) %>% 
        mutate(freq= n/sum(n)) %>% 
        drop_na() %>%
        ggbarplot(x='interaction_category', y='freq', fill='term', 
                  color=NA, palette='simpsons') +
        labs(x='Interaction Category', y='Proportion', fill='Protein Impact') +
        theme_pubr(x.text.angle = 45, legend = 'right')
        
    
    plts[['interactions-overall_protein_impact']] = X %>% 
        filter(interaction_category!='no_regulation') %>% 
        group_by(term,interaction_category) %>% 
        summarize(n=n()) %>% 
        ungroup() %>% 
        mutate(freq=n/sum(n)) %>% 
        drop_na() %>% 
        pivot_wider(id_cols='term', 
                    names_from='interaction_category', 
                    values_from='freq') %>% replace(is.na(.), 0) %>% 
        column_to_rownames('term') %>% 
        t() %>% 
        pheatmap(border_color = 'white', silent = TRUE)
    
    return(plts)
}


plot_enrichments = function(result, 
                            pattern='',
                            palette='lancet', 
                            legend_label='Interaction Category'){
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


make_plots = function(models, results_enrich){
    plts = list(
        plot_interactions(models),
        plot_enrichments(results_enrich[['hallmarks']], 'hallmarks'),
        plot_enrichments(results_enrich[['oncogenic_signatures']], 'oncogenic_signatures'),
        plot_enrichments(results_enrich[['GO_BP']], 'GO_BP'),
        plot_enrichments(results_enrich[['protein_impact']], 'protein_impact')
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
    # interactions
    save_plt(plts, 'interactions-counts_overview', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'interactions-counts_subcategories', '.pdf', figs_dir, width=8, height=5)
    save_plt(plts, 'interactions-rel_protein_impact', '.pdf', figs_dir, width=7, height=5)
    save_plt(plts, 'interactions-overall_protein_impact', '.pdf', figs_dir, width=10, height=10)
    
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
    ## protein impact
    save_plt(plts, 'protein_impact-enrichment-dotplot', '.pdf', figs_dir, width=15, height=10)
    save_plt(plts, 'protein_impact-enrichment-cnetplot', '.png', figs_dir, width=20, height=20)
}


main = function(){
    args = getParsedArgs()
    models_file = args$models_file
    possible_interactions_file = args$possible_interactions_file
    protein_impact_file = args$protein_impact_file
    msigdb_dir = args$msigdb_dir
    figs_dir = args$figs_dir
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    possible_interactions = read_tsv(possible_interactions_file)
    models = read_tsv(models_file) %>%
        filter(n_obs>MIN_OBS) %>%
        mutate(event_gene = paste0(EVENT,'_',GENE),
               event_type = gsub('Hsa','',gsub("[^a-zA-Z]", "",EVENT))) 
    ontologies = list(
        "hallmarks" = read.gmt(file.path(msigdb_dir,'h.all.v7.4.symbols.gmt')),
        "oncogenic_signatures" = read.gmt(file.path(msigdb_dir,'c6.all.v7.4.symbols.gmt')),
        "GO_BP" = read.gmt(file.path(msigdb_dir,'c5.go.bp.v7.4.symbols.gmt')),
        "protein_impact" = read_tsv(protein_impact_file) %>%
                            dplyr::rename(EVENT=EventID, term=ONTO) %>%
                            dplyr::select(term,EVENT)
        )
    
    # add protein impact
    models = models %>% left_join(ontologies[['protein_impact']], by='EVENT')
    
    # add interaction categories to models
    intcats = get_interaction_categories(models, possible_interactions)
    models = models %>% left_join(intcats,by='event_gene')
    
    # run enrichments
    genes_lists = get_genes_lists(models)
    events_lists = get_events_lists(models)
    universe = get_universe(models)
    enrichments = run_enrichments(genes_lists, events_lists, universe, ontologies)
    results_enrich = get_enrichment_result(enrichments)
 
    # plot
    plts = make_plots(models, results_enrich)

    # save
    save_plots(plts, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}