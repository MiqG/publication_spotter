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
require(scattermore)
require(cowplot)
require(ggrepel)
require(extrafont)
require(cluster)
require(ComplexHeatmap)
require(ggplotify)
require(gridExtra)
require(clusterProfiler)

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

# variables
THRESH_FDR = 0.1
THRESH_PVALUE = 0.05
RANDOM_SEED = 1234

# Development
# -----------
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# RESULTS_DIR = file.path(ROOT,'results','splicing_dependency_drugs')
# MODELS_DIR = file.path(ROOT,'results','model_splicing_dependency')
# models_file = file.path(RESULTS_DIR,'files','model_summaries_drug_response-EX.tsv.gz')
# drug_targets_file = file.path(RAW_DIR,'GDSC','screened_compunds_rel_8.2.csv')
# figs_dir = file.path(RESULTS_DIR,'figures','model_drug_screens')
# embedding_file = file.path(RESULTS_DIR,'files','embedded_drug_associations-EX.tsv.gz')
# estimated_response_file = file.path(RESULTS_DIR,'files','estimated_drug_response_by_drug-EX.tsv.gz')
# drug_screen_file = file.path(RAW_DIR,'DepMap','gdsc','sanger-dose-response.csv')
# msigdb_dir = file.path(ROOT,'data','raw','MSigDB','msigdb_v7.4','msigdb_v7.4_files_to_download_locally','msigdb_v7.4_GMTs')
# clusters_file = file.path(RESULTS_DIR,'files','cluster_estimated_drug_response-merged-EX.tsv.gz')
# metadata_file = file.path(PREP_DIR,'metadata','CCLE.tsv.gz')
# spldep_ccle_file = file.path(MODELS_DIR,'files','splicing_dependency-EX','mean.tsv.gz')


##### FUNCTIONS #####
get_sets = function(df, set_names, set_values){
    sets = df[,c(set_values,set_names)] %>%
        distinct() %>%
        with(., split(get(set_values),get(set_names)))
    sets = sapply(sets, unique, simplufy=FALSE)
    return(sets)
}


run_enrichment = function(genes, universe, ontologies){
    enrichments = list()
    enrichments[['hallmarks']] = enricher(genes, TERM2GENE=ontologies[['hallmarks']], universe=universe)
    enrichments[['GO_BP']] = enricher(genes, TERM2GENE=ontologies[['GO_BP']], universe=universe)
    enrichments[['oncogenic_signatures']] = enricher(genes, TERM2GENE=ontologies[['oncogenic_signatures']], universe=universe)
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


plot_associations = function(models, drug_targets, embedding, ontologies, rankings){
    top_n = 25
    X = models %>% 
        left_join(drug_targets %>% distinct(DRUG_ID,DRUG_NAME), by="DRUG_ID")
    
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
        count(drug_screen, DRUG_NAME) %>% 
        group_by(DRUG_NAME) %>%
        mutate(total=sum(n)) %>%
        arrange(total) %>%
        ggbarplot(x='DRUG_NAME', y='n', fill='drug_screen', palette='jco', color=NA) + 
        labs(x='Drug', y='N. Significant Associations')
    
    plts[['associations-top_drug_counts']] = X %>% 
        filter(lr_padj < THRESH_FDR) %>% 
        count(drug_screen, DRUG_NAME) %>% 
        group_by(DRUG_NAME) %>%
        mutate(total=sum(n)) %>%
        arrange(-total) %>%
        ungroup() %>%
        mutate(index=as.numeric(factor(DRUG_NAME, levels=unique(DRUG_NAME)))) %>%
        filter(index <= top_n) %>%
        ggbarplot(x='DRUG_NAME', y='n', fill='drug_screen', palette='jco', color=NA) + 
        labs(x='Drug', y='N. Significant Associations', 
             title=sprintf('Top %s', top_n)) +
        coord_flip()
    
    plts[['associations-spldep_counts']] = X %>% 
        filter(lr_padj < THRESH_FDR) %>% 
        count(drug_screen, event_gene) %>% 
        group_by(event_gene) %>%
        mutate(total=sum(n)) %>%
        arrange(total) %>%
        ggbarplot(x='event_gene', y='n', fill='drug_screen', palette='jco', color=NA) + 
        labs(x='Event & Gene', y='N. Significant Associations')
    
    plts[['associations-top_spldep_counts']] = X %>% 
        filter(lr_padj < THRESH_FDR) %>% 
        count(drug_screen, event_gene) %>% 
        group_by(event_gene) %>%
        mutate(total=sum(n)) %>%
        arrange(-total) %>%
        ungroup() %>%
        mutate(index=as.numeric(factor(event_gene, levels=unique(event_gene)))) %>%
        filter(index <= top_n) %>%
        ggbarplot(x='event_gene', y='n', fill='drug_screen', palette='jco', color=NA) + 
        labs(x='Event & Gene', y='N. Significant Associations', 
             title=sprintf('Top %s', top_n)) +
        coord_flip()
    
    # is there any drug that is significantly associated to the splicing
    # dependency of an event of a gene that it targets?
    found_targets = X %>% 
        left_join(drug_targets, by=c("DRUG_ID","DRUG_NAME","GENE"="TARGET")) %>%
        drop_na(TARGET_PATHWAY) %>%
        mutate(is_target=TRUE)
    drugs_oi = found_targets %>% filter(lr_padj < THRESH_FDR) %>% pull(DRUG_NAME) %>% unique()
    events_oi = found_targets %>% filter(lr_padj < THRESH_FDR) %>% pull(EVENT) %>% unique()
    
    plts[['associations-top_found_targets']] = X %>% 
        filter(lr_padj<THRESH_FDR & DRUG_NAME%in%drugs_oi) %>% 
        # add is_target only to some genes
        left_join(found_targets %>% distinct(DRUG_ID,GENE,is_target), 
                  by=c('DRUG_ID','GENE')) %>%
        # add drug pathway info to all entries
        left_join(found_targets %>% distinct(DRUG_ID,TARGET_PATHWAY), 
                  by="DRUG_ID") %>%
        ungroup() %>%
        group_by(DRUG_NAME) %>% 
        slice_max(order_by=spldep_coefficient, n=15) %>%
        mutate(is_target=replace_na(is_target,FALSE), 
               drug_name_clean=sprintf('%s | %s',DRUG_NAME,TARGET_PATHWAY), 
               name=reorder_within(event_gene, spldep_coefficient, DRUG_NAME)) %>%
        ggplot(aes(x=name, y=spldep_coefficient)) +
        geom_col(aes(fill=is_target, color=drug_screen, group=drug_screen), 
                     position=position_dodge(1)) +
        color_palette(c("white","black")) + 
        fill_palette("Set2") +    
        facet_wrap(~drug_name_clean, scales='free', ncol=4) + 
        coord_flip() +
        scale_x_reordered() +     
        labs(x='Drug', y='Effect Size', color='Drug Screen', fill='Is Drug Target',
             title=sprintf('Top 15 Significant Effect Sizes | %s out of %s', 
                           found_targets %>% filter(lr_padj<THRESH_FDR) %>% 
                               distinct(DRUG_NAME) %>% nrow(),
                           found_targets %>% distinct(DRUG_NAME) %>% nrow())) + 
        theme_pubr(border=TRUE) +
        theme(strip.text.x = element_text(size=6, family='Arial'))
    
    # do significantly associated events in gene targets rank to the top of
    # the association effect size?
    ranking_real = X %>% 
        filter(lr_padj<THRESH_FDR & DRUG_NAME%in%drugs_oi) %>%
        distinct(DRUG_NAME,EVENT,event_gene,spldep_coefficient) %>%
        group_by(DRUG_NAME) %>%
        arrange(-abs(spldep_coefficient)) %>%
        mutate(ranking = row_number() / n()) %>%
        filter(EVENT %in% events_oi) %>%
        ungroup() %>%
        arrange(-ranking) %>%
        mutate(ranking_type="real", 
               ranking_cdf = cumsum(ranking) / sum(ranking),
               index = row_number())
    
    set.seed(RANDOM_SEED)
    ranking_random = X %>% 
        filter(lr_padj<THRESH_FDR & DRUG_NAME%in%drugs_oi) %>%
        distinct(DRUG_NAME,EVENT,event_gene,spldep_coefficient) %>%
        group_by(DRUG_NAME) %>%
        arrange(-abs(spldep_coefficient)) %>%
        mutate(EVENT = sample(EVENT)) %>% # randomize
        mutate(ranking = row_number() / n()) %>%
        filter(EVENT %in% events_oi) %>%
        ungroup() %>%
        arrange(-ranking) %>%
        mutate(ranking_type="random", 
               ranking_cdf = cumsum(ranking) / sum(ranking),
               index = row_number())
    
    ranking = rbind(ranking_real, ranking_random)
    
    plts[['associations-target_ranking-cdf']] = ranking %>%
        ggscatter(x="index", y="ranking_cdf", 
                  color="ranking_type", palette=c("grey","darkred")) +
        labs(x='Index', y='Ranking CDF', color='Ranking Type',
             title=sprintf('n = %s', nrow(ranking_real)))
    
    plts[['associations-target_ranking-vio']] = ranking %>%
        ggviolin(x="ranking_type", y="ranking", color=NA, 
                 fill="ranking_type", palette=c("darkred","grey"), trim=TRUE) +
        geom_boxplot(width=0.1) + 
        stat_compare_means(method="wilcox.test") +
        guides(fill="none") +
        labs(x='Ranking Type', y='Ranking Ratio', 
             title=sprintf('n = %s', nrow(ranking_real)))
    
    # best associations of drugs without targets
    drugs_oi = drug_targets %>% filter(is.na(TARGET)) %>% pull(DRUG_ID) %>% unique()
    
    plts[['associations-top_no_targets']] = X %>% 
        filter(lr_padj<THRESH_FDR & DRUG_ID%in%drugs_oi) %>%
        group_by(DRUG_NAME) %>% 
        slice_max(order_by=spldep_coefficient, n=15) %>% 
        left_join(drug_targets %>% distinct(DRUG_ID,TARGET_PATHWAY), by="DRUG_ID") %>%
        mutate(is_target=FALSE, 
               drug_name_clean=sprintf('%s | %s',DRUG_NAME,TARGET_PATHWAY), 
               name=reorder_within(event_gene, spldep_coefficient, DRUG_NAME)) %>% 
        ggplot(aes(x=name, y=spldep_coefficient)) +
        geom_col(aes(fill=is_target, color=drug_screen, group=drug_screen), 
                     position=position_dodge(1)) +
        color_palette(c("white","black")) + 
        fill_palette("Set2") +
        facet_wrap(~drug_name_clean, scales='free', ncol=4) +
        scale_x_reordered() +     
        labs(x='Event & Gene', y='Effect Size', fill='Is Drug Target', 
             title='Top 15 Significant Effect Sizes') + 
        coord_flip() +
        theme_pubr(border=TRUE) +
        theme(strip.text.x = element_text(size=6, family='Arial'))
    
    # are association profiles informative of drug mechanism of action?
    plts[['associations-target_pathway-counts']] = drug_targets %>% 
        distinct(DRUG_ID,TARGET_PATHWAY) %>%
        count(TARGET_PATHWAY) %>%
        # filter(!(TARGET_PATHWAY %in% c("Other", "Other, kinases", "Unclassified"))) %>%
        ggbarplot(x='TARGET_PATHWAY', y="n",
                  fill="TARGET_PATHWAY", color=NA, 
                  palette=get_palette("jco", 24)) +
        labs(x='Target Pathway', y='Count') +
        coord_flip() +
        theme_pubr(legend="none")
    
    plts[['associations-target_pathway-umap']] = embedding %>%
        left_join(drug_targets %>% distinct(DRUG_ID,TARGET_PATHWAY), 
                  by=c("index"="DRUG_ID")) %>%
        # filter(!(TARGET_PATHWAY %in% c("Other", "Other, kinases", "Unclassified"))) %>%
        ggplot(aes(x=UMAP0, y=UMAP1, color=TARGET_PATHWAY)) +
        geom_scattermore(pixels=c(1000,1000), pointsize=8) +
        theme_pubr() +
        color_palette(palette=get_palette("jco", 24)) +
        theme(aspect.ratio=1)
    
    ## check silhouettes
    X = embedding %>%
        left_join(drug_targets %>% distinct(DRUG_ID,TARGET_PATHWAY), 
                  by=c("index"="DRUG_ID")) %>%
        # filter(!(TARGET_PATHWAY %in% c("Other", "Other, kinases", "Unclassified"))) %>%
        mutate(lab = as.numeric(as.factor(TARGET_PATHWAY)))
        
    
    dis = X %>%
        column_to_rownames("index") %>%
        dplyr::select(starts_with("UMAP")) %>%
        dist() %>%
        as.matrix()
    
    sil = silhouette(X %>% pull(lab), dis^2)
    
    ## Are silhouettes good only when there is a significant association 
    ## with event(s) in the gene target(s)?
    pathways_w_events = drug_targets %>% 
        mutate(is_found = as.factor(TARGET %in% (found_targets %>% 
                                     filter(lr_padj<THRESH_FDR) %>% 
                                     pull(GENE)))) %>%
        count(TARGET_PATHWAY, is_found, .drop=FALSE) %>%
        filter(is_found=="TRUE") %>%
        mutate(label = sprintf("%s (%s)", TARGET_PATHWAY, n))
    
    plts[['associations-target_pathway-silhouettes']] = X %>%
        left_join(as.data.frame(sil[,]), by=c("lab"="cluster")) %>% 
        left_join(pathways_w_events, by=c("TARGET_PATHWAY")) %>%
        arrange(TARGET_PATHWAY) %>%
        dplyr::select(sil_width, label, TARGET_PATHWAY) %>% 
        ggboxplot(x="label", y="sil_width", fill="TARGET_PATHWAY", 
                  outlier.size=0.1, palette=get_palette("jco", 24)) + 
        guides(fill="none") +
        labs(x='Target Pathway (N. Targets)', y='Silhouette Score') +
        coord_flip()
    
    
    # how many drugs of each pathway fall in each leiden cluster?
    plts[['associations-leiden-vs_target_pathway']] = embedding %>%
        left_join(drug_targets %>% distinct(DRUG_ID,TARGET_PATHWAY), 
                  by=c("index"="DRUG_ID")) %>%
        group_by(leiden_labels,TARGET_PATHWAY) %>%
        summarize(n = n()) %>%
        mutate(prop = n / sum(n),
               leiden_labels = as.character(leiden_labels)) %>%
        ggballoonplot(x="TARGET_PATHWAY", y="leiden_labels",
                      size="prop", fill="aquamarine3") +
        scale_size_continuous(range=c(0.5,2.5))

    # are there other biological reasons for the umap clusters?
    plts[['associations-leiden-umap']] = embedding %>%
        left_join(drug_targets %>% distinct(DRUG_ID,TARGET_PATHWAY), 
                  by=c("index"="DRUG_ID")) %>%
        # filter(!(TARGET_PATHWAY %in% c("Other", "Other, kinases", "Unclassified"))) %>%
        mutate(leiden_labels=as.factor(leiden_labels)) %>%
        ggplot(aes(x=UMAP0, y=UMAP1, color=leiden_labels)) +
        geom_scattermore(pixels=c(1000,1000), pointsize=8) +
        theme_pubr() +
        theme(aspect.ratio=1)
    
    gene_sets = embedding %>%
        left_join(
            drug_targets %>%
            distinct(DRUG_ID,TARGET),
          by=c("index"="DRUG_ID")) %>%
        distinct(leiden_labels,TARGET) %>%
        get_sets("leiden_labels","TARGET")
    
    universe = unique(unlist(gene_sets))
    results = sapply(names(gene_sets), function(gene_set){
        enrichment = run_enrichment(gene_sets[[gene_set]], universe, ontologies)
        return(enrichment)
    }, simplify=FALSE)
    results = get_enrichment_result(results)
    
    plts_enrichment = sapply(names(results), function(onto){
        result = results[[onto]]
        res = new("compareClusterResult", compareClusterResult = result)
        plt = dotplot(res) + labs(title=onto, x='Leiden Cluster')
        return(plt)
    }, simplify=FALSE)
    names(plts_enrichment) = sprintf('associations-leiden-enrichment-%s',
                                     names(plts_enrichment))
    plts = c(plts,plts_enrichment)
    
    # overall best drug-exon rankings per drug
    X = rankings %>% 
        group_by(DRUG_NAME) %>% 
        slice_min(combined_ranking, n=1) %>% 
        ungroup() %>% 
        arrange(index) %>%
        left_join(drug_targets %>% 
                      distinct(DRUG_ID,TARGET) %>% 
                      mutate(is_target=TRUE), 
                  by=c("DRUG_ID", "GENE"="TARGET")) %>%
        mutate(is_target=!is.na(is_target))
    
    ## all drugs    
    plts[['associations-rankings-all']] = X %>%
        slice_min(index, n=top_n) %>%
        ggbarplot(x="DRUG_NAME", y="combined_ranking", 
                  fill="is_target", color="drug_screen") +
        color_palette(c("white","black")) + 
        fill_palette("Set2") + 
        geom_text(aes(label=index), X %>% slice_min(index, n=top_n), 
                  size=1, family='Arial') +
        geom_text(aes(y=0.005, label=event_gene), 
                  X %>% slice_min(index, n=top_n), size=1, family='Arial') +
        labs(x="Drug", y="Ranking Ratio Sum", 
             fill="Is Drug Target", color="Drug Screen") +
        coord_flip()
        
    ## without targets
    plts[['associations-rankings-no_target']] = X %>%
        filter(DRUG_ID %in% drugs_oi) %>%
        ggbarplot(x="DRUG_NAME", y="combined_ranking", 
                  fill="is_target", color="drug_screen") +
        color_palette(c("white","black")) + 
        fill_palette("Set2") + 
        geom_text(aes(label=index), X %>% filter(DRUG_ID %in% drugs_oi), 
                  size=1, family='Arial') +
        geom_text(aes(y=0.05, label=event_gene), 
                  X %>% filter(DRUG_ID %in% drugs_oi), size=1, family='Arial') +
        labs(x="Drug", y="Ranking Ratio Sum", 
             fill="Is Drug Target", color="Drug Screen") +
        coord_flip()
    
    return(plts)
}


plot_drug_rec = function(estimated_response, drug_screen, drug_targets){
    X = estimated_response %>%
        left_join(drug_screen %>% mutate(real_ic50 = log(IC50_PUBLISHED)), 
                  by=c("drug_screen"="DATASET","DRUG_ID", "sample"="ARXSPAN_ID")) %>%
        drop_na(real_ic50, predicted_ic50)
    
    corrs = X %>% 
        group_by(drug_screen, sample) %>%
        summarize(correlation = cor(real_ic50, predicted_ic50, method="spearman"))
    
    plts = list()
    plts[["drug_rec-spearmans"]] = corrs %>% 
        ggviolin(x="drug_screen", y="correlation", trim=TRUE,
                 color=NA, fill="drug_screen", palette="jco") +
        geom_boxplot(width=0.1, outlier.size=0.1) +
        geom_hline(yintercept=0, linetype="dashed") +
        guides(fill="none") + 
        labs(x="Drug Screen", y="Spearman Correlation") +
        geom_text(aes(y=1, label=n), 
                  corrs %>% count(drug_screen), 
                  size=1, family='Arial')
    
    # best and worse correlations
    samples_oi = corrs %>% 
        group_by(drug_screen) %>%
        arrange(correlation) %>% 
        filter(row_number()==1 | row_number()==n()) %>%
        left_join(X, by=c("drug_screen","sample"))
    
    plts[['drug_rec-best_worse']] = samples_oi %>% 
        ggscatter(x="real_ic50", y="predicted_ic50", size=0.5) + 
        facet_wrap(drug_screen~sample, scales='free') +
        stat_cor(method="spearman", size=2) + 
        geom_abline(intercept=0, slope=1, linetype="dashed") +
        labs(x="Real log(IC50)", y="Predicted log(IC50)") + 
        theme_pubr(border=TRUE) +
        theme(strip.text.x = element_text(size=6, family='Arial'))
    
    # influence of pathway in ranking capacity?
    corrs_bypath = X %>%
        left_join(drug_targets %>% distinct(DRUG_NAME, TARGET_PATHWAY), 
                  by="DRUG_NAME") %>% 
        group_by(drug_screen, sample, TARGET_PATHWAY) %>%
        summarize(correlation = cor(real_ic50, predicted_ic50, method="spearman"))
    
    plts[["drug_rec-spearmans_by_pathway"]] = corrs_bypath %>%
        drop_na() %>%
        # filter(!(TARGET_PATHWAY %in% c("Other", "Other, kinases", "Unclassified"))) %>%
        arrange(TARGET_PATHWAY) %>%
        ggboxplot(x="TARGET_PATHWAY", y="correlation", outlier.size=0.1,
                  fill="TARGET_PATHWAY", palette=get_palette("jco", 24)) + 
        ylim(-1,1) + 
        facet_wrap(~drug_screen) +
        labs(x="Target Pathway", y="Spearman Correlation") +
        guides(fill="none") + 
        coord_flip() +
        theme(strip.text.x = element_text(size=6, family='Arial'))
    
    
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
        coord_flip() +
        theme(strip.text.x = element_text(size=6, family='Arial'))
    
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
        coord_flip() +
        theme(strip.text.x = element_text(size=6, family='Arial'))
    
    return(plts)
}


plot_drugs_common_targets = function(models, drug_targets){
    require(proxy)
    X = drug_targets %>% distinct(DRUG_NAME,TARGET) %>% mutate(is_target=1) %>% pivot_wider(id_cols=TARGET, names_from=DRUG_NAME, values_from=is_target, values_fill=FALSE)
    sim = simil(X[,-1], method='Jaccard', by_rows=FALSE) %>% as.matrix()
    
    ind = which(upper.tri(sim, diag=FALSE), arr.ind=TRUE)
    sim = data.frame(
        col = dimnames(sim)[[2]][ind[,2]],
        row = dimnames(sim)[[1]][ind[,1]],
        jaccard = sim[ind]
    )
    sim %>% filter(jaccard==1)
    
    drugs_oi = c('AZD6482','TGX-221') # actually TGX221
    models %>% 
        mutate(drug_name=gsub("(.*),.*", "\\1",drug_name)) %>% 
        filter(drug_name %in% drugs_oi) %>% #& lr_padj<THRESH_FDR) %>% 
        pivot_wider(id_cols=event_gene, names_from=drug_name, values_from=spldep_coefficient) %>%
        ggscatter(x='TGX-221', y='AZD6482') +
        geom_abline(slope=1, intercept=0, linetype='dashed') +
        stat_cor()
        count(drug_name)
}


plot_drug_clusters = function(clusters, metadata, drug_screen, estimated_response, spldep_ccle){
    X = drug_screen %>%
        distinct(DRUG_NAME, DRUG_ID, DATASET) %>%
        filter(str_detect(DRUG_NAME,"THZ-2-49")) %>%
        left_join(clusters, by=c("DRUG_ID", "DATASET"="drug_screen")) %>%
        left_join(metadata, by=c("index"="DepMap_ID")) %>%
        left_join(drug_screen, by=c("index"="ARXSPAN_ID", "DRUG_NAME", 
                                    "DRUG_ID", "DATASET")) %>%
        left_join(estimated_response, by=c("DRUG_ID", "index"="sample",
                                           "DATASET"="drug_screen")) %>%
        mutate(log_ic50 = log(IC50_PUBLISHED),
               in_training_set = !is.na(log_ic50),
               leiden_labels = as.factor(leiden_labels)) %>%
        left_join(spldep_ccle %>% 
                      column_to_rownames("index") %>% 
                      t() %>% as.data.frame() %>% 
                      rownames_to_column("index"), by="index")

    X %>% ggscatter(x="UMAP0", y="UMAP1", color="log_ic50") + scale_color_gradient2(low="blue",mid="white",high="red") + facet_wrap(~DRUG_NAME) 
    X %>% ggscatter(x="UMAP0", y="UMAP1", color="in_training_set") + facet_wrap(~DRUG_NAME)
    X %>% drop_na(HsaEX6007325) %>% ggscatter(x="UMAP0", y="UMAP1", color="HsaEX6007325") + scale_color_gradient2(low="blue", mid="white", high="red") + facet_wrap(~DRUG_NAME)
    X %>% ggscatter(x="UMAP0", y="UMAP1", color="leiden_labels") + facet_wrap(~DRUG_NAME)
     X %>% ggscatter(x="UMAP0", y="UMAP1", color="primary_disease", palette=get_palette("Paired", length(unique(X[["primary_disease"]])))) + facet_wrap(~DRUG_NAME)

}


make_plots = function(models, drug_targets, embedding, estimated_response, drug_screen, ontologies, rankings){
    plts = list(
        plot_associations(models, drug_targets, embedding, ontologies, rankings),
        plot_drug_rec(estimated_response, drug_screen, drug_targets)
    )
    plts = do.call(c,plts)
    return(plts)
}


save_plt = function(plts, plt_name, extension='.pdf', 
                    directory='', dpi=350, format=TRUE,
                    width = par("din")[1], height = par("din")[2]){
    plt = plts[[plt_name]]
    if (format){
        plt = ggpar(plt, font.title=10, font.subtitle=10, font.caption=10, 
                    font.x=8, font.y=8, font.legend=8,
                    font.tickslab=6, font.family='Arial')    
    }
    filename = file.path(directory,paste0(plt_name,extension))
    save_plot(filename, plt, base_width=width, base_height=height, dpi=dpi, units='cm')
}


save_plots = function(plts, figs_dir){
    # drug-event associations
    save_plt(plts, 'associations-lr_pvalues', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'associations-lr_fdr', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'associations-drug_counts', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'associations-top_drug_counts', '.pdf', figs_dir, width=8, height=8)
    save_plt(plts, 'associations-spldep_counts', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'associations-top_spldep_counts', '.pdf', figs_dir, width=8, height=8)
    save_plt(plts, 'associations-top_found_targets', '.pdf', figs_dir, width=25, height=19)
    save_plt(plts, 'associations-top_no_targets', '.pdf', figs_dir, width=25, height=25)
    save_plt(plts, 'associations-target_pathway-counts', '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, 'associations-target_pathway-umap', '.pdf', figs_dir, width=10, height=10)
    save_plt(plts, 'associations-target_pathway-silhouettes', '.pdf', figs_dir, width=6.5, height=6)
    save_plt(plts, 'associations-leiden-vs_target_pathway', '.pdf', figs_dir, width=10.5, height=5)
    save_plt(plts, 'associations-leiden-umap', '.pdf', figs_dir, width=10, height=10)
    save_plt(plts, 'associations-leiden-enrichment-hallmarks', '.pdf', figs_dir, width=8, height=6)
    save_plt(plts, 'associations-leiden-enrichment-GO_BP', '.pdf', figs_dir, width=15, height=14)
    
    save_plt(plts, 'associations-target_ranking-vio', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'associations-target_ranking-cdf', '.pdf', figs_dir, width=5, height=5)
    
    save_plt(plts, 'associations-rankings-all', '.pdf', figs_dir, width=8, height=8)
    save_plt(plts, 'associations-rankings-no_target', '.pdf', figs_dir, width=8, height=8)
    
    # drug recommendations
    save_plt(plts, 'drug_rec-spearmans', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'drug_rec-best_worse', '.pdf', figs_dir, width=8, height=8)
    save_plt(plts, 'drug_rec-spearmans_by_pathway', '.pdf', figs_dir, width=8, height=6)
    
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
    embedding_file = args$embedding_file
    estimated_response_file = args$estimated_response_file
    drug_screen_file = args$drug_screen_file
    figs_dir = args$figs_dir
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    metadata = read_tsv(metadata_file)
    models = read_tsv(models_file) %>% 
        mutate(event_gene = paste0(EVENT,'_',GENE),
               event_type = gsub('Hsa','',gsub("[^a-zA-Z]", "",EVENT)))
    drug_targets = read_csv(drug_targets_file) %>% 
        mutate(DRUG_NAME=toupper(DRUG_NAME)) %>%
        dplyr::select(DRUG_ID,DRUG_NAME,TARGET,TARGET_PATHWAY) %>%
        separate_rows(TARGET) %>%
        distinct()
    embedding = read_tsv(embedding_file)
    estimated_response = read_tsv(estimated_response_file)
    drug_screen = read_csv(drug_screen_file)
    ontologies = list(
        "hallmarks" = read.gmt(file.path(msigdb_dir,'h.all.v7.4.symbols.gmt')),
        "oncogenic_signatures" = read.gmt(file.path(msigdb_dir,'c6.all.v7.4.symbols.gmt')),
        "GO_BP" = read.gmt(file.path(msigdb_dir,'c5.go.bp.v7.4.symbols.gmt'))
    )
    clusters = read_tsv(clusters_file)
    spldep_ccle = read_tsv(spldep_ccle_file) %>%
        filter(index %in% unique(models[["EVENT"]]))
    
    rankings = models %>% 
        filter(lr_padj<0.1) %>% 
        group_by(DRUG_ID) %>% 
        arrange(DRUG_ID, -spldep_coefficient) %>% 
        mutate(ranking = row_number(), 
               ranking_ratio = ranking/n(), 
               ranking_drug = ranking_ratio) %>% 
        ungroup() %>% 
        group_by(event_gene) %>% 
        arrange(event_gene, -spldep_coefficient) %>% 
        mutate(ranking = row_number(), 
               ranking_ratio = ranking/n(), 
               ranking_event = ranking_ratio) %>%
        ungroup() %>%
        dplyr::select(event_gene, EVENT, GENE, DRUG_ID, spldep_coefficient,
                      drug_screen, ranking_drug, ranking_event) %>%
        # shortest distance to diagonal 
        mutate(combined_ranking = ranking_event + ranking_drug) %>% 
        arrange(combined_ranking) %>%
        mutate(index = row_number()) %>%
        left_join(drug_targets %>% distinct(DRUG_ID,DRUG_NAME), by="DRUG_ID")
    
    # make plots
    plts = make_plots(models, drug_targets, embedding, estimated_response, drug_screen, ontologies, rankings)
    
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