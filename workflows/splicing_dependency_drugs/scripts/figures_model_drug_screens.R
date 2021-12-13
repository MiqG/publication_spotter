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
require(extrafont)
require(cluster)
require(ComplexHeatmap)
require(ggplotify)
require(gridExtra)

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
# models_file = file.path(RESULTS_DIR,'files','model_summaries_drug_response-EX.tsv.gz')
# drug_targets_file = file.path(RAW_DIR,'GDSC','screened_compunds_rel_8.2.csv')
# figs_dir = file.path(RESULTS_DIR,'figures','model_drug_screens')
# embedding_file = file.path(RESULTS_DIR,'files','embedded_drug_associations-EX.tsv.gz')
# estimated_response_file = file.path(RESULTS_DIR, 'files', 'estimated_drug_response-GDSC1-EX', 'estimated_drug_response_by_drug.tsv.gz')
# drug_screen_file = file.path(PREP_DIR,'drug_screens','GDSC1.tsv.gz')

##### FUNCTIONS #####
get_sets = function(df, set_names, set_values){
    sets = df[,c(set_values,set_names)] %>%
        distinct() %>%
        with(., split(get(set_values),get(set_names)))
    sets = sapply(sets, unique, simplufy=FALSE)
    return(sets)
}

plot_associations = function(models, drug_targets, embedding){
    top_n = 25
    X = models %>% 
        left_join(drug_targets[,c("DRUG_ID","DRUG_NAME")], by="DRUG_ID")
    
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
        count(DRUG_NAME) %>% 
        arrange(n) %>% 
        ggbarplot(x='DRUG_NAME', y='n', fill='#2a9d8f', color=NA) + 
        labs(x='Drug', y='N. Significant Associations')
    
    plts[['associations-top_drug_counts']] = X %>% 
        filter(lr_padj < THRESH_FDR) %>% 
        count(DRUG_NAME) %>% 
        arrange(n) %>% 
        tail(top_n) %>%
        ggbarplot(x='DRUG_NAME', y='n', fill='#2a9d8f', color=NA) + 
        labs(x='Drug', y='N. Significant Associations', 
             title=sprintf('Top %s', top_n)) +
        coord_flip()
    
    plts[['associations-spldep_counts']] = X %>% 
        filter(lr_padj < THRESH_FDR) %>% 
        count(event_gene) %>% 
        arrange(n) %>% 
        ggbarplot(x='event_gene', y='n', fill='orange', color=NA) + 
        labs(x='Event & Gene', y='N. Significant Associations')
    
    plts[['associations-top_spldep_counts']] = X %>% 
        filter(lr_padj < THRESH_FDR) %>% 
        count(event_gene) %>% 
        arrange(n) %>% 
        tail(top_n) %>%
        ggbarplot(x='event_gene', y='n', fill='orange', color=NA) + 
        labs(x='Event & Gene', y='N. Significant Associations', 
             title=sprintf('Top %s', top_n)) +
        coord_flip()
    
    # is there any drug that is significantly associated to the splicing
    # dependency of an event of a gene that it targets?
    found_targets = X %>% 
        filter(lr_padj < THRESH_FDR) %>% 
        left_join(drug_targets, by=c("DRUG_ID","DRUG_NAME","GENE"="TARGET")) %>%
        drop_na(TARGET_PATHWAY) %>%
        mutate(is_target=TRUE)
    drugs_oi = found_targets %>% pull(DRUG_NAME) %>% unique()
    events_oi = found_targets %>% pull(EVENT) %>% unique()
    
    plts[['associations-top_found_targets']] = X %>% 
        filter(lr_padj<THRESH_FDR & DRUG_NAME%in%drugs_oi) %>% 
        group_by(DRUG_NAME) %>% 
        slice_max(order_by=abs(spldep_coefficient), n=15) %>% 
        left_join(distinct(found_targets[,c('DRUG_NAME','GENE','is_target')]), 
                  by=c('DRUG_NAME','GENE')) %>% 
        left_join(distinct(found_targets[,c('DRUG_NAME','TARGET_PATHWAY')]), 
                  by=c('DRUG_NAME')) %>% 
        mutate(is_target=replace_na(is_target,FALSE), 
               log10_pvalue=-log10(lr_pvalue), 
               drug_name_clean=sprintf('%s | %s',DRUG_NAME,TARGET_PATHWAY), 
               event_gene=reorder_within(event_gene, spldep_coefficient, DRUG_NAME)) %>% 
        ggbarplot(x='event_gene', y='spldep_coefficient', fill='is_target', 
                  palette='lancet', color=NA) + 
        facet_wrap(~drug_name_clean, scales='free') + 
        scale_x_reordered() +     
        labs(x='Event & Gene', y='Effect Size', fill='Is Drug Target', 
             title='Top 15 Significant Effect Sizes') + 
        coord_flip() +
        theme_pubr(border=TRUE)
    
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
    
    # are association profiles informative of drug mechanism of action?
    plts[['associations-target_pathway-counts']] = drug_targets %>% 
        distinct(DRUG_ID,TARGET_PATHWAY) %>%
        count(TARGET_PATHWAY) %>%
        filter(!(TARGET_PATHWAY %in% c("Other", "Other, kinases", "Unclassified"))) %>%
        ggbarplot(x='TARGET_PATHWAY', y="n",
                  fill="TARGET_PATHWAY", color=NA, 
                  palette=get_palette("jco", 21)) +
        labs(x='Target Pathway', y='Count') +
        coord_flip() +
        theme_pubr(legend="none")
    
    plts[['associations-target_pathway-umap']] = embedding %>%
        left_join(drug_targets %>% distinct(DRUG_ID,TARGET_PATHWAY), 
                  by=c("index"="DRUG_ID")) %>%
        filter(!(TARGET_PATHWAY %in% c("Other", "Other, kinases", "Unclassified"))) %>%
        ggscatter(x="UMAP0", y="UMAP1", color="TARGET_PATHWAY", 
                  palette=get_palette("jco", 21)) +
        theme(aspect.ratio=1)
    
    ## check silhouettes
    X = embedding %>%
        left_join(drug_targets %>% distinct(DRUG_ID,TARGET_PATHWAY), 
                  by=c("index"="DRUG_ID")) %>%
        filter(!(TARGET_PATHWAY %in% c("Other", "Other, kinases", "Unclassified"))) %>%
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
        mutate(is_found = as.factor(TARGET %in% (found_targets %>% pull(GENE)))) %>%
        count(TARGET_PATHWAY, is_found, .drop=FALSE) %>%
        filter(is_found=="TRUE") %>%
        mutate(label = sprintf("%s (%s)", TARGET_PATHWAY, n))
    
    plts[['associations-target_pathway-silhouettes']] = X %>%
        left_join(as.data.frame(sil[,]), by=c("lab"="cluster")) %>% 
        left_join(pathways_w_events, by=c("TARGET_PATHWAY")) %>%
        arrange(TARGET_PATHWAY) %>%
        dplyr::select(sil_width, label, TARGET_PATHWAY) %>% 
        ggboxplot(x="label", y="sil_width", fill="TARGET_PATHWAY", 
                  outlier.size=0.5, palette=get_palette("jco", 21)) + 
        guides(fill="none") +
        labs(x='Target Pathway (N. Targets)', y='Silhouette Score') +
        coord_flip()
    
    # in some cases, there seems to be an association between dependency profiles
    # and drug mechanism of action / target pathway, but in others not.
    # Does this have something to do with the number of specific splicing events
    # whose dependency is associated with drugs of a certain target pathway?
#     X = models %>% 
#         left_join(drug_targets, by="DRUG_ID")
    
#     events_oi = X %>% 
#         filter(lr_padj < THRESH_FDR) %>% 
#         distinct(EVENT,TARGET_PATHWAY) %>% 
#         filter(!(TARGET_PATHWAY %in% c("Other", "Other, kinases", "Unclassified"))) %>% 
#         get_sets('TARGET_PATHWAY', 'EVENT')
#     m = events_oi %>% list_to_matrix() %>% make_comb_mat()
#     UpSet(m, comb_order = order(comb_size(m)))
#     as.ggplot(grid.grabExpr(draw()))
    
    
    return(plts)
}


plot_drug_rec = function(estimated_response, drug_screen, drug_targets){
    X = estimated_response %>%
        left_join(drug_screen %>% mutate(real_ic50 = log(IC50_PUBLISHED)), 
                  by=c("DRUG_ID", "sample"="ARXSPAN_ID")) %>%
        drop_na(real_ic50, predicted_ic50)
    
    corrs = X %>% 
        group_by(sample) %>%
        summarize(correlation = cor(real_ic50, predicted_ic50, method="spearman"))
    
    plts = list()
    plts[["drug_rec-spearmans"]] = corrs %>% 
        gghistogram(x="correlation", color=NA, fill="gold4") + 
        xlim(-1,1) + 
        geom_vline(xintercept=median(corrs[['correlation']]), # 0.64
                   linetype="dashed") + 
        labs(x="Spearman Correlation", y="Count")
    
    # best and worse correlations
    samples_oi = corrs %>% 
        arrange(correlation) %>% 
        filter(row_number()==1 | row_number()==n()) %>% 
        pull(sample)
    plts[['drug_rec-best_worse']] = X %>% 
        filter(sample %in% samples_oi) %>% 
        ggscatter(x="real_ic50", y="predicted_ic50", size=1) + 
        facet_wrap(~sample, scales='free') +
        stat_cor(method="spearman") + 
        geom_abline(intercept=0, slope=1, linetype="dashed") +
        labs(x="Real log(IC50)", y="Predicted log(IC50)") + 
        theme_pubr(border=TRUE)
    
    # influence of pathway in ranking capacity?
    corrs_bypath = X %>%
        left_join(drug_targets %>% distinct(DRUG_NAME, TARGET_PATHWAY), 
                  by="DRUG_NAME") %>% 
        group_by(sample, TARGET_PATHWAY) %>%
        summarize(correlation = cor(real_ic50, predicted_ic50, method="spearman"))
    
    plts[["drug_rec-spearmans_by_pathway"]] = corrs_bypath %>%
        drop_na() %>%
        filter(!(TARGET_PATHWAY %in% c("Other", "Other, kinases", "Unclassified"))) %>%
        arrange(TARGET_PATHWAY) %>%
        ggboxplot(x="TARGET_PATHWAY", y="correlation", outlier.size=0.5,
                  fill="TARGET_PATHWAY", palette=get_palette("jco", 21)) + 
        ylim(-1,1) + 
        labs(x="Target Pathway", y="Spearman Correlation") +
        guides(fill="none") + 
        coord_flip()
    
    
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

make_plots = function(models, drug_targets, embedding, estimated_response, drug_screen){
    plts = list(
        plot_associations(models, drug_targets, embedding),
        plot_drug_rec(estimated_response, drug_screen, drug_targets)
    )
    plts = do.call(c,plts)
    return(plts)
}


save_plt = function(plts, plt_name, extension='.pdf', 
                    directory='', dpi=350, 
                    width = par("din")[1], height = par("din")[2]){
    plt = plts[[plt_name]]
    plt = ggpar(plt, font.title=11, font.subtitle=10, font.caption=10, 
                font.x=10, font.y=10, font.legend=10,
                font.tickslab=8, font.family='Arial')
    filename = file.path(directory,paste0(plt_name,extension))
    save_plot(filename, 
              plt, 
              base_width=width, base_height=height, dpi=dpi)
}


save_plots = function(plts, figs_dir){
    # drug-event associations
    save_plt(plts, 'associations-lr_pvalues', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'associations-lr_fdr', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'associations-drug_counts', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'associations-top_drug_counts', '.pdf', figs_dir, width=8, height=5)
    save_plt(plts, 'associations-spldep_counts', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'associations-top_spldep_counts', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'associations-top_found_targets', '.pdf', figs_dir, width=20, height=20)
    save_plt(plts, 'associations-target_pathway-counts', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'associations-target_pathway-umap', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'associations-target_pathway-silhouettes', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'associations-target_ranking-vio', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'associations-target_ranking-cdf', '.pdf', figs_dir, width=5, height=5)
    
    # drug recommendations
    save_plt(plts, 'drug_rec-spearmans', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'drug_rec-best_worse', '.pdf', figs_dir, width=5, height=2.7)
    save_plt(plts, 'drug_rec-spearmans_by_pathway', '.pdf', figs_dir, width=5, height=5)
    
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
    drug_screen = read_tsv(drug_screen_file)
    
    # make plots
    plts = make_plots(models, drug_targets, embedding)

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