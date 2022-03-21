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
# EDA associations
# - pvalues
# - n drug/events significantly associated
# - association agreement between GDSC1 and GDSC2
# - overlaps with cell lines used to model splicing dependency
# - overlaps with cell lines used in GDSC1 and GDSC2
# 
# Rankings
# - significant associations with drug target(s)
# - significant associations distance from drug target(s)
# - significant associations better than with drug target(s)
# - drugs that share target(s) ? 
# - significant associations with drug whose target is not significantly associated
#   at similar overall ranking levels
#
# Association profiles
# - splicing dependency coefficients separate drugs in clusters that separate 
#   their targets (splicing dependencies that differentiate the clusters?)

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
require(writexl)

ROOT = here::here()
source(file.path(ROOT,"src","R","utils.R"))

# variables
THRESH_FDR = 0.1
THRESH_PVALUE = 0.05
THRESH_NOBS = 20
RANDOM_SEED = 1234

# Development
# -----------
# RAW_DIR = file.path(ROOT,"data","raw")
# PREP_DIR = file.path(ROOT,"data","prep")
# RESULTS_DIR = file.path(ROOT,"results","splicing_dependency_drugs")
# MODELS_DIR = file.path(ROOT,"results","model_splicing_dependency")
# models_file = file.path(RESULTS_DIR,"files","model_summaries_drug_response-EX.tsv.gz")
# drug_targets_file = file.path(PREP_DIR,'drug_screens','drug_targets.tsv.gz')
# figs_dir = file.path(RESULTS_DIR,"figures","model_drug_screens")
# embedding_file = file.path(RESULTS_DIR,"files","embedded_drug_associations-EX.tsv.gz")
# estimated_response_file = file.path(RESULTS_DIR,"files","estimated_drug_response_by_drug-EX.tsv.gz")
# drug_screens_dir = file.path(PREP_DIR,'drug_screens')
# msigdb_dir = file.path(ROOT,"data","raw","MSigDB","msigdb_v7.4","msigdb_v7.4_files_to_download_locally","msigdb_v7.4_GMTs")
# clusters_file = file.path(RESULTS_DIR,"files","cluster_estimated_drug_response-merged-EX.tsv.gz")
# metadata_file = file.path(PREP_DIR,"metadata","CCLE.tsv.gz")
# spldep_models_file = file.path(MODELS_DIR,"files","models_gene_dependency-EX","model_summaries.tsv.gz")
# spldep_ccle_file = file.path(MODELS_DIR,"files","splicing_dependency-EX","mean.tsv.gz")
# paths_real_file = file.path(RESULTS_DIR,"files","ppi","shortest_path_lengths_to_drug_targets-EX.tsv.gz")
# paths_random_file = file.path(RESULTS_DIR,"files","ppi","shortest_path_lengths_to_drug_targets-random.tsv.gz")
# rnai_file = file.path(PREP_DIR,"demeter2","CCLE.tsv.gz")


##### FUNCTIONS #####
load_drug_screens = function(drug_screens_dir){
    filenames = cbind(
        drug_screens_dir,
        expand_grid(c("train","test"),c("GDSC1.tsv.gz","GDSC2.tsv.gz")
                   )
    )
    dfs = apply(filenames, 1, function(x){
        x = as.vector(x)
        f = file.path(x[1],x[2],x[3])
        df = read_tsv(f)
        return(df)
    })
    dfs = do.call(rbind,dfs)
    return(dfs)
}

compute_rankings = function(models, drug_targets){
    rankings = models %>% 
        filter(lr_padj<THRESH_FDR  & n_obs>THRESH_NOBS) %>% 
        group_by(ID) %>% 
        arrange(ID, -spldep_coefficient) %>% 
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
        dplyr::select(event_gene, EVENT, GENE, ID, DRUG_ID, spldep_coefficient,
                      drug_screen, ranking_drug, ranking_event, n_obs) %>%
        mutate(combined_ranking = ranking_event + ranking_drug) %>% 
        arrange(combined_ranking) %>%
        mutate(index = row_number()) %>%
        left_join(drug_targets %>% distinct(DRUG_ID,DRUG_NAME), by="DRUG_ID") %>%
        mutate(DRUG_NAME_CLEAN = ID) %>%
        separate(DRUG_NAME_CLEAN, c("drug_id","min_conc","max_conc"), sep="_") %>%
        mutate(DRUG_NAME_CLEAN = sprintf("%s (%s-%s)", DRUG_NAME, min_conc, max_conc))
    
    return(rankings)
}


load_ontologies = function(msigdb_dir){
    ontologies = list(
        "hallmarks" = read.gmt(file.path(msigdb_dir,"h.all.v7.4.symbols.gmt")),
        "oncogenic_signatures" = read.gmt(file.path(msigdb_dir,"c6.all.v7.4.symbols.gmt")),
        "GO_BP" = read.gmt(file.path(msigdb_dir,"c5.go.bp.v7.4.symbols.gmt"))
    )
    return(ontologies)
}

get_sets = function(df, set_names, set_values){
    sets = df[,c(set_values,set_names)] %>%
        distinct() %>%
        with(., split(get(set_values),get(set_names)))
    sets = sapply(sets, unique, simplufy=FALSE)
    return(sets)
}


run_enrichment = function(genes, universe, ontologies){
    enrichments = list()
    enrichments[["hallmarks"]] = enricher(genes, TERM2GENE=ontologies[["hallmarks"]], universe=universe)
    enrichments[["GO_BP"]] = enricher(genes, TERM2GENE=ontologies[["GO_BP"]], universe=universe)
    enrichments[["oncogenic_signatures"]] = enricher(genes, TERM2GENE=ontologies[["oncogenic_signatures"]], universe=universe)
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


plot_associations = function(models, spldep_ccle, drug_screen){
    top_n = 25
    X = models %>% 
        left_join(drug_screen %>% distinct(DRUG_ID,DRUG_NAME) %>% drop_na(), by="DRUG_ID") %>%
        mutate(DRUG_NAME_CLEAN = ID) %>%
        separate(DRUG_NAME_CLEAN, c("drug_id","min_conc","max_conc"), sep="_") %>%
        mutate(DRUG_NAME_CLEAN = sprintf("%s (%s-%s)", DRUG_NAME, min_conc, max_conc))
    
    plts = list()
    
    # what are the distributions of p-values and FDR?
    plts[["associations-lr_pvalues"]] = X %>% 
        gghistogram(x="lr_pvalue", fill="lightblue", color=NA, bins=100) + 
        geom_vline(xintercept=median(X[["lr_pvalue"]], na.rm=TRUE), 
                   linetype="dashed") +
        labs(x="LR Test p-value")
    
    plts[["associations-lr_fdr"]] = X %>% 
        gghistogram(x="lr_padj", fill="darkred", color=NA, bins=100) + 
        geom_vline(xintercept=median(X[["lr_padj"]], na.rm=TRUE), 
                   linetype="dashed") +
        labs(x="LR Test FDR")
    
    # how do number of observations relate to association coefficient?
    plts[["associations-nobs_vs_coef"]] = X %>% 
        filter(lr_padj < THRESH_FDR) %>% 
        ggplot(aes(x=n_obs, y=spldep_coefficient)) + 
        geom_scattermore(pixels=c(1000,1000), pointsize=2, alpha=0.5) + 
        theme_pubr() + 
        geom_vline(xintercept=THRESH_NOBS, linetype="dashed") + 
        labs(x="N. Observations", y="Assoc. Coefficient")
    
    # how many drugs and events are significantly associated?
    plts[["associations-drug_counts"]] = X %>% 
        filter(lr_padj < THRESH_FDR & n_obs > THRESH_NOBS) %>% 
        count(drug_screen, DRUG_NAME_CLEAN) %>% 
        group_by(DRUG_NAME_CLEAN) %>%
        mutate(total=sum(n)) %>%
        arrange(total) %>%
        ggbarplot(x="DRUG_NAME_CLEAN", y="n", fill="drug_screen", palette="jco", color=NA) + 
        labs(x="Drug", y="N. Significant Associations")
    
    plts[["associations-top_drug_counts"]] = X %>% 
        filter(lr_padj < THRESH_FDR & n_obs > THRESH_NOBS) %>% 
        count(drug_screen, DRUG_NAME_CLEAN) %>% 
        group_by(DRUG_NAME_CLEAN) %>%
        mutate(total=sum(n)) %>%
        arrange(-total) %>%
        ungroup() %>%
        mutate(index=as.numeric(factor(DRUG_NAME_CLEAN, levels=unique(DRUG_NAME_CLEAN)))) %>%
        filter(index <= top_n) %>%
        ggbarplot(x="DRUG_NAME_CLEAN", y="n", fill="drug_screen", palette="jco", color=NA) + 
        labs(x="Drug", y="N. Significant Associations", 
             title=sprintf("Top %s", top_n)) +
        coord_flip()
    
    plts[["associations-spldep_counts"]] = X %>% 
        filter(lr_padj < THRESH_FDR & n_obs > THRESH_NOBS) %>% 
        count(drug_screen, event_gene) %>% 
        group_by(event_gene) %>%
        mutate(total=sum(n)) %>%
        arrange(total) %>%
        ggbarplot(x="event_gene", y="n", fill="drug_screen", palette="jco", color=NA) + 
        labs(x="Event & Gene", y="N. Significant Associations")
    
    plts[["associations-top_spldep_counts"]] = X %>% 
        filter(lr_padj < THRESH_FDR & n_obs > THRESH_NOBS) %>% 
        count(drug_screen, event_gene) %>% 
        group_by(event_gene) %>%
        mutate(total=sum(n)) %>%
        arrange(-total) %>%
        ungroup() %>%
        mutate(index=as.numeric(factor(event_gene, levels=unique(event_gene)))) %>%
        filter(index <= top_n) %>%
        ggbarplot(x="event_gene", y="n", fill="drug_screen", palette="jco", color=NA) + 
        labs(x="Event & Gene", y="N. Significant Associations", 
             title=sprintf("Top %s", top_n)) +
        coord_flip()
    
    # - agreement within GDSC1 and GDSC2
    drugs_oi = X %>%
        distinct(DRUG_ID, DRUG_NAME_CLEAN, drug_screen) %>%
        count(DRUG_NAME_CLEAN, drug_screen) %>%
        filter(n>1)
    correls = drugs_oi %>% 
        left_join(X, by=c("DRUG_NAME_CLEAN","drug_screen")) %>% 
        group_by(DRUG_NAME_CLEAN) %>% 
        mutate(replicate = paste0("rep",as.numeric(factor(DRUG_ID)))) %>% 
        pivot_wider(id_cols=c("DRUG_NAME_CLEAN","drug_screen","event_gene"), 
                    names_from="replicate", 
                    values_from="spldep_coefficient") %>% 
        group_by(drug_screen,DRUG_NAME_CLEAN) %>%
        summarize(correl=cor(rep1,rep2,method="pearson",use="pairwise.complete.obs"))
    
    plts[["associations-agreement_within"]] = correls %>%
        ggviolin(x="drug_screen", y="correl", trim=TRUE,
                 color=NA, fill="drug_screen", palette="jco") + 
        geom_boxplot(width=0.1, outlier.size=0.1) +
        geom_text_repel(aes(label=DRUG_NAME_CLEAN), correls, size=1, family="Arial") +
        geom_text(aes(y=1, label=n), correls %>% count(drug_screen), size=1, family="Arial") +
        guides(fill="none") +
        labs(x="Drug Screen", y="Pearson Correlation")
            
    # - overlaps with cell lines used to model splicing dependency
    # - overlaps with cell lines used in GDSC1 and GDSC2
    samples_oi = list(
        "Demeter2" = rnai %>% dplyr::select(-index) %>% colnames(),
        "GDSC1" = drug_screen %>% filter(DATASET=="GDSC1") %>% pull(ARXSPAN_ID) %>% unique() %>% na.omit(),
        "GDSC2" = drug_screen %>% filter(DATASET=="GDSC2") %>% pull(ARXSPAN_ID) %>% unique() %>% na.omit()
    )
    m = samples_oi %>% 
        list_to_matrix() %>% 
        make_comb_mat()
    plts[["associations-overlaps_screens-upset"]] = m %>%
    UpSet(comb_order = order(comb_size(m))) %>%
        draw() %>%
        grid.grabExpr() %>%
        as.ggplot()
    
    return(plts)
}


plot_rankings = function(shortest_paths, rankings, drug_targets, models, spldep_models){
    plts = list()
    
    # to extract mechanisms, we may only want to focus on drug-exon associations
    # in which splicing is more associated than gene expression
    spldep_models %>% ggscatter(x="event_coefficient_mean", y="gene_coefficient_mean", alpha=0.5)
    spldep_models %>% 
        mutate(event_coefficient_mean = abs(event_coefficient_mean),
               gene_coefficient_mean = abs(gene_coefficient_mean)) %>%
        ggscatter(x="event_coefficient_mean", y="gene_coefficient_mean", alpha=0.5)
    
    spldep_models %>% 
        mutate(diff = abs(event_coefficient_mean) - abs(gene_coefficient_mean)) %>%
        gghistogram(x="diff")
    
    
    # do significantly associated events in gene targets rank to the top of
    # the association effect size?
    ranking_real = rankings %>% 
        left_join(drug_targets %>% distinct(DRUG_ID,TARGET) %>% mutate(is_target=TRUE), 
                  by=c("DRUG_ID","GENE"="TARGET")) %>%
        arrange(combined_ranking) %>% 
        mutate(ranking_type = "real")
    
    set.seed(RANDOM_SEED)
    ranking_random = ranking_real %>%
        mutate(is_target = sample(is_target),
               ranking_type = "random")
    ranking_targets = rbind(ranking_real, ranking_random)
    
    # significant drug-event associations of drug targets tend to rank at the top
    plts[["rankings-index_vs_known_targets"]] = ranking_targets %>%
        filter(is_target) %>%
        ggviolin(x="ranking_type", y="index", color=NA, 
                 fill="ranking_type", palette=c("darkred","grey"), trim=TRUE) +
        geom_boxplot(width=0.1, outlier.size=0.1) + 
        stat_compare_means(method="wilcox.test") +
        guides(fill="none") +
        labs(x="Ranking Type", y="Combined Ranking", 
             title=sprintf("n = %s", nrow(ranking_real %>% filter(is_target))))

    X = shortest_paths %>% 
        dplyr::select(source, shortest_path_length, DRUG_ID) %>% 
        left_join(rankings, by=c("source"="GENE", "DRUG_ID")) %>%
        left_join(models %>% distinct(ID, drug_screen, event_gene, 
                                      lr_padj, pearson_correlation), 
                  by=c("ID","event_gene","drug_screen"))
    
    best_scores = X %>% 
        group_by(DRUG_NAME_CLEAN) %>%
        slice_min(combined_ranking, n=1) %>%
        mutate(is_best=TRUE) %>%
        distinct(ID, DRUG_ID, event_gene, combined_ranking, is_best) %>%
        ungroup()

    X = X %>%
        left_join(best_scores, by=c("DRUG_NAME_CLEAN","ID","DRUG_ID",
                                    "event_gene","combined_ranking")) %>%
        mutate(is_best=replace_na(is_best,FALSE)) %>%
        distinct()
    
    thresh_targets = ranking_targets %>% 
        filter(ranking_type=="real" & is_target) %>% 
        pull(index) %>% 
        median()
    
    X = X %>%
        separate(ID, c("drug_id","min_conc","max_conc"), sep="_") %>%
        mutate(DRUG_NAME_CLEAN = sprintf("%s (%s-%s)", DRUG_NAME, min_conc, max_conc))
    
    # when significant event-drug association with drug targets,
    # there can also be other good associations of non-targets
    plts[["rankings-index_vs_paths_all"]] = X %>%
        group_by(DRUG_NAME_CLEAN) %>%
        filter(any(shortest_path_length==0)) %>%
        ungroup() %>%
        ggplot(aes(x=as.factor(shortest_path_length), y=index)) +
        geom_violin(color=NA, width=1.2, fill="darkred") + 
        geom_boxplot(width=0.1, outlier.color=NA) + 
        geom_text(aes(x=as.factor(shortest_path_length), y=13000, label=n), 
                  X %>%
                    group_by(DRUG_NAME_CLEAN) %>%
                    filter(any(shortest_path_length==0)) %>%
                    ungroup() %>%
                    count(shortest_path_length),
                  size=1, family="Arial") +
        geom_hline(yintercept = thresh_targets, linetype="dashed") +
        theme_pubr() + 
        labs(x="Shortest Path Length to Drug Target", y="Combined Raking")
    
    # some indirect associations are better than those with target
    x = X %>%
        group_by(DRUG_NAME) %>%
        filter(any(shortest_path_length==0) & is_best) %>%
        ungroup() %>% 
        ## when drug has multiple targets, 
        ## consider only shortest path length to one of them
        group_by(DRUG_ID, source) %>% 
        slice_min(shortest_path_length, n=1) %>%
        ungroup()
    
    plts[["rankings-index_vs_paths_known_best"]] = x %>%
        ggbarplot(x="DRUG_NAME_CLEAN", y="index", fill="event_gene", 
                  color="drug_screen", palette="Dark3",
                  position=position_dodge(0.9)) + 
        geom_text(aes(y=2, label=event_gene), x, size=1, family="Arial") +
        geom_text(aes(label=shortest_path_length), x, size=1, family="Arial",
                  position=position_dodge(0.9)) +
        geom_hline(yintercept = thresh_targets, linetype="dashed") +
        yscale("log10", .format=TRUE) +
        coord_flip() +
        guides(fill="none") + 
        color_palette(c("white","black")) +
        labs(x="Drug", y="Combined Raking", 
             color="Drug Screen")
    
    # for drugs whose target is not significantly associated,
    # we find events that are best associated
    x = X %>%
        group_by(DRUG_NAME) %>%
        filter(!any(shortest_path_length==0) & is_best) %>%
        # slice_min(lr_padj, n=1) %>%
        ungroup()
    
    plts[["rankings-index_vs_paths_unknown_best"]] = x %>%
        ggplot(aes(x=as.factor(shortest_path_length), y=index)) +
        geom_violin(color=NA, width=1.2, fill="darkred") + 
        geom_boxplot(width=0.1, outlier.color=NA) + 
        geom_text(aes(x=as.factor(shortest_path_length), y=19000, label=n), 
                  x %>% count(shortest_path_length),
                  size=1, family="Arial") +
        geom_hline(yintercept = thresh_targets, linetype="dashed") +
        theme_pubr() + 
        labs(x="Shortest Path Length to Drug Target", y="Combined Raking")
    
    ## top associations
    x = X %>%
        group_by(DRUG_NAME) %>%
        filter(!any(shortest_path_length==0)) %>%
        filter(is_best) %>%
        filter(index<thresh_targets) %>%
        slice_min(index, n=1) %>%
        ungroup() %>%
        group_by(shortest_path_length) %>%
        ungroup() %>%
        arrange(index) %>%
        mutate(DRUG_NAME_CLEAN = as.factor(DRUG_NAME_CLEAN),
               name = reorder_within(DRUG_NAME_CLEAN, combined_ranking, 
                                     shortest_path_length))
    
    plts[["rankings-index_vs_paths_unknown_top"]] = x %>%
        ggbarplot(x="name", y="index", fill="event_gene", 
                  position=position_dodge(0.9), color="drug_screen", 
                  palette=get_palette("Dark2",length(unique(x[["DRUG_NAME"]])))) + 
        color_palette(c("black","white")) + 
        facet_wrap(~shortest_path_length, scales="free") +
        geom_text(aes(y=0.01, label=event_gene), x, size=1, family="Arial") + 
        guides(fill="none") +    
        scale_x_reordered() +
        coord_flip() +
        theme(strip.text.x = element_text(size=6, family="Arial")) +
        labs(x="Drug", y="Combined Ranking", color="Drug Screen")
    
    ## sensitivity to nutlin
    ### PUF60 splicing may indirectly control MDM2 splicing
#     drug_screen %>% 
#     x = splicing %>% filter(EVENT %in% c("HsaEX0038400","HsaEX0005606","HsaEX0051262")) %>% column_to_rownames("EVENT") %>% t() %>% as.data.frame() #%>% rownames_to_column("sampleID") %>% pivot_longer(cols = -sampleID, names_to = "EVENT", values_to = "spldep")
#     x %>% ggscatter(x="HsaEX0038400", y="HsaEX0005606") + stat_cor(method="spearman")
#     x %>% ggscatter(x="HsaEX0038400", y="HsaEX0051262") + stat_cor(method="spearman")
    
    ## top rankings with HsaEX0034998_KRAS, known event-drug synergy
    event_oi = "HsaEX0034998_KRAS"
    x = rankings %>% filter(event_gene==event_oi) %>%
        separate(ID, c("drug_id","min_conc","max_conc"), sep="_") %>%
        mutate(DRUG_NAME_CLEAN = sprintf("%s (%s-%s)", DRUG_NAME, min_conc, max_conc))
    
    plts[["rankings-event_oi-KRAS"]] = x %>%
        ggbarplot(x="DRUG_NAME_CLEAN", y="index", fill="event_gene", 
                  position=position_dodge(0.9), color="drug_screen", 
                  palette=get_palette("Dark2",length(unique(x[["DRUG_NAME"]])))) + 
        color_palette(c("black","white")) +
        coord_flip() +
        labs(x="Drug", y="Combined Ranking", color="Drug Screen")
    
    return(plts)
}


plot_preds_ic50 = function(clusters, metadata, drug_screen, estimated_response, spldep_ccle){
    # show prediction differences between seen and unseen data for a drug
    
    # NUTLIN-3A(-), TANESPIMYCIN, DOCETAXEL
    drugs_oi = c(1047)
    X = drug_screen %>%
        mutate(DRUG_NAME_CLEAN = ID) %>%
        separate(DRUG_NAME_CLEAN, c("drug_id","min_conc","max_conc"), sep="_") %>%
        mutate(DRUG_NAME_CLEAN = sprintf("%s (%s-%s)", DRUG_NAME, min_conc, max_conc)) %>%
        distinct(DRUG_NAME_CLEAN, DRUG_ID, DATASET) %>%
        filter(DRUG_ID%in%drugs_oi) %>%
        left_join(clusters, by=c("DRUG_ID", "DATASET"="drug_screen")) %>%
        left_join(metadata, by=c("index"="DepMap_ID")) %>%
        left_join(drug_screen, 
                  by=c("index"="ARXSPAN_ID", "DRUG_NAME_CLEAN", 
                       "DRUG_ID", "DATASET")) %>%
        left_join(estimated_response, 
                  by=c("DRUG_ID", "index"="sample","DATASET"="drug_screen")) %>%
        mutate(log_ic50 = log(IC50_PUBLISHED),
               in_training_set = !is.na(log_ic50),
               leiden_labels = as.factor(leiden_labels)) %>%
        left_join(spldep_ccle %>% 
                  column_to_rownames("index") %>% 
                  t() %>% as.data.frame() %>% 
                  rownames_to_column("index"), by="index") %>%
        group_by(DRUG_NAME_CLEAN) %>%
        mutate(log_ic50 = scale(log_ic50),
               max_conc = as.factor(MAX_CONC))
    
    plts = list()
    plts[["preds_ic50-real_response"]] = X %>% 
        ggscatter(x="UMAP0", y="UMAP1", color="log_ic50", alpha=0.5, size=1) + 
        scale_color_gradient2(low="blue",mid="white",high="red") + 
        facet_wrap(~DRUG_NAME_CLEAN) +
        labs(color="Scaled log(IC50)") +
        theme(strip.text.x = element_text(size=6, family="Arial"))
    
    plts[["preds_ic50-unseen"]] = X %>% 
        ggscatter(x="UMAP0", y="UMAP1", color="in_training_set", alpha=0.5, size=1) + 
        facet_wrap(~DRUG_NAME_CLEAN) +
        labs(color="In Training Set") +
        theme(strip.text.x = element_text(size=6, family="Arial"))
    
    plts[["preds_ic50-leiden"]] = X %>% 
        ggscatter(x="UMAP0", y="UMAP1", color="leiden_labels", palette="Dark3", alpha=0.5) + 
        facet_wrap(~DRUG_NAME_CLEAN) +
        theme(strip.text.x = element_text(size=6, family="Arial"))
    
     plts[["preds_ic50-primary_disease"]] = X %>% 
        ggscatter(x="UMAP0", y="UMAP1", color="primary_disease", alpha=0.5, size=1, 
                  palette=get_palette("Paired", length(unique(X[["primary_disease"]])))) + 
        facet_wrap(~DRUG_NAME_CLEAN) + 
        labs(color="Primary Disease") +
        theme(strip.text.x = element_text(size=6, family="Arial"))
    
    plts[["preds_ic50-drug_screen"]] = X %>% 
        ggscatter(x="UMAP0", y="UMAP1", color="DATASET", palette="jco", size=1, alpha=0.5) + 
        facet_wrap(~DRUG_NAME_CLEAN) +
        labs(color="Drug Screen") +
        theme(strip.text.x = element_text(size=6, family="Arial"))
    
    plts[["preds_ic50-max_conc"]] = X %>% 
        ggscatter(x="UMAP0", y="UMAP1", color="max_conc", alpha=0.5, size=1) + 
        facet_wrap(~DRUG_NAME_CLEAN) +
        labs(color="Max. Conc.") +
        theme(strip.text.x = element_text(size=6, family="Arial"))
        
    return(plts)
}


plot_moa_clusters = function(embedding, enrichment){
    plts = list()

    # are there other biological reasons for the umap clusters?
    plts[["moa_clusters-leiden-umap"]] = embedding %>%
        mutate(leiden_labels=as.factor(leiden_labels)) %>%
        ggplot(aes(x=UMAP0, y=UMAP1, color=leiden_labels)) +
        geom_scattermore(pixels=c(1000,1000), pointsize=8) +
        theme_pubr() +
        theme(aspect.ratio=1)

    plts_enrichment = sapply(names(enrichment), function(onto){
        result = enrichment[[onto]]
        res = new("compareClusterResult", compareClusterResult = result)
        plt = dotplot(res) + labs(title=onto, x="Leiden Cluster")
        return(plt)
    }, simplify=FALSE)
    names(plts_enrichment) = sprintf("moa_clusters-leiden-enrichment-%s",
                                     names(plts_enrichment))
    plts = c(plts,plts_enrichment)
    
    return(plts)
}


plot_drug_rec = function(estimated_response, drug_screen, models){
    
    # all models
    models_oi = models %>%
        filter(lr_padj<THRESH_FDR  & n_obs>THRESH_NOBS) %>% 
        distinct(ID)
    
    X = models_oi %>%
        left_join(estimated_response, by="ID") %>%
        left_join(drug_screen %>% mutate(real_ic50 = log(IC50_PUBLISHED)), 
                  by=c("drug_screen"="DATASET", "in_demeter2",
                       "DRUG_ID","ID", "sample"="ARXSPAN_ID")) %>%
        drop_na(real_ic50, predicted_ic50)
    
    corrs = X %>% 
        group_by(drug_screen, sample, in_demeter2) %>%
        summarize(n_pos = n(), # number of possible associations
                  correlation = cor(real_ic50, predicted_ic50, method="pearson")) %>%
        ungroup()
    
    plts = list()
    plts[["drug_rec-n_screens"]] = X %>% 
        count(drug_screen, sample, in_demeter2) %>%
        gghistogram(x="n") +
        facet_wrap(~drug_screen+in_demeter2) +
        theme(strip.text.x = element_text(size=6, family="Arial"))
    
    plts[["drug_rec-correlations"]] = corrs %>% 
        ggplot(aes(x=drug_screen, y=correlation, 
                   group=interaction(drug_screen,in_demeter2))) +
        geom_violin(aes(fill=in_demeter2), color=NA) +
        geom_boxplot(width=0.1, outlier.size=0.1, position=position_dodge(0.9)) +
        geom_hline(yintercept=0, linetype="dashed") +
        labs(x="Drug Screen", y="Pearson Correlation", fill="In Demeter2") +
        geom_text(aes(y=1.00, label=n), 
                  corrs %>% count(in_demeter2, drug_screen), 
                  position=position_dodge(0.9),
                  size=1, family="Arial") +
        theme_pubr() +
        fill_palette("Dark3")
    
    plts[["drug_rec-npos_vs_corrs"]] = corrs %>% 
        ggscatter(x="n_pos", y="correlation", alpha=0.5, size=1,
                  color="in_demeter2", palette="Dark3") + 
        facet_wrap(~drug_screen+in_demeter2, scales="free_x") +
        labs(x="N. Associations per Sample", y="Pearson Correlation") +
        guides(color="none") +
        theme(strip.text.x = element_text(size=6, family="Arial"))

    # best and worse correlations
    samples_oi = corrs %>% 
        group_by(drug_screen) %>%
        drop_na(correlation) %>%
        arrange(correlation) %>% 
        filter(row_number()==1 | row_number()==n()) %>%
        left_join(X, by=c("drug_screen","sample"))
    
    plts[["drug_rec-best_worse"]] = samples_oi %>% 
        ggscatter(x="real_ic50", y="predicted_ic50", size=0.5) + 
        facet_wrap(drug_screen~sample, scales="free") +
        stat_cor(method="pearson", size=2) + 
        labs(x="Real log(IC50)", y="Predicted log(IC50)") + 
        theme_pubr(border=TRUE) +
        theme(strip.text.x = element_text(size=6, family="Arial"))
    
    return(plts)
}


make_plots = function(models, spldep_ccle, drug_screen,
                      shortest_paths, rankings, drug_targets,
                      clusters, metadata, estimated_response,
                      embedding, enrichment){
    plts = list(
        plot_associations(models, spldep_ccle, drug_screen),
        plot_rankings(shortest_paths, rankings, drug_targets, models),
        #plot_preds_ic50(clusters, metadata, drug_screen, estimated_response, spldep_ccle),
        plot_moa_clusters(embedding, enrichment),
        plot_drug_rec(estimated_response, drug_screen, models)
    )
    plts = do.call(c,plts)
    return(plts)
}

make_figdata = function(models,
                        shortest_paths, rankings, drug_targets,
                        clusters, metadata, estimated_response,
                        embedding, enrichment){
    # prep
    df_enrichs = do.call(rbind,
        lapply(names(enrichment), function(e){
            res = enrichment[[e]]
            res[["ontology"]] = e
            return(res)
        })
    )
    
    figdata = list(
        "drug_event_assoc" = list(
            "model_summaries" = models,
            "shortest_paths" = shortest_paths,
            "association_rankings" = rankings,
            "drug_targets" = drug_targets,
            "pred_ic50_by_event_drug_clustering" = clusters,
            "pred_ic50_by_drug" = estimated_response,
            "assoc_coeff_clustering" = embedding,
            "assoc_coeff_gsea" = df_enrichs
        )
    )
    return(figdata)
}


save_plt = function(plts, plt_name, extension=".pdf", 
                    directory="", dpi=350, format=TRUE,
                    width = par("din")[1], height = par("din")[2]){
    print(plt_name)
    plt = plts[[plt_name]]
    if (format){
        plt = ggpar(plt, font.title=10, font.subtitle=10, font.caption=10, 
                    font.x=8, font.y=8, font.legend=8,
                    font.tickslab=6, font.family="Arial")    
    }
    filename = file.path(directory,paste0(plt_name,extension))
    save_plot(filename, plt, base_width=width, base_height=height, dpi=dpi, units="cm")
}


save_plots = function(plts, figs_dir){
    # drug-event associations
    save_plt(plts, "associations-lr_pvalues", ".pdf", figs_dir, width=5, height=5)
    save_plt(plts, "associations-lr_fdr", ".pdf", figs_dir, width=5, height=5)
    save_plt(plts, "associations-nobs_vs_coef", ".pdf", figs_dir, width=5, height=5)
    save_plt(plts, "associations-drug_counts", ".pdf", figs_dir, width=5, height=5)
    save_plt(plts, "associations-top_drug_counts", ".pdf", figs_dir, width=8, height=8)
    save_plt(plts, "associations-spldep_counts", ".pdf", figs_dir, width=5, height=5)
    save_plt(plts, "associations-top_spldep_counts", ".pdf", figs_dir, width=8, height=8)
    save_plt(plts, "associations-agreement_within", ".pdf", figs_dir, width=5, height=5)
    save_plt(plts, "associations-overlaps_screens-upset", ".pdf", figs_dir, width=8, height=8)
    
    # rankings
    save_plt(plts, "rankings-index_vs_known_targets", ".pdf", figs_dir, width=5, height=5)
    save_plt(plts, "rankings-index_vs_paths_all", ".pdf", figs_dir, width=6, height=5)
    save_plt(plts, "rankings-index_vs_paths_known_best", ".pdf", figs_dir, width=7, height=7)
    save_plt(plts, "rankings-index_vs_paths_unknown_best", ".pdf", figs_dir, width=6, height=5)
    save_plt(plts, "rankings-index_vs_paths_unknown_top", ".pdf", figs_dir, width=14, height=10)
    save_plt(plts, "rankings-event_oi-KRAS", ".pdf", figs_dir, width=7, height=5)
    
    # predictions
#     n_ = 1
#     save_plt(plts, "preds_ic50-real_response", ".pdf", figs_dir, width=4*n_, height=6)
#     save_plt(plts, "preds_ic50-unseen", ".pdf", figs_dir, width=4*n_, height=5.5)
#     save_plt(plts, "preds_ic50-leiden", ".pdf", figs_dir, width=4*n_, height=7)
#     save_plt(plts, "preds_ic50-primary_disease", ".pdf", figs_dir, width=4*n_, height=8.5)
#     save_plt(plts, "preds_ic50-drug_screen", ".pdf", figs_dir, width=4*n_, height=5.5)
#     save_plt(plts, "preds_ic50-max_conc", ".pdf", figs_dir, width=4*n_, height=5.5)
    
    # MoA
    save_plt(plts, "moa_clusters-leiden-umap", ".pdf", figs_dir, width=10, height=10)
    save_plt(plts, "moa_clusters-leiden-enrichment-GO_BP", ".pdf", figs_dir, width=18, height=15)
    
    # drug recommendations
    save_plt(plts, "drug_rec-correlations", ".pdf", figs_dir, width=5, height=5)
    save_plt(plts, "drug_rec-npos_vs_corrs", ".pdf", figs_dir, width=7, height=8)
    save_plt(plts, "drug_rec-best_worse", ".pdf", figs_dir, width=7, height=8)
}


save_figdata = function(figdata, dir){
    lapply(names(figdata), function(x){
        d = file.path(dir,'figdata',x)
        dir.create(d, recursive=TRUE)
        lapply(names(figdata[[x]]), function(nm){
            df = figdata[[x]][[nm]]
            filename = file.path(d, paste0(nm,'.tsv.gz'))
            write_tsv(df, filename)
            
            print(filename)
        })
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
        mutate(event_gene = paste0(EVENT,"_",GENE),
               event_type = gsub("Hsa","",gsub("[^a-zA-Z]", "",EVENT)))
    drug_targets = read_tsv(drug_targets_file) %>% 
        mutate(DRUG_NAME=toupper(DRUG_NAME)) %>%
        distinct(DRUG_ID,DRUG_NAME,TARGET,TARGET_PATHWAY)
    embedding = read_tsv(embedding_file) %>% mutate(DRUG_ID=as.numeric(gsub("_.*","",index)))
    estimated_response = read_tsv(estimated_response_file)
    drug_screen = load_drug_screens(drug_screens_dir)
    ontologies = load_ontologies(msigdb_dir)
    clusters = read_tsv(clusters_file)
    spldep_ccle = read_tsv(spldep_ccle_file) %>%
        filter(index %in% unique(models[["EVENT"]]))
    paths_real = read_tsv(paths_real_file) %>% drop_na()
    rnai = read_tsv(rnai_file)
    spldep_models = read_tsv(spldep_models_file) %>%
        filter(EVENT %in% unique(models[["EVENT"]]))
    
    # prep inputs
    rankings = compute_rankings(models, drug_targets)
    shortest_paths = paths_real
    estimated_response = estimated_response %>% mutate(in_demeter2 = sample %in% colnames(rnai))
    
    ## GSEA
    gene_sets = embedding %>%
        left_join(
            drug_targets %>%
            distinct(DRUG_ID,TARGET),
          by="DRUG_ID") %>%
        distinct(leiden_labels,TARGET) %>%
        get_sets("leiden_labels","TARGET")
    universe = unique(unlist(gene_sets))
    enrichment = sapply(names(gene_sets), function(gene_set){
        res = run_enrichment(gene_sets[[gene_set]], universe, ontologies)
        return(res)
    }, simplify=FALSE)
    enrichment = get_enrichment_result(enrichment)
    enrichment = enrichment[sapply(enrichment,nrow)>0]
    
    # make plots
    plts = make_plots(models, spldep_ccle, drug_screen,
                      shortest_paths, rankings, drug_targets,
                      clusters, metadata, estimated_response,
                      embedding, enrichment)
    
    # make figdata
    figdata = make_figdata(models,
                           shortest_paths, rankings, drug_targets,
                           clusters, metadata, estimated_response,
                           embedding, enrichment)
    
    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}