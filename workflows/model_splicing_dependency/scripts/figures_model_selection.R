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
# model selection
# - pvalue distributions tend to recover TP exons
# - TPR vs FDR with p-value thresholds
# - Spearman correlation with Pearson thresholds
# 
# models properties
# - exon inclusion variation
# - exon length distribution
# - protein impact
# - GSEA exons
# - correlations with mitotic and stemness indices
# - tumorigenesis
#
# models validation
# - mutation frequencies at the gene level
# - mutation frequencies at the exon level
# - selected and not selected exons from the TP set

require(optparse)
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
source(file.path(ROOT,"src","R","utils.R"))

# variables
#MIN_N_OBS = 20
THRESH_LR_PVALUE = 0.025
THRESH_CORR = 0.2
SIZE_CTL = 100
THRESH_INDICES = 0.3 # correlation with sample indices

# formatting
PAL_SINGLE_ACCENT = "orange"
PAL_SINGLE_LIGHT = "#efb300ff"#"#6AC2BF"
PAL_SINGLE_DARK = "#007e67ff"
PAL_SINGLE_NEUTRAL = "#716454"
PAL_DUAL = c(PAL_SINGLE_DARK, PAL_SINGLE_LIGHT)
PAL_FDR_DARK = "#005AB5"
PAL_FDR_LIGHT = "#DC3220"

LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

# Development
# -----------
# PREP_DIR = file.path(ROOT,"data","prep")
# RAW_DIR = file.path(ROOT,"data","raw")
# RESULTS_DIR = file.path(ROOT,"results","model_splicing_dependency")
# models_file = file.path(RESULTS_DIR,"files","models_gene_dependency-EX","model_summaries.tsv.gz")
# ccle_stats_file = file.path(PREP_DIR,"stats","CCLE.tsv.gz")
# indices_file = file.path(RESULTS_DIR,"files","correlation_spldep_indices-EX.tsv.gz")
# event_mut_file = file.path(PREP_DIR,'event_snv','CCLE-EX.tsv.gz')
# gene_mut_freq_file = file.path(ROOT,"data","prep","gene_mutation_freq","CCLE.tsv.gz")
# event_mut_freq_file = file.path(ROOT,"data","prep","event_mutation_freq","CCLE-EX.tsv.gz")
# msigdb_dir = file.path(ROOT,"data","raw","MSigDB","msigdb_v7.4","msigdb_v7.4_files_to_download_locally","msigdb_v7.4_GMTs")
# protein_impact_file = file.path(ROOT,"data","raw","VastDB","PROT_IMPACT-hg38-v3.tab.gz")
# cosmic_genes_file = file.path(ROOT,"data","raw","COSMIC","cancer_gene_census.tsv")
# spldep_file = file.path(RESULTS_DIR,"files","splicing_dependency-EX","mean.tsv.gz")
# rnai_file = file.path(PREP_DIR,"demeter2","CCLE.tsv.gz")
# genexpr_file = file.path(PREP_DIR,"genexpr_tpm","CCLE.tsv.gz")
# splicing_file = file.path(PREP_DIR,"event_psi","CCLE-EX.tsv.gz")
# cancer_events_file = file.path(ROOT,"support","cancer_events.tsv")
# ascanceratlas_file = file.path(RAW_DIR,"ASCancerAtlas","CASE_all-VastDB_mapped.tsv.gz")
# event_info_file = file.path(RAW_DIR,"VastDB","EVENT_INFO-hg38_noseqs.tsv")
# metadata_file = file.path(PREP_DIR,"metadata","CCLE.tsv.gz")
# ppi_closeness_file = file.path(RESULTS_DIR,'files','COSMIC','ppi_closeness-EX','merged.tsv.gz')
# randsel_events_file = file.path(RESULTS_DIR,'files','random_model_selection-EX-1000its','events-merged.tsv.gz')
# randsel_genes_file = file.path(RESULTS_DIR,'files','random_model_selection-EX-1000its','genes-merged.tsv.gz')
# figs_dir = file.path(RESULTS_DIR,"figures","model_selection")


##### FUNCTIONS #####
load_ontologies = function(msigdb_dir, protein_impact_file, cosmic_genes_file){
    ontologies = list(
        "hallmarks" = read.gmt(file.path(msigdb_dir,"h.all.v7.4.symbols.gmt")),
        "oncogenic_signatures" = read.gmt(file.path(msigdb_dir,"c6.all.v7.4.symbols.gmt")),
        "GO_BP" = read.gmt(file.path(msigdb_dir,"c5.go.bp.v7.4.symbols.gmt")),
        "GO_CC" = read.gmt(file.path(msigdb_dir,"c5.go.cc.v7.4.symbols.gmt")),
        "protein_impact" = read_tsv(protein_impact_file) %>%
            dplyr::rename(EVENT=EventID, term=ONTO) %>%
            dplyr::select(term,EVENT) %>%
            mutate(term_clean=gsub(" \\(.*","",term),
                   term_clean=gsub("ORF disruption upon sequence exclusion",
                                   "ORF disruption (exclusion)",term_clean),
                   term_clean=gsub("ORF disruption upon sequence inclusion",
                                   "ORF disruption (inclusion)",term_clean),
                   term_clean=gsub("In the CDS, with uncertain impact",
                                   "In the CDS (uncertain)",term_clean)),
        "cosmic" = read_tsv(cosmic_genes_file) %>%
            dplyr::select("Gene Symbol") %>%
            rename(gene = `Gene Symbol`) %>%
            mutate(term = "COSMIC_CENSUS") %>%
            dplyr::select(term,gene)
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


run_enrichment = function(genes, events, universe, ontologies){
    enrichments = list()
    enrichments[["hallmarks"]] = enricher(genes, TERM2GENE=ontologies[["hallmarks"]], universe=universe[["genes"]])
    enrichments[["oncogenic_signatures"]] = enricher(genes, TERM2GENE=ontologies[["oncogenic_signatures"]], universe=universe[["genes"]])
    enrichments[["GO_BP"]] = enricher(genes, TERM2GENE=ontologies[["GO_BP"]], universe=universe[["genes"]])
    enrichments[["GO_CC"]] = enricher(genes, TERM2GENE=ontologies[["GO_CC"]], universe=universe[["genes"]])
    enrichments[["protein_impact"]] = enricher(events, TERM2GENE=ontologies[["protein_impact"]] %>% dplyr::select(term_clean, EVENT), universe=universe[["events"]], maxGSSize=1e6)
    enrichments[["cosmic"]] = enricher(genes, TERM2GENE=ontologies[["cosmic"]], universe=universe[["genes"]], maxGSSize=1000)
    
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


run_enrichments_indices = function(indices, selected_events, event_info){
    X = indices %>% 
        left_join(event_info %>% distinct(EVENT,GENE), by=c("index"="EVENT")) %>% 
        filter(index%in%selected_events) %>%
        mutate(corr_sign = ifelse(sign(correlation)>0,"positive","negative"),
               corr_sign = ifelse(sign(correlation)==0,NA,corr_sign),
               is_sel = abs(correlation)>THRESH_INDICES)
    
    index_names = indices %>% pull(index_name) %>% unique()
    results = lapply(index_names, function(index_oi){
        result = lapply(c("positive","negative"), function(sign_oi){
            genes = X %>% 
                filter(corr_sign==sign_oi & index_name==index_oi & is_sel) %>% 
                pull(GENE) %>% 
                unique()
            
            print(length(genes))
            universe = X %>% filter(index_name==index_oi) %>% pull(GENE) %>% unique()
            res = enricher(genes, TERM2GENE=ontologies[["GO_BP"]], universe=universe)
            res = res %>% 
                as.data.frame() %>% 
                mutate(corr_sign=sign_oi, 
                       index_name=index_oi, 
                       Cluster=paste0(index_oi,"\n&\n",corr_sign))
            return(res)
        })
        result = do.call(rbind, result)
        return(result)
    })
    results = do.call(rbind, results)
    return(results)
}


thresh_eval_pvalue = function(models, ctl_neg, ctl_pos){
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
        df[["tpr"]] = df[["total_ctl_pos"]] / length(ctl_pos_corrected)
        df[["fpr"]] = df[["total_ctl_neg"]] / length(ctl_neg_corrected)
        df[["tpr_vs_fpr"]] = df[["tpr"]] - df[["fpr"]]
        return(df)
    })
    eval_pvalue = do.call(rbind,eval_pvalue)
    eval_pvalue = eval_pvalue %>% mutate(cumsums=cumsum(tpr_vs_fpr))
    return(eval_pvalue)
}


thresh_eval_corr = function(models, rnai, spldep){
    # prep
    rnai = rnai %>% filter(index%in%models[["GENE"]]) %>% column_to_rownames("index")
    spldep = spldep %>% filter(index%in%models[["EVENT"]]) %>% column_to_rownames("index")
    
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
                       method="spearman", use="pairwise.complete.obs")
            df = data.frame(
                id = sample_oi,
                corr = test[["estimate"]][["rho"]],
                pvalue = test[["p.value"]],
                thresh = thresh,
                total_events = nrow(models_filtered),
                total_genes = length(genes),
                total_ctl_neg = sum(genes %in% ctl_neg_corrected),
                total_ctl_pos = sum(models_filtered %>% pull(EVENT) %in% ctl_pos_corrected)
            )
            df[["tpr"]] = df[["total_ctl_pos"]] / length(ctl_pos_corrected)
            df[["fpr"]] = df[["total_ctl_neg"]] / length(ctl_neg_corrected)
            return(df)
        })
        tmp = do.call(rbind,tmp)
        return(tmp)
    })
    eval_corr = do.call(rbind, eval_corr)
    eval_corr = eval_corr %>% mutate(thresh_fct = as.factor(thresh))
    return(eval_corr)
}


get_rnai_stats = function(rnai, models){
    rnai = rnai %>% filter(index%in%models[["GENE"]]) %>% column_to_rownames("index")
    rnai_stats = data.frame(
            rnai_med = apply(rnai,1,median, na.rm=TRUE),
            rnai_std = apply(rnai,1,sd, na.rm=TRUE),
            n_missing = rowSums(is.na(rnai))
        ) %>% rownames_to_column("GENE")   
    
    return(rnai_stats)
}


get_spldep_stats = function(spldep, models){
    X = spldep %>% 
        filter(index %in% (models %>% filter(is_selected) %>% pull(EVENT))) %>% 
        column_to_rownames("index")
    spldep_stats = data.frame(
            avg = apply(X,1,mean, na.rm=TRUE),
            med = apply(X,1,median, na.rm=TRUE),
            std = apply(X,1,sd, na.rm=TRUE),
            min = apply(X,1,min, na.rm=TRUE),
            max = apply(X,1,max, na.rm=TRUE),
            q75 = apply(X,1,quantile, na.rm=TRUE, probs=0.75),
            q25 = apply(X,1,quantile, na.rm=TRUE, probs=0.25),
            range = abs(apply(X,1,max, na.rm=TRUE)) - abs(apply(X,1,min, na.rm=TRUE)),
            n_missing = rowSums(is.na(X))
        ) %>% rownames_to_column("EVENT")
    return(spldep_stats)
}


plot_model_selection = function(models, rnai_stats, cancer_events, 
                                eval_pvalue, eval_corr){
    plts = list()         
    
    # - pvalue distributions tend to recover TP exons    
    # as a positive control, we got ~300 exons whose modulation affects cell
    # proliferation. Of course, we need them to belong to genes that vary to
    # be likely to have an effect across cancers. So we took the top 100.
    ctl_pos_genes = cancer_events %>% pull(GENE) %>% unique()
    plts[["model_sel-deps_sorted_vs_std_ctl_pos"]] = rnai_stats %>% 
        arrange(-rnai_std) %>% 
        mutate(index=row_number(), 
               is_ctl=GENE %in% ctl_pos_genes) %>%
        ggscatter(x="index", y="rnai_std", color="is_ctl", 
                  size="is_ctl", palette=PAL_DUAL) +
        scale_size_discrete(range=c(0.5,1)) +
        labs(x="Index", y="Gene Std. Demeter2", 
             color="In Positive Control", size="In Positive Control", 
             title=sprintf("Total Control Genes = %s", length(ctl_pos_genes))) +
        theme(aspect.ratio=1)
    
    plts[["model_sel-deps_med_vs_std_ctl_pos"]] = rnai_stats %>% 
        mutate(is_ctl=GENE %in% ctl_pos_genes) %>%
        ggplot(aes(x=rnai_med, y=rnai_std, color=is_ctl)) +
        geom_scattermore(pixels = c(1000,1000), pointsize=5, alpha=0.5) +
        theme_pubr() +
        color_palette(palette = PAL_DUAL) +
        labs(x="Gene Median Demeter2", y="Gene Std. Demeter2", 
             title=sprintf("Total Control Genes = %s", length(ctl_pos_genes)),
             color="In Positive Control") +
        theme(aspect.ratio=1)

    ctl_pos_events = cancer_events %>% pull(EVENT) %>% unique()
    plts[["model_sel-lr_pvalue_ctl_pos"]] = models %>%
        mutate(is_ctl_pos = EVENT %in% ctl_pos_events) %>%
        ggviolin(x="is_ctl_pos", y="lr_pvalue", trim = TRUE, 
                 fill="is_ctl_pos", color=NA, palette = PAL_DUAL) + 
        geom_boxplot(width=0.1) +
        stat_compare_means(method="wilcox.test", family=FONT_FAMILY, size=FONT_SIZE) +
        guides(fill="none") + 
        labs(x="Is Positive Control", y="LR Test p-value")
    
    # as a negative control, we selected 100 uniform gene dependencies
    ctl_neg = rnai_stats %>% 
        filter(!(GENE %in% ctl_pos_genes)) %>% # do not include positive controls
        slice_min(order_by = rnai_std, n=SIZE_CTL) %>% # selected top 100
        pull(GENE)
    
    plts[["model_sel-deps_sorted_vs_std_ctl_neg"]] = rnai_stats %>% 
        arrange(-rnai_std) %>% 
        mutate(index=row_number(), 
               is_ctl=GENE %in% ctl_neg) %>%
        ggscatter(x="index", y="rnai_std", color="is_ctl", 
                  size="is_ctl", palette=PAL_DUAL) +
        scale_size_discrete(range=c(0.5,1)) +
        labs(x="Index", y="Gene Std. Demeter2", 
             title=sprintf("Total Control Genes = %s", length(ctl_neg)),
             color="In Negative Control", size="In Negative Control") +
        theme(aspect.ratio=1)
    
    plts[["model_sel-deps_med_vs_std_ctl_neg"]] = rnai_stats %>% 
        mutate(is_ctl=GENE %in% ctl_neg) %>%
        ggplot(aes(x=rnai_med, y=rnai_std, color=is_ctl)) +
        geom_scattermore(pixels = c(1000,1000), pointsize=5, alpha=0.5) +
        theme_pubr() +
        color_palette(palette = PAL_DUAL) +
        labs(x="Gene Median Demeter2", y="Gene Std. Demeter2", 
             title=sprintf("Total Control Genes = %s", length(ctl_neg)),
             color="In Negative Control") +
        theme(aspect.ratio=1)
        
    # we selected the ensemble of models by their ability to rank dependencies
    # using likelihood ratio tests" p-values as thresholds
    plts[["model_sel-lr_pvalue"]] = models %>%
        gghistogram(x="lr_pvalue", bins=100, fill=PAL_SINGLE_LIGHT, color=NA) +
        geom_vline(xintercept=median(models[["lr_pvalue"]], na.rm=TRUE),
                   linetype="dashed", size=LINE_SIZE) +
        labs(x="LR Test p-value", y="Count")
    
    # - TPR vs FDR with p-value thresholds
    plts[["model_sel-tpr_vs_fpr"]] = eval_pvalue %>% 
        ggbarplot(x="thresh", y="cumsums", fill=PAL_SINGLE_LIGHT, color=NA) +
        labs(x="Thresholds LR Test p-value", y="CumSum(TPR - FPR)") +
        geom_text(aes(y=0.8+0.10, label=total_events), 
                  size=FONT_SIZE, family=FONT_FAMILY, color="black") +
        geom_text(aes(y=0.8+0.05, label=total_genes), 
                  size=FONT_SIZE, family=FONT_FAMILY, color="black") + 
        theme_pubr(x.text.angle=45)
    
    # - Spearman correlation with Pearson threshold
    # find the threshold for the mean pearson correlation
    plts[["model_sel-pearson_corr"]] = models %>%
        gghistogram(x="pearson_correlation_mean", bins=100, fill=PAL_SINGLE_LIGHT, color=NA) +
        geom_vline(xintercept=median(models[["pearson_correlation_mean"]], na.rm=TRUE),
                   linetype="dashed", size=LINE_SIZE) +
        labs(x="Pearson Correlation Mean (Test Set)", y="Count")
    
    plts[["model_sel-pearson_corr_vs_spearman"]] = eval_corr %>%
        ggviolin(x="thresh_fct", y="corr", fill=PAL_SINGLE_LIGHT, color=NA, trim=TRUE) + 
        geom_boxplot(fill=NA, outlier.size=0.1, width=0.2) +
        labs(x="Thresholds Single-Model Avg. Pearson Correlation", y="Spearman Correlation") +
        theme_pubr(x.text.angle=70) +
        geom_text(
            aes(x=thresh_fct, y=1.03, label=total_events), 
            eval_corr %>% 
                distinct(thresh_fct,total_events), 
            size=FONT_SIZE-0.5, color="black", family=FONT_FAMILY) +
        geom_text(
            aes(x=thresh_fct, y=1.01, label=total_genes), 
            eval_corr %>% 
                distinct(thresh_fct,total_genes), 
            size=FONT_SIZE-0.5, color="black", family=FONT_FAMILY)
    
    return(plts)
}


plot_model_properties = function(models, enrichment, indices, indices_enrich, 
                                 spldep_stats, harm_stats, ppi_closeness){
    
    plts=list()
    
    # - exon inclusion variation    
    # with this set of models, we expect to see certain properties
    X = models %>%
        drop_na(is_selected)
    
    ## selected models come from variant events
    plts[["model_prop-selected_vs_event_std"]] = X %>%
        ggviolin(x="is_selected", y="event_std", color=NA, trim=TRUE,
                 fill="is_selected", palette=PAL_DUAL) +
        geom_boxplot(fill=NA, outlier.size = 0.1) +
        stat_compare_means(method="wilcox.test", size=FONT_SIZE, family=FONT_FAMILY) +
        guides(color="none", fill="none") +
        labs(x="Selected Model", y="Event Std.", fill="Selected Model")
    
    # - exon length distribution
    plts[["model_prop-selected_vs_event_length"]] = X %>%
        ggviolin(x="is_selected", y="LE_o", color=NA, trim=TRUE,
                 fill="is_selected", palette=PAL_DUAL) +
        geom_boxplot(fill=NA, outlier.size = 0.1) +
        stat_compare_means(method="wilcox.test", size=FONT_SIZE, family=FONT_FAMILY) +
        guides(color="none", fill="none") +
        labs(x="Selected Model", y="log10(Event Length)", fill="Selected Model") +
        geom_text(aes(y=10e3, label=lab), 
                  models %>% 
                      group_by(is_selected) %>% 
                      summarize(med=median(LE_o, na.rm=TRUE)) %>% 
                      mutate(lab=paste0('med=',med)), 
                  size=FONT_SIZE, family=FONT_FAMILY) +
        yscale("log10", .format=TRUE)
    
    # - n. selected exons per gene
    plts[["model_prop-n_exons_gene"]] = X %>%
        count(is_selected, GENE) %>%
        gghistogram(x="n", fill="is_selected", color="NA", palette=PAL_DUAL, alpha=0.9) + 
        labs(x="N. Exons per Gene", y="Count", fill="Selected Model")
           
    # protein impact of selected exons/genes
    prot_imp = models %>%
        count(term, is_selected) %>%
        group_by(is_selected) %>%
        mutate(freq = n / sum(n)) %>%
        ungroup() %>%
        group_by(term) %>%
        mutate(rel_freq = freq / sum(freq)) %>%
        ungroup() %>%
        arrange(n)
    
    prot_imp_clean = models %>%
        count(term_clean, is_selected) %>%
        group_by(is_selected) %>%
        mutate(freq = n / sum(n)) %>%
        ungroup() %>%
        group_by(term_clean) %>%
        mutate(rel_freq = freq / sum(freq)) %>%
        ungroup() %>%
        arrange(n)
    
    plts[["model_prop-protein_impact-counts"]] = prot_imp %>%
        filter(is_selected) %>%
        ggbarplot(x="term", y="n", label=TRUE, lab.size=FONT_SIZE, lab.family=FONT_FAMILY,
                  fill="is_selected", color=NA, palette=PAL_SINGLE_LIGHT) + 
        theme_pubr(x.text.angle = 45, legend="right") +
        labs(x="Splicing Impact", y="Count", fill="Selected Model")
    
    plts[["model_prop-protein_impact-freqs"]] = prot_imp %>%
        ggbarplot(x="term", y="freq", label = prot_imp %>% pull(n), 
                  lab.size=FONT_SIZE, lab.family=FONT_FAMILY,
                  fill="is_selected", color=NA, palette=PAL_DUAL,
                  position=position_dodge(0.9)) + 
        theme_pubr(x.text.angle = 45, legend="right") +
        labs(x="Splicing Impact", y="Proportion", fill="Selected Model")
    
    plts[["model_prop-protein_impact-rel_freqs"]] = prot_imp %>%
        ggbarplot(x="term", y="rel_freq", fill="is_selected", color=NA, palette=PAL_DUAL) + 
        geom_hline(yintercept=0.5, linetype="dashed", size=LINE_SIZE) + 
        theme_pubr(x.text.angle = 45, legend="right") +
        labs(x="Splicing Impact", y="Norm. Proportion", fill="Selected Model")
    
    plts[["model_prop-protein_impact_clean-counts"]] = prot_imp_clean %>%
        filter(is_selected) %>%
        ggbarplot(x="term_clean", y="n", label=TRUE, lab.size=FONT_SIZE, lab.family=FONT_FAMILY,
                  fill="is_selected", color=NA, palette=PAL_SINGLE_LIGHT) + 
        theme_pubr(x.text.angle = 45, legend="right") +
        labs(x="Splicing Impact", y="Count", fill="Selected Model")
    
    plts[["model_prop-protein_impact_clean-freqs"]] = prot_imp_clean %>%
        ggbarplot(x="term_clean", y="freq", label = prot_imp_clean %>% pull(n), 
                  lab.size=FONT_SIZE, lab.family=FONT_FAMILY,
                  fill="is_selected", color=NA, palette=PAL_DUAL,
                  position=position_dodge(0.9)) + 
        theme_pubr(x.text.angle = 45, legend="right") +
        labs(x="Splicing Impact", y="Proportion", fill="Selected Model")
    
    plts[["model_prop-protein_impact_clean-rel_freqs"]] = prot_imp_clean %>%
        ggbarplot(x="term_clean", y="rel_freq", fill="is_selected", color=NA, palette=PAL_DUAL) + 
        geom_hline(yintercept=0.5, linetype="dashed", size=LINE_SIZE) + 
        theme_pubr(x.text.angle = 45, legend="right") +
        labs(x="Splicing Impact", y="Norm. Proportion", fill="Selected Model")
    
    x = models %>%
        # the more included, the more it proliferates upon exclusion (betaPSI > 0) ~ +Prolif.
        # the more included, the less it proliferates upon exclusion (betaPSI < 0) ~ -Prolif.
        filter(is.finite(event_coefficient_mean) & is_selected) %>%
        mutate(event_coef_sign = ifelse(sign(event_coefficient_mean)>0, "+Prolif.", "-Prolif.")) %>%
        count(event_coef_sign, term_clean) %>%
        group_by(event_coef_sign) %>%
        mutate(freq = n / sum(n)) %>%
        ungroup() %>%
        group_by(term_clean) %>%
        mutate(rel_freq = freq / sum(freq)) %>%
        ungroup() %>%
        arrange(n)
    plts[["model_prop-protein_impact_clean_vs_sign_coef-freqs"]] = x %>%
        ggbarplot(x="term_clean", y="freq", label = x %>% pull(n), 
                  lab.size=FONT_SIZE, lab.family=FONT_FAMILY, 
                  fill="event_coef_sign", color=NA, palette="uchicago",
                  position=position_dodge(0.9)) + 
        theme_pubr(x.text.angle = 45, legend="right") +
        labs(x="Splicing Impact", y="Proportion", fill="Upon Exclusion")
    
    # GSEA of selected exons/genes
    plts_enrichment = lapply(names(enrichment), function(e_name){
        res = enrichment[[e_name]]
        plts = list()
        plts[["dotplot"]] = dotplot(res) + 
            scale_size(range=c(0.5,3)) + 
            scale_color_continuous(
                low=PAL_FDR_LIGHT, high=PAL_FDR_DARK, 
                name="FDR", guide=guide_colorbar(reverse=TRUE)) +
            theme_pubr()
        
        plts[["cnetplot"]] = cnetplot(res, cex_label_category=0.5, 
                                      cex_label_gene=0.5)
        names(plts) = sprintf("%s-%s",e_name,names(plts))
        return(plts)
    })
    plts_enrichment = do.call(c,plts_enrichment)
    names(plts_enrichment) = sprintf("model_prop-enrichment-%s",
                                     names(plts_enrichment))
    plts = c(plts, plts_enrichment)
    
    # are selected exons associated with transcriptomic indices?
    X = models %>%
        left_join(indices, by=c("EVENT"="index")) %>% 
        mutate(corr_sign = ifelse(sign(correlation)>0,"Positive","Negative")) %>%
        drop_na(correlation, is_selected)
    
    plts[["model_prop-indices-violin"]] = X %>% 
        filter(is_selected) %>%
        ggplot(aes(x=corr_sign, y=correlation)) + 
        geom_violin(aes(fill=is_selected), color=NA) + 
        geom_boxplot(fill=NA, outlier.size=0.1, 
                     width=0.1, position=position_dodge(0.9)) +
        fill_palette(PAL_SINGLE_LIGHT) + 
        geom_hline(yintercept=c(-THRESH_INDICES,0,THRESH_INDICES), linetype="dashed", size=LINE_SIZE) + 
        theme_pubr() + 
        labs(x="Correlation Sign", y="Spearman Correlation", fill="Selected Model")  +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        geom_text(aes(y=0.8, label=n), 
                  X %>% 
                    filter(is_selected) %>%
                    count(corr_sign),
                  size=FONT_SIZE, family=FONT_FAMILY)
    
    res = new("compareClusterResult", compareClusterResult = indices_enrich)
    plts[["model_prop-indices-enrichment-GO_BP-dotplot"]] = res %>% 
        dotplot() + 
        scale_size(range=c(0.5,3)) + 
            scale_size(range=c(0.5,3)) + 
            scale_color_continuous(
                low=PAL_FDR_LIGHT, high=PAL_FDR_DARK, 
                name="FDR", guide=guide_colorbar(reverse=TRUE)) +
            theme_pubr()
    
    plts[["model_prop-indices-top_pos"]] = X %>% 
        filter(is_selected) %>%
        group_by(index_name) %>%
        slice_max(correlation, n=15) %>%
        ungroup() %>%
        mutate(event_gene = as.factor(event_gene),
               name = reorder_within(event_gene, correlation, index_name)) %>%
        ggbarplot(x="name", y="correlation", fill=PAL_SINGLE_LIGHT, color=NA) +
        facet_wrap(~index_name, scales="free") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        scale_x_reordered() +
        labs(x="Event & Gene", y="Spearman Correlation") +
        coord_flip()
    
    plts[["model_prop-indices-top_neg"]] = X %>% 
        filter(is_selected) %>%
        group_by(index_name) %>%
        slice_max(-correlation, n=15) %>%
        ungroup() %>%
        mutate(event_gene = as.factor(event_gene),
               name = reorder_within(event_gene, correlation, index_name)) %>%
        ggbarplot(x="name", y="correlation", fill=PAL_SINGLE_LIGHT, color=NA) +
        facet_wrap(~index_name, scales="free") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        scale_x_reordered() +
        labs(x="Event & Gene", y="Spearman Correlation") +
        coord_flip()
    
    # - tumorigenesis
    #     there are exons with different behaviors dependencing 
    #     on the splicing dependency:
    #       - oncoexons: SplDep<0
    #       - tumor suppressor exons: SplDep>0
    #       - double agent exons: SplDep ambiguous
    X = spldep_stats %>%
        left_join(models[c("EVENT","event_gene","term","term_clean","is_selected")], 
                  by="EVENT")

    plts[["model_prop-tumorigenesis-scatter"]] = X %>% 
        ggplot(aes(x=q25, y=q75)) +
        geom_scattermore(aes(color=med), pointsize=5, pixels=c(1000,1000), alpha=1) +
        scale_color_gradient2(low='#0073C2FF', mid="lightgray", high='#EFC000FF', midpoint=0) +
        theme_pubr() +
        theme(aspect.ratio=1) +
        geom_hline(yintercept=0, linetype="dashed", size=LINE_SIZE) + 
        geom_vline(xintercept=0, linetype="dashed", size=LINE_SIZE) +
        labs(x="0.25 Quantile", y="0.75 Quantile") +
        geom_text_repel(aes(label=event_gene), size=FONT_SIZE, 
                        segment.size=0.1, family=FONT_FAMILY, 
                        X %>% slice_min(order_by=q75, n=10), max.overlaps=50) +
        geom_text_repel(aes(label=event_gene), size=FONT_SIZE, 
                        segment.size=0.1, family=FONT_FAMILY, 
                        X %>% slice_max(order_by=q75, n=10), max.overlaps=50)
    
    X = X %>% 
        left_join(indices, by=c("EVENT"="index"))
    
    plts[["model_prop-tumorigenesis_vs_indices"]] = X %>% 
        ggscatter(x="correlation", y="med", alpha=0.5, size=1, color=PAL_SINGLE_LIGHT) + 
        stat_cor(method="spearman", label.y.npc = "bottom", 
                 size=FONT_SIZE, family=FONT_FAMILY) + 
        geom_hline(yintercept=0, linetype="dashed", size=LINE_SIZE) + 
        geom_vline(xintercept=0, linetype="dashed", size=LINE_SIZE) +
        facet_wrap(~index_name) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) + 
        geom_text_repel(aes(label=event_gene),
                        X %>% 
                            group_by(index_name) %>% 
                            slice_max(abs(correlation)*abs(med), n=5),
                        size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1) +
        labs(x="Spearman Correlation", y="Median(Spl. Dep.)")
    
    ord = prot_imp_clean %>% 
        filter(is_selected) %>% 
        pull(term_clean)
    plts[["model_prop-tumorigenesis_vs_prot_imp"]] = X %>%
        ggboxplot(x="term_clean", y="med", fill=PAL_SINGLE_LIGHT, outlier.size=0.1, order=ord) +
        theme_pubr(x.text.angle = 45, legend="right") +
        labs(x="Splicing Impact", y="med(Spl. Dep.)", fill="Selected Model")
    
    # harm scores
    X = harm_stats %>%
        left_join(spldep_stats, by="EVENT", suffix=c("_harm","_spldep")) %>%
        left_join(models[c("EVENT","event_gene","term","term_clean","is_selected")], by="EVENT") %>%
        mutate(harm_type = ifelse(med_spldep<0,"Exclusion","Inclusion"))

    plts[["model_prop-harm_scores-scatter"]] = X %>% 
        ggplot(aes(x=q25_harm, y=q75_harm)) +
        geom_scattermore(aes(color=med_harm), pointsize=5, pixels=c(1000,1000), alpha=1) +
        scale_color_gradient2(low='darkred', mid="white", midpoint=0) +
        theme_pubr() +
        theme(aspect.ratio=1) +
        geom_hline(yintercept=0, linetype="dashed", size=LINE_SIZE) + 
        geom_vline(xintercept=0, linetype="dashed", size=LINE_SIZE) +
        labs(x="0.25 Quantile", y="0.75 Quantile") +
        geom_text_repel(aes(label=event_gene), size=FONT_SIZE, 
                        segment.size=0.1, family=FONT_FAMILY, 
                        X %>% group_by(harm_type) %>% slice_min(order_by=q75_harm, n=5), 
                        max.overlaps=50) +
        facet_wrap(~harm_type, scales="free_y") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    plts[["model_prop-harm_scores_vs_prot_imp"]] = X %>% 
        ggboxplot(x="term_clean", y="med_harm", fill=PAL_SINGLE_LIGHT, outlier.size=0.1, order=ord) +
        theme_pubr(x.text.angle = 45, legend="right") +
        labs(x="Splicing Impact", y="med(Spl. Dep.)", fill="Selected Model") +
        facet_wrap(~harm_type, scales="free_y") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    # are genes bearing selected exons closer to oncogenes?
    X = ppi_closeness %>%
        count(type, shortest_path_length) %>%
        group_by(type) %>%
        mutate(freq = n / sum(n)) %>%
        ungroup()
    
    test = compare_means(
        shortest_path_length ~ type, data=ppi_closeness,
        method = "kruskal.test"
    ) %>%
    mutate(label = paste0(method,", = ",p.adj))
    
    plts[["model_prop-ppi_closeness"]] = X %>% 
        ggbarplot(x="shortest_path_length", y="freq", fill="type", 
                  color=NA, position=position_dodge(0.7), 
                  palette=c("grey", PAL_SINGLE_LIGHT)) +
        labs(x="Shortest Path Length", y="Rel. Proportion", fill="Dataset Type") +
        geom_text(aes(x=1.7, y=0.7, label=label), test, size=FONT_SIZE, family=FONT_FAMILY)
    
    # consider type of splicing impact
    X = ppi_closeness %>%
        left_join(models %>% distinct(GENE, term_clean), by=c("sink"="GENE")) %>%
        count(term_clean, type, shortest_path_length) %>%
        group_by(term_clean, type) %>%
        mutate(freq = n / sum(n)) %>%
        ungroup() %>%
        drop_na()
    
    test = lapply(unique(models[["term_clean"]]), function(term_i){
        
        df = tryCatch(
                ppi_closeness %>%
                left_join(models %>% distinct(GENE, term_clean), by=c("sink"="GENE")) %>%
                filter(term_clean == term_i) %>%
                compare_means(
                    shortest_path_length ~ type, data=.,
                    method = "kruskal.test"
                ) %>%
                mutate(label = paste0(method,", = ",p.adj),
                       term_clean = term_i),
            
                error = function(e){
                    df = data.frame(
                        ".y." = NA,
                        "p" = NA,
                        "p.adj" = NA,
                        "p.format" = NA,
                        "p.signif" = NA,
                        "method" = NA,
                        "label" = NA,
                        "term_clean" = term_i
                    )
                    
                    return(df)
                }
            )
        
        return(df)
    })
    test = do.call(rbind, test)
    
    plts[["model_prop-ppi_closeness_by_impact"]] = X %>% 
        ggbarplot(x="shortest_path_length", y="freq", fill="type", 
                  color=NA, position=position_dodge(0.7), 
                  palette=c("grey", PAL_SINGLE_LIGHT)) +
        labs(x="Shortest Path Length", y="Rel. Proportion", fill="Dataset Type") +
        geom_text(aes(x=1.7, y=0.7, label=label), test, size=FONT_SIZE, family=FONT_FAMILY) + 
        facet_wrap(~term_clean)
    
    return(plts)
}


plot_model_validation = function(models, gene_mut_freq, event_mut_freq, randsel_genes, randsel_events){
    plts = list()
    
    # - mutation frequencies at the gene level
    ## selected models are in genes less prone to have deleterious mutations
    X = models %>%
        left_join(gene_mut_freq, by="GENE") %>%
        group_by(Variant_Classification,GENE,mut_freq_per_kb) %>%
        summarize(is_selected = any(is_selected)) %>%
        drop_na() %>%
        ungroup()
    
    ## filter out variant types without less than 10 observations per group
    variants_oi = X %>%
        count(Variant_Classification, is_selected) %>%
        filter(n>10) %>%
        count(Variant_Classification) %>%
        filter(n>1) %>%
        pull(Variant_Classification)
    
    plts[["model_val-mutation_gene_count"]] = X %>% 
        filter(Variant_Classification %in% variants_oi) %>%
        count(Variant_Classification, is_selected) %>%
        ggbarplot(x="Variant_Classification", y="n", 
                  label=TRUE, lab.size=FONT_SIZE, lab.family=FONT_FAMILY, 
                  palette=PAL_DUAL, fill="is_selected", color=NA, 
                  position=position_dodge(0.9)) + 
        yscale("log10", .format=TRUE) + 
        labs(x="Mutation Effect", y="No. Genes", fill="Selected Model") +
        theme_pubr(x.text.angle=70)
    
    plts[["model_val-mutation_gene_frequency"]] = X %>% 
        filter(Variant_Classification %in% variants_oi) %>%
        ggplot(aes(x=Variant_Classification, y=mut_freq_per_kb, 
                   group=interaction(Variant_Classification,is_selected))) +
        geom_boxplot(aes(fill=is_selected), outlier.size=0.1, 
                     position=position_dodge(0.7)) +
        stat_compare_means(aes(group=is_selected), method="wilcox.test", 
                           label="p.signif", size=FONT_SIZE, family=FONT_FAMILY) +
        yscale("log10", .format=TRUE) + 
        fill_palette(PAL_DUAL) +
        labs(x="Mutation Effect", y="log10(Mut. Freq. per Kb)", fill="Selected Model") +
        theme_pubr(x.text.angle=70)
    
    # normalize with respect to silent mutations
    gene_mut_freq_null = X %>%
        distinct(GENE, Variant_Classification, is_selected, mut_freq_per_kb) %>%
        filter(Variant_Classification == "Silent") %>%
        group_by(is_selected) %>%
        mutate(null_mut_freq = median(mut_freq_per_kb, na.rm=TRUE)) %>%
        ungroup() %>%
        distinct(is_selected, null_mut_freq)
    
    x = X %>%
        left_join(gene_mut_freq_null, by="is_selected") %>%
        mutate(mut_freq_norm = mut_freq_per_kb / null_mut_freq)
    
    plts[["model_val-mutation_gene_frequency_silent_norm"]] = x %>% 
        filter(Variant_Classification %in% variants_oi) %>%
        ggplot(aes(x=Variant_Classification, y=mut_freq_norm, 
                   group=interaction(Variant_Classification,is_selected))) +
        geom_boxplot(aes(fill=is_selected), outlier.size=0.1, 
                     position=position_dodge(0.7)) +
        stat_compare_means(aes(group=is_selected), method="wilcox.test", 
                           label="p.signif", size=FONT_SIZE, family=FONT_FAMILY) +
        yscale("log10", .format=TRUE) + 
        fill_palette(PAL_DUAL) +
        labs(x="Mutation Effect", y="log10(Mut. Freq. per Kb) Norm.", fill="Selected Model") +
        theme_pubr(x.text.angle=70)
    
    # normalize by random selection selected and not selected
#     gene_mut_freq_null_selected = randsel_genes %>% 
#         left_join(X %>% distinct(Variant_Classification, mut_freq_per_kb, GENE), by="GENE") %>%
#         group_by(random_iteration, Variant_Classification) %>%
#         summarize(null_mut_freq = median(mut_freq_per_kb),
#                   n = n()) %>%
#         mutate(is_selected = TRUE) %>%
#         ungroup()
    
#     gene_mut_freq_null_notsel = X %>% 
#         distinct(GENE) %>%
#         mutate(count = 1000) %>%
#         uncount(count) %>%
#         group_by(GENE) %>%
#         mutate(random_iteration = paste0("it",row_number()-1)) %>%
#         ungroup() %>%
#         left_join(randsel_genes %>% mutate(is_selected = TRUE), 
#                   by=c("GENE","random_iteration")) %>%
#         filter(is.na(is_selected)) %>%
#         left_join(X %>% distinct(Variant_Classification, mut_freq_per_kb, GENE), by="GENE") %>%
#         group_by(random_iteration, Variant_Classification) %>%
#         summarize(null_mut_freq = median(mut_freq_per_kb),
#                   n = n()) %>%
#         mutate(is_selected = FALSE)
    
#     gene_mut_freq_null = rbind(gene_mut_freq_null_selected, gene_mut_freq_null_notsel)
    
#     x = X %>%
#         left_join(gene_mut_freq_null, by="is_selected") %>%
#         mutate(mut_freq_norm = mut_freq_per_kb / null_mut_freq)
    
#     plts[["model_val-mutation_gene_frequency_random_norm"]] = x %>% 
#         ggplot(aes(x=Variant_Classification, y=mut_freq_norm, 
#                    group=interaction(Variant_Classification,is_selected))) +
#         geom_boxplot(aes(fill=is_selected), outlier.size=0.1, 
#                      position=position_dodge(0.7)) +
#         stat_compare_means(aes(group=is_selected), method="wilcox.test", 
#                            label="p.signif", size=FONT_SIZE, family=FONT_FAMILY) +
#         yscale("log10", .format=TRUE) + 
#         fill_palette(PAL_DUAL) +
#         labs(x="Mutation Effect", y="log10(Mut. Freq. per Kb) Norm.", fill="Selected Model") +
#         theme_pubr(x.text.angle=70)
    
    # - mutation frequencies at the exon level
    # how often do selected exons get hit when the gene is mutated?
    X = models %>%
        left_join(event_mut_freq, by=c("EVENT","GENE")) %>%
        # keep only genes with a selected exon
        group_by(GENE) %>%
        filter(any(is_selected)) %>%
        ungroup() %>%
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
        left_join(event_mut_freq, by=c("EVENT","GENE")) %>%
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
    
    plts[["model_val-mutation_event_count"]] = X %>% 
        filter(Variant_Classification %in% variants_oi) %>%
        count(dataset, Variant_Classification, is_selected) %>%
        ggbarplot(x="Variant_Classification", y="n", 
                  label=TRUE, palette=PAL_DUAL, 
                  lab.size=FONT_SIZE, lab.family=FONT_FAMILY,
                  fill="is_selected", color=NA, position=position_dodge(0.9)) + 
        yscale("log10", .format=TRUE) + 
        labs(x="Mutation Effect", y="No. Events", fill="Selected Model") +
        theme_pubr(x.text.angle=70) +
        facet_wrap(~dataset, ncol=1) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    
    plts[["model_val-mutation_event_frequency"]] = X %>% 
        filter(Variant_Classification %in% variants_oi) %>%
        filter(is_selected) %>%
        ggplot(aes(x=Variant_Classification, y=fc_mut_freq, 
                   group=interaction(Variant_Classification,dataset))) +
        geom_boxplot(aes(fill=dataset), outlier.size=0.1, 
                     position=position_dodge(0.7)) +
        fill_palette(PAL_DUAL) + 
        geom_hline(yintercept=0, linetype="dashed", size=LINE_SIZE) +
        stat_compare_means(aes(group=dataset), method="wilcox.test", 
                           label="p.signif", size=FONT_SIZE, family=FONT_FAMILY) +
        labs(x="Mutation Effect", y="log2(FC Mut. Freq. per Kb)", 
             fill="Selected Model") +
        theme_pubr(x.text.angle=70)
    
    
    plts[["model_val-mutation_event_frequency-by_protein_impact"]] = X %>% 
        filter(Variant_Classification %in% variants_oi) %>%
        drop_na(fc_mut_freq) %>%
        filter(is_selected) %>%
        ggplot(aes(x=Variant_Classification, y=fc_mut_freq, 
                   group=interaction(Variant_Classification,dataset))) +
        geom_boxplot(aes(fill=dataset), outlier.size=0.1, 
                     position=position_dodge(0.7)) +
        fill_palette(PAL_DUAL) + 
        geom_hline(yintercept=0, linetype="dashed", size=LINE_SIZE) +
        stat_compare_means(aes(group=dataset), method="wilcox.test", 
                           label="p.signif", size=FONT_SIZE, family=FONT_FAMILY) +
        labs(x="Mutation Effect", y="log2(FC Mut. Freq. per Kb)", 
             fill="Selected Model") +
        theme_pubr(x.text.angle=70) +
        facet_wrap(~term_clean) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    ## alternative protein isoforms have a higher mutation frequency than expected,
    ## what exons are causing that?
#     X %>% 
#         drop_na(fc_mut_freq) %>% 
#         filter(dataset=="Real" & 
#                term_clean=="Alternative protein isoforms" & 
#                Variant_Classification=="Frame_Shift_Del" & 
#                is_selected) %>% 
#         arrange(fc_mut_freq)
    
    # Are selected events and genes mutated more/less frequently than by random chance?
    ## genes
#     X = randsel_genes %>%
#         # the random dataset
#         left_join(
#             gene_mut_freq %>% distinct(GENE, mut_freq_per_kb, Variant_Classification), 
#         by="GENE") %>%
#         mutate(type = "Random") %>%
#         bind_rows(
#             # the real dataset
#             gene_mut_freq %>%
#             distinct(GENE, mut_freq_per_kb, Variant_Classification) %>%
#             left_join(models %>% distinct(GENE, is_selected), by="GENE") %>%
#             filter(is_selected) %>%
#             mutate(type = "Real")
#         ) %>%
#         drop_na(Variant_Classification)
    
#     plts[["model_val-mutation_gene_frequency_vs_random"]] = X %>% 
#         ggboxplot(x="Variant_Classification", y="mut_freq_per_kb", 
#                   fill="type", outlier.size=0.1, palette=c("grey",PAL_SINGLE_LIGHT)) +
#         yscale("log10", .format=TRUE) +
#         stat_compare_means(aes(group=type), label.y=log10(20), 
#                            method="wilcox.test", label="p.signif", 
#                            size=FONT_SIZE, family=FONT_FAMILY) +
#         geom_text(aes(y=75, label=n), 
#                   X %>% filter(type=="Real") %>% count(Variant_Classification), 
#                   size=FONT_SIZE, family=FONT_FAMILY) +
#         labs(x="Mutation Effect", y="log10(Mut. Freq. per Gene Kb)", fill="Dataset Type") +
#         theme_pubr(x.text.angle=70)
    
#     ## events
#     X = randsel_events %>%
#         # the random dataset
#         left_join(
#             event_mut_freq %>% distinct(EVENT, event_mut_freq_per_kb, Variant_Classification), 
#         by="EVENT") %>%
#         mutate(type = "Random") %>%
#         bind_rows(
#             # the real dataset
#             event_mut_freq %>%
#             distinct(EVENT, event_mut_freq_per_kb, Variant_Classification) %>%
#             left_join(models %>% distinct(EVENT, is_selected), by="EVENT") %>%
#             filter(is_selected) %>%
#             mutate(type = "Real")
#         ) %>%
#         drop_na(Variant_Classification)
    
#     plts[["model_val-mutation_event_frequency_vs_random"]] = X %>% 
#         ggboxplot(x="Variant_Classification", y="event_mut_freq_per_kb", 
#                   fill="type", outlier.size=0.1, palette=c("grey",PAL_SINGLE_LIGHT)) +
#         yscale("log10", .format=TRUE) +
#         stat_compare_means(aes(group=type), label.y=log10(500), 
#                            method="wilcox.test", label="p.signif", 
#                            size=FONT_SIZE, family=FONT_FAMILY) +
#         geom_text(aes(y=1500, label=n), 
#                   X %>% filter(type=="Real") %>% count(Variant_Classification), 
#                   size=FONT_SIZE, family=FONT_FAMILY) +
#         labs(x="Mutation Effect", y="log10(Mut. Freq. per Event Kb)", fill="Dataset Type") +
#         theme_pubr(x.text.angle=70)
    
#     x = X %>% 
#         mutate(random_iteration = replace_na(random_iteration, "real")) %>% 
#         group_by(Variant_Classification, random_iteration) %>% 
#         summarize(med = median(event_mut_freq_per_kb, na.rm=TRUE)) %>% 
#         ungroup()
    
#     x %>% ggviolin(x="Variant_Classification", y="med") + yscale("log10") + geom_point(data=x %>% filter(random_iteration=="real")) + theme_pubr(x.text.angle = 70)
    
    return(plts)
}


plot_mutation_distances = function(event_mut, margin=500){
    X = event_mut %>%
        group_by(Variant_Classification) %>%
        mutate(
            on_event = sign(distance_to_3ss) != sign(distance_to_5ss),
            var_class_lab = sprintf("%s (n=%s)", Variant_Classification, n())
        ) %>%
        filter(sum(is_selected)>5) %>%
        ungroup() 
    
    plts = list()
    
    plts[["mutation_dists-closest_ss-distrs"]] = X %>%
        filter(abs(distance_to_closest_ss)<margin) %>% 
        ggdensity(x="distance_to_closest_ss", color="is_selected", fill=NA) + 
        geom_text(
            aes(x=-250, y=0.004, label=label, color=is_selected), 
            . %>% count(Variant_Classification, is_selected) %>%
            mutate(label = sprintf("%s (n=%s)", is_selected, n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        color_palette(PAL_DUAL) + 
        facet_wrap(~Variant_Classification, scales="free_y") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY), aspect.ratio=1) +
        labs(x="Distance to Closest Splice Site", y="Density", color="Selected Model")
    
    plts[["mutation_dists-5ss-distrs"]] = X %>% 
        filter(abs(distance_to_5ss)<margin) %>% 
        ggdensity(x="distance_to_5ss", color="is_selected", fill=NA) + 
        geom_text(
            aes(x=-250, y=0.004, label=label, color=is_selected), 
            . %>% count(Variant_Classification, is_selected) %>%
            mutate(label = sprintf("%s (n=%s)", is_selected, n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        color_palette(PAL_DUAL) + 
        facet_wrap(~Variant_Classification, scales="free_y") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY), aspect.ratio=1) +
        labs(x="Distance to Closest Splice Site", y="Density", color="Selected Model")
    
    plts[["mutation_dists-3ss-distrs"]] = X %>% 
        filter(abs(distance_to_3ss)<margin) %>% 
        ggdensity(x="distance_to_3ss", color="is_selected", fill=NA) + 
        geom_text(
            aes(x=-250, y=0.004, label=label, color=is_selected), 
            . %>% count(Variant_Classification, is_selected) %>%
            mutate(label = sprintf("%s (n=%s)", is_selected, n)),
            size=FONT_SIZE, family=FONT_FAMILY
        ) +
        color_palette(PAL_DUAL) + 
        facet_wrap(~Variant_Classification, scales="free_y") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY), aspect.ratio=1) +
        labs(x="Distance to Closest Splice Site", y="Density", color="Selected Model")
    
    return(plts)
}

plot_event_oi = function(event_oi, gene_oi, ensembl_oi, 
                         rnai, spldep, splicing, genexpr, patt){
    
    spldep_oi = spldep %>% filter(index%in%event_oi) %>% column_to_rownames("index") %>% t() %>% as.data.frame() %>% rownames_to_column("sampleID")
    rnai_oi = rnai %>% filter(index%in%gene_oi) %>% column_to_rownames("index") %>% t() %>% as.data.frame() %>% rownames_to_column("sampleID")
    splicing_oi = splicing %>% filter(EVENT%in%event_oi) %>% column_to_rownames("EVENT") %>% t() %>% as.data.frame() %>% rownames_to_column("sampleID")
    genexpr_oi = genexpr %>% filter(ID%in%ensembl_oi) %>% column_to_rownames("ID") %>% t() %>% as.data.frame() %>% rownames_to_column("sampleID")
    
    X = rnai_oi %>%
        left_join(splicing_oi, by="sampleID") %>%
        left_join(genexpr_oi, by="sampleID") %>%
        left_join(spldep_oi, by="sampleID")
    colnames(X) = c("sampleID","rnai","splicing","genexpr","spldep")
    
    plt_title = sprintf("%s_%s (%s)", event_oi, gene_oi, ensembl_oi)
    plts = list()
    plts[["spldep_vs_rnai"]] = X %>% 
        ggscatter(x="spldep", y="rnai", size=1, alpha=0.5, color=PAL_SINGLE_NEUTRAL) + 
        stat_cor(method="pearson", label.y.npc="top", 
                 label.x.npc = "middle", family=FONT_FAMILY, size=FONT_SIZE) + 
        labs(title=plt_title, x="Predicted Dep.", y="Real Dep.")
    
    plts[["splicing_vs_rnai"]] = X %>% 
        ggscatter(x="splicing", y="rnai", size=1, alpha=0.5, color=PAL_SINGLE_NEUTRAL) + 
        stat_cor(method="spearman", label.y.npc="top", 
                     label.x.npc = "middle", family=FONT_FAMILY, size=FONT_SIZE) + 
        labs(title=plt_title, x="Splicing (PSI)", y="Real Dep.")
    
    plts[["genexpr_vs_rnai"]] = X %>% 
        ggscatter(x="genexpr", y="rnai", size=1, alpha=0.5, color=PAL_SINGLE_NEUTRAL) + 
        stat_cor(method="spearman", label.y.npc="top", 
                 label.x.npc = "middle", family=FONT_FAMILY, size=FONT_SIZE) + 
        labs(title=plt_title, x="mRNA levels (log2(TPM+1))", y="Real Dep.")
    
    plts[["genexpr_vs_splicing"]] = X %>% 
        ggscatter(x="genexpr", y="splicing", size=1, alpha=0.5, color=PAL_SINGLE_NEUTRAL) + 
        stat_cor(method="spearman", label.y.npc="top", 
                 label.x.npc = "middle", family=FONT_FAMILY, size=FONT_SIZE) + 
        labs(title=plt_title, x="mRNA levels (log2(TPM+1))", y="Splicing (PSI)")
    
    
    names(plts) = paste0(patt,"-",names(plts))
    
    return(plts)
}


plot_events_oi = function(models, cancer_events, rnai, spldep, splicing, genexpr, metadata){
    # which events did we select from the positive set?
    events_oi = models %>% filter(is_selected) %>% pull(EVENT)
    X = cancer_events %>% mutate(is_selected=EVENT%in%events_oi) 
    X %>% count(is_ctl_pos, is_selected) %>% print()
    X %>% filter(is_selected) %>% print()
    
    avail_samples = intersect(colnames(rnai),colnames(splicing))
    samples_oi = metadata %>% filter(primary_disease=="Lung Cancer") %>% pull(DepMap_ID)
    samples_oi = intersect(samples_oi, avail_samples)
    
    # Examples from selected and not selected exons from the TP set
    plts = list(
        # selected
        ## KRAS
        plot_event_oi("HsaEX0034998", "KRAS", "ENSG00000133703", 
                      rnai, spldep, splicing, genexpr, "KRAS"),
        ## SMNDC1
        plot_event_oi("HsaEX0060482", "SMNDC1", "ENSG00000119953", 
                      rnai, spldep, splicing, genexpr, "SMNDC1"),
        # not selected
        ## NUMB (not included in positive set)
        plot_event_oi("HsaEX0044216", "NUMB", "ENSG00000133961", 
                      rnai, spldep, splicing, genexpr, "NUMB"),
        ## NUMB (only lung)
        plot_event_oi("HsaEX0044216", "NUMB", "ENSG00000133961", 
                      rnai[,c("index",samples_oi)], 
                      spldep[,c("index",samples_oi)], 
                      splicing[,c("EVENT",samples_oi)], 
                      genexpr[,c("ID",samples_oi)], "NUMB_lung")
    )
    plts = do.call(c,plts)
    names(plts) = paste0("events_oi-",names(plts))

    return(plts)
}


make_plots = function(models, rnai_stats, cancer_events, 
                      eval_pvalue, eval_corr, 
                      enrichment, indices, indices_enrich, spldep_stats, harm_stats, ppi_closeness,
                      gene_mut_freq, event_mut_freq, randsel_genes, randsel_events,
                      rnai, spldep, splicing, genexpr, metadata){
    plts = list(
        plot_model_selection(models, rnai_stats, cancer_events, eval_pvalue, eval_corr),
        plot_model_properties(models, enrichment, indices, indices_enrich, 
                              spldep_stats, harm_stats, ppi_closeness),
        plot_model_validation(models, gene_mut_freq, event_mut_freq, randsel_genes, randsel_events),
        plot_mutation_distances(event_mut),
        plot_events_oi(models, cancer_events, rnai, spldep, splicing, genexpr, metadata)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(models, rnai_stats, cancer_events, 
                      eval_pvalue, eval_corr, 
                      enrichment, indices, indices_enrich, spldep_stats, harm_stats,
                      gene_mut_freq, event_mut_freq, 
                      ppi_closeness, randsel_genes, randsel_events,
                      rnai, spldep, splicing, genexpr, metadata){
    # prep enrichments
    df_enrichs = do.call(rbind,
        lapply(names(enrichment), function(e){
            res = as.data.frame(enrichment[[e]])
            res[["ontology"]] = e
            return(res)
        })
    ) %>% filter(ontology=="GO_BP")
    
    figdata = list(
#         "datasets_CCLE" = list(
#             "demeter2" = rnai,
#             "gene_tpm" = genexpr,
#             "splicing_psi" = splicing,
#             "splicing_dependency" = spldep
#         ),
        "model_selection" = list(
            "model_summaries"= models,
            "rnai_stats" = rnai_stats,
            "cancer_events" = cancer_events,
            "evaluation_pvalue" = eval_pvalue,
            "evaluation_correlation" = eval_corr
        ),
        "model_properties" = list(
            "gsoa_selected" = df_enrichs,
            "correlations_transcriptomic_indices" = indices,
            "gsea_corrs_transcriptomic_indices" = indices_enrich,
            "splicing_dependecy_stats" = spldep_stats
        ),
        "model_validation" = list(
            "gene_mutation_frequency" = gene_mut_freq,
            "event_mutation_frequency" = event_mut_freq
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
        plt = ggpar(plt, font.title=8, font.subtitle=8, font.caption=8, 
                    font.x=8, font.y=8, font.legend=6,
                    font.tickslab=6, font.family=FONT_FAMILY, device=cairo_pdf)
    }
    filename = file.path(directory,paste0(plt_name,extension))
    save_plot(filename, plt, base_width=width, base_height=height, dpi=dpi, units="cm")
}


save_plots = function(plts, figs_dir){
    # model selection
    ## controls
    save_plt(plts, "model_sel-deps_sorted_vs_std_ctl_neg", ".pdf", figs_dir, width=5, height=5)
    save_plt(plts, "model_sel-deps_med_vs_std_ctl_neg", ".pdf", figs_dir, width=5, height=5)
    save_plt(plts, "model_sel-deps_sorted_vs_std_ctl_pos", ".pdf", figs_dir, width=5, height=5)
    save_plt(plts, "model_sel-deps_med_vs_std_ctl_pos", ".pdf", figs_dir, width=5, height=5)
    ## p-value selection
    save_plt(plts, "model_sel-lr_pvalue_ctl_pos", ".pdf", figs_dir, width=4, height=4)
    save_plt(plts, "model_sel-lr_pvalue", ".pdf", figs_dir, width=4, height=4)
    save_plt(plts, "model_sel-tpr_vs_fpr", ".pdf", figs_dir, width=5, height=5)
    ## pearson selection
    save_plt(plts, "model_sel-pearson_corr", ".pdf", figs_dir, width=5, height=5)
    save_plt(plts, "model_sel-pearson_corr_vs_spearman", ".pdf", figs_dir, width=5, height=5)
    
    # model properties
    ## std
    save_plt(plts, "model_prop-selected_vs_event_std", ".pdf", figs_dir, width=4, height=5)
    ## lengths
    save_plt(plts, "model_prop-selected_vs_event_length", ".pdf", figs_dir, width=4, height=5)
    ## n. exons
    save_plt(plts, "model_prop-n_exons_gene", ".pdf", figs_dir, width=4, height=5)
    ## protein impact
    save_plt(plts, "model_prop-protein_impact-counts", ".pdf", figs_dir, width=8, height=8)
    save_plt(plts, "model_prop-protein_impact-freqs", ".pdf", figs_dir, width=8, height=8)
    save_plt(plts, "model_prop-protein_impact-rel_freqs", ".pdf", figs_dir, width=8, height=8)
    save_plt(plts, "model_prop-protein_impact_clean-counts", ".pdf", figs_dir, width=7, height=6)
    save_plt(plts, "model_prop-protein_impact_clean-freqs", ".pdf", figs_dir, width=7, height=6)
    save_plt(plts, "model_prop-protein_impact_clean-rel_freqs", ".pdf", figs_dir, width=7, height=6)
    save_plt(plts, "model_prop-protein_impact_clean_vs_sign_coef-freqs", ".pdf", figs_dir, width=7, height=6)
    ## GSEA
    save_plt(plts, "model_prop-enrichment-hallmarks-dotplot", ".pdf", figs_dir, width=5, height=5)
    save_plt(plts, "model_prop-enrichment-hallmarks-cnetplot", ".pdf", figs_dir, width=8, height=8)
    save_plt(plts, "model_prop-enrichment-GO_BP-dotplot", ".pdf", figs_dir, width=9, height=6)
    save_plt(plts, "model_prop-enrichment-GO_BP-cnetplot", ".pdf", figs_dir, width=20, height=20)
    ## transcriptomic indices
    save_plt(plts, "model_prop-indices-violin", ".pdf", figs_dir, width=4, height=6)
    save_plt(plts, "model_prop-indices-enrichment-GO_BP-dotplot", ".pdf", figs_dir, width=8.2, height=6)
    save_plt(plts, "model_prop-indices-top_pos", ".pdf", figs_dir, width=5.5, height=5)
    save_plt(plts, "model_prop-indices-top_neg", ".pdf", figs_dir, width=5.5, height=5)
    ## tumorigenesis
    save_plt(plts, "model_prop-tumorigenesis-scatter", ".pdf", figs_dir, width=6, height=6)
    save_plt(plts, "model_prop-tumorigenesis_vs_indices", ".pdf", figs_dir, width=5, height=5)
    save_plt(plts, "model_prop-tumorigenesis_vs_prot_imp", ".pdf", figs_dir, width=5, height=5)
    save_plt(plts, "model_prop-harm_scores-scatter", ".pdf", figs_dir, width=8, height=7)
    save_plt(plts, "model_prop-harm_scores_vs_prot_imp", ".pdf", figs_dir, width=8, height=7)
    # ppi closeness
    save_plt(plts, "model_prop-ppi_closeness", ".pdf", figs_dir, width=5, height=6)

    # model validation
    save_plt(plts, "model_val-mutation_gene_count", ".pdf", figs_dir, width=8, height=8)
    save_plt(plts, "model_val-mutation_event_count", ".pdf", figs_dir, width=8, height=10)
    save_plt(plts, "model_val-mutation_gene_frequency", ".pdf", figs_dir, width=8, height=8)
    save_plt(plts, "model_val-mutation_gene_frequency_silent_norm", ".pdf", figs_dir, width=8, height=8)
    #save_plt(plts, "model_val-mutation_gene_frequency_random_norm", ".pdf", figs_dir, width=8, height=8)
    save_plt(plts, "model_val-mutation_event_frequency", ".pdf", figs_dir, width=8, height=8)
    save_plt(plts, "model_val-mutation_event_frequency-by_protein_impact", ".pdf", figs_dir, width=14, height=8)
    #save_plt(plts, "model_val-mutation_gene_frequency_vs_random", ".pdf", figs_dir, width=8, height=8)
    #save_plt(plts, "model_val-mutation_event_frequency_vs_random", ".pdf", figs_dir, width=8, height=8)
    
    # event mutation distances
    save_plt(plts, "mutation_dists-closest_ss-distrs", ".pdf", figs_dir, width=12, height=14)
    save_plt(plts, "mutation_dists-5ss-distrs", ".pdf", figs_dir, width=12, height=14)
    save_plt(plts, "mutation_dists-3ss-distrs", ".pdf", figs_dir, width=12, height=14)
    
    # events oi
    genes_oi = c("KRAS","SMNDC1","NUMB","NUMB_lung")
    lapply(genes_oi, function(gene_oi){
        lapply(c("spldep_vs_rnai","genexpr_vs_rnai","genexpr_vs_splicing","splicing_vs_rnai"), 
               function(comparison){
            plt_name = sprintf("events_oi-%s-%s",gene_oi,comparison)
            save_plt(plts, plt_name, ".pdf", figs_dir, width=4, height=4)            
        })
    })
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

parseargs = function(){
    
    option_list = list( 
        make_option("--models_file", type="character"),
        make_option("--ccle_stats_file", type="character"),
        make_option("--indices_file", type="character"),
        make_option("--event_mut_file", type="character"),
        make_option("--gene_mut_freq_file", type="character"),
        make_option("--event_mut_freq_file", type="character"),
        make_option("--msigdb_dir", type="character"),
        make_option("--spldep_file", type="character"),
        make_option("--rnai_file", type="character"),
        make_option("--genexpr_file", type="character"),
        make_option("--splicing_file", type="character"),
        make_option("--cancer_events_file", type="character"),
        make_option("--ascanceratlas_file", type="character"),
        make_option("--event_info_file", type="character"),
        make_option("--metadata_file", type="character"),
        make_option("--ppi_closeness_file", type="character"),
        make_option("--randsel_events_file", type="character"),
        make_option("--randsel_genes_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}


main = function(){
    args = parseargs()
    
    models_file = args[["models_file"]]
    ccle_stats_file = args[["ccle_stats_file"]]
    indices_file = args[["indices_file"]]
    event_mut_file = args[["event_mut_file"]]
    gene_mut_freq_file = args[["gene_mut_freq_file"]]
    event_mut_freq_file = args[["event_mut_freq_file"]]
    msigdb_dir = args[["msigdb_dir"]]
    protein_impact_file = args[["protein_impact_file"]]
    cosmic_genes_file = args[["cosmic_genes_file"]]
    spldep_file = args[["spldep_file"]]
    rnai_file = args[["rnai_file"]]
    genexpr_file = args[["genexpr_file"]]
    splicing_file = args[["splicing_file"]]
    cancer_events_file = args[["cancer_events_file"]]
    ascanceratlas_file = args[["ascanceratlas_file"]]
    event_info_file = args[["event_info_file"]]
    metadata_file = args[["metadata_file"]]
    ppi_closeness_file = args[["ppi_closeness_file"]]
    randsel_events_file = args[["randsel_events_file"]]
    randsel_genes_file = args[["randsel_genes_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    models = read_tsv(models_file)
    ccle_stats = read_tsv(ccle_stats_file)
    indices = read_tsv(indices_file) %>% filter(index_name %in% c("stemness"))
    event_mut = read_tsv(event_mut_file)
    gene_mut_freq = read_tsv(gene_mut_freq_file) %>% dplyr::rename(GENE=Hugo_Symbol)
    event_mut_freq = read_tsv(event_mut_freq_file)
    ontologies = load_ontologies(msigdb_dir, protein_impact_file, cosmic_genes_file)
    spldep = read_tsv(spldep_file)
    rnai = read_tsv(rnai_file)
    genexpr = read_tsv(genexpr_file)
    splicing = read_tsv(splicing_file)
    cancer_events = read_tsv(cancer_events_file)
    ascanceratlas = read_tsv(ascanceratlas_file)
    event_info = read_tsv(event_info_file)
    metadata = read_tsv(metadata_file)
    ppi_closeness = read_tsv(ppi_closeness_file) %>% 
        mutate(type = gsub("_.*","",dataset_id))
    
    randsel_events = read_tsv(randsel_events_file)
    randsel_genes = read_tsv(randsel_genes_file)
    
    gc()
    
    # prep cancer events
    ascanceratlas = ascanceratlas %>%
        dplyr::rename(EVENT = EVENT_perf) %>%
        left_join(
            event_info %>% distinct(EVENT,GENE),
            by=c("EVENT")
        ) %>%
        drop_na(EVENT) %>%
        mutate(source="ASCancerAtlas")
    
    cols_oi = c('EVENT','ENSEMBL','GENE','length','source_authors','source_year','source_doi','description','event_id','source')
    cancer_events = cancer_events %>%
        mutate(source="HandCurated") %>%
        bind_rows(ascanceratlas) %>%
        dplyr::select(cols_oi) %>%
        distinct()
    
    # log normalize gene expression
    genexpr = genexpr %>% mutate_at(vars(-("ID")), function(x){ log2(x+1) })
    
    # add info to models
    models = models %>% 
        # filter(n_obs > MIN_N_OBS) %>%
        mutate(event_gene = paste0(EVENT,"_",GENE),
               event_type = gsub("Hsa","",gsub("[^a-zA-Z]", "",EVENT)),
               is_selected = pearson_correlation_mean > THRESH_CORR & 
                             lr_pvalue < THRESH_LR_PVALUE,
               is_selected = replace_na(is_selected, FALSE)) %>% 
        dplyr::select(-c(event_mean,event_std,gene_mean,gene_std)) %>% 
        left_join(ccle_stats, by=c("EVENT","ENSEMBL","GENE")) %>%
        left_join(event_info %>% distinct(EVENT,LE_o), by="EVENT") %>%
        left_join(ontologies[["protein_impact"]], by="EVENT")
    
    selected_events = models %>% filter(is_selected) %>% pull(EVENT)
    
    # prep for plotting
    ## model selection
    rnai_stats = get_rnai_stats(rnai, models)
    spldep_stats = get_spldep_stats(spldep, models)
    splicing_stats = get_spldep_stats(splicing %>% rename(index=EVENT), models)
    
    ### get controls
    ctl_pos_genes = cancer_events %>% pull(GENE) %>% unique()
    ctl_pos = rnai_stats %>% 
        filter(GENE %in% ctl_pos_genes) %>% 
        slice_max(order_by=rnai_std, n=SIZE_CTL) %>% 
        left_join(cancer_events, by="GENE") %>%
        pull(EVENT) %>%
        unique()
    ctl_neg = rnai_stats %>% 
        filter(!(GENE %in% ctl_pos_genes)) %>% # do not include positive controls
        slice_min(order_by = rnai_std, n=SIZE_CTL) %>% # selected top 100
        pull(GENE)
    cancer_events = cancer_events %>%
        mutate(is_ctl_pos = EVENT %in% ctl_pos,
               is_ctl_neg = GENE %in% ctl_neg)
    ### evaluate thresholds
    eval_pvalue = thresh_eval_pvalue(models, ctl_neg, ctl_pos)
    eval_corr = thresh_eval_corr(models, rnai, spldep)
    
    ## model properties
    ### GSEA of selected events
    genes_oi = models %>% 
            filter(is_selected) %>% 
            pull(GENE) %>% 
            unique()
    events_oi = models %>% 
            filter(is_selected) %>% 
            pull(EVENT)
    universe = list(
        "genes" = models %>% pull(GENE) %>% unique(),
        "events" = models %>% pull(EVENT) %>% unique()
    )
    enrichment = run_enrichment(genes_oi, events_oi, universe, ontologies)
    enrichment[sapply(enrichment, nrow)<1] = NULL
    
    ## GSEA of selected events correlating with transcriptomic indices
    indices_enrich = run_enrichments_indices(indices, selected_events, event_info)
    
    ## mutations affecting selected events
    event_mut = event_mut %>%
        mutate(is_selected = EVENT %in% selected_events)
    
    # compute harm score
    ## H.S. = (-1) * SplDep * DeltaPSI
    ## where DeltaPSI is in the direction of inclusion if SplDep>0, or exclusion if SplDep<0
    common_samples = intersect(colnames(spldep), colnames(splicing))
    common_events = intersect(spldep[["index"]], splicing[["EVENT"]])
    common_events = intersect(common_events, selected_events)
       
    spldep_mat = spldep %>% column_to_rownames("index")
    spldep_mat = spldep_mat[common_events, common_samples]
    splicing_mat = splicing %>% column_to_rownames("EVENT")
    splicing_mat = splicing_mat[common_events, common_samples]
    
    psi_final = spldep_mat
    psi_final[psi_final<0] = 0 # remove oncoexons
    psi_final[psi_final>0] = 100 # include tumor-suppressor exons
    harm = (-1) * spldep_mat * (psi_final - splicing_mat)
    
    harm_stats = get_spldep_stats(harm %>% rownames_to_column("index"), models)
    
    # number of selected genes in COSMIC
    available_genes = models %>% pull(GENE) %>% unique() 
    selected_genes = models %>% filter(is_selected) %>% pull(GENE) %>% unique()
    cosmic_test = data.frame(gene = available_genes) %>%
        mutate(
            has_cancer_driver_exon = gene %in% selected_genes,
            is_cancer_driver = gene %in% ontologies[["cosmic"]][["gene"]]
        ) %>%
        with(table(has_cancer_driver_exon, is_cancer_driver)) %>%
        fisher.test()
    # 74 cancer driver genes have cancer-driver exons
    
    # plot
    plts = make_plots(
        models, rnai_stats, cancer_events, 
        eval_pvalue, eval_corr, 
        enrichment, indices, indices_enrich, spldep_stats, harm_stats, ppi_closeness,
        gene_mut_freq, event_mut_freq, randsel_genes, randsel_events,
        rnai, spldep, splicing, genexpr, metadata
    )

    # make figdata
    figdata = make_figdata(
        models, rnai_stats, cancer_events, 
        eval_pvalue, eval_corr, 
        enrichment, indices, indices_enrich, spldep_stats, harm_stats, ppi_closeness,
        gene_mut_freq, event_mut_freq, randsel_genes, randsel_events,
        rnai, spldep, splicing, genexpr, metadata
    )
    
    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}