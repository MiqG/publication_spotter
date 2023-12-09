#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# EDA of gene dependencies regressed on event PSI and gene TPMs.
# 

require(optparse)
require(ensembldb)
require(tidyverse)
require(ggpubr)
require(cowplot)
require(scattermore)
require(ggrepel)
require(clusterProfiler)
require(ComplexHeatmap)
require(ggplotify)
require(grid)
require(extrafont)
require(gtools)
require(tidytext)
require(ggvenn)
require(pROC)

# variables
#MIN_N_OBS = 20
THRESH_LR_PVALUE = 0.025
THRESH_CORR = 0.2
SIZE_CTL = 100
RANDOM_SEED = 1234

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

PAL_COSMIC_COMP = setNames(c("#d5b1d8ff","orange"),c("is_cancer_driver","has_cancer_driver_exon"))
PAL_DRIVER_COMP = c(
    "in_cosmic" = "#d5b1d8ff",
    "is_known_driver" = PAL_SINGLE_ACCENT,
    "is_pan_essential" = PAL_SINGLE_NEUTRAL
)

# Development
# -----------
# ROOT = here::here()
# PREP_DIR = file.path(ROOT,"data","prep")
# RAW_DIR = file.path(ROOT,"data","raw")
# RESULTS_DIR = file.path(ROOT,"results","model_splicing_dependency")
# SUPPORT_DIR = file.path(ROOT,"support")
# models_file = file.path(RESULTS_DIR,"files","models_gene_dependency-EX","model_summaries.tsv.gz")
# ccle_stats_file = file.path(PREP_DIR,"stats","CCLE.tsv.gz")
# gene_mut_freq_file = file.path(ROOT,"data","prep","gene_mutation_freq","CCLE.tsv.gz")
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
# #(TODO: preprocessing)shrna_raw_file = file.path(RAW_DIR,"DepMap","demeter2","")
# shrna_seqs_file = file.path(RAW_DIR,"DepMap","demeter2","shRNA-mapping.csv")
# shrna_mapping_ensembl_file = file.path(PREP_DIR,"demeter2","shRNA_to_gencode.v44.transcripts-mapping_coords.tsv.gz")
# shrna_mapping_vastdb_file = file.path(PREP_DIR,"demeter2","shRNA-mapping_to_vastdb_exons.tsv.gz")
# models_achilles_file = file.path(RESULTS_DIR,'files','achilles','models_gene_dependency-EX','model_summaries.tsv.gz')
# crispr_essentials_file = file.path(SUPPORT_DIR,"CRISPRInferredCommonEssentials.csv")
# ccle_info_file = file.path(SUPPORT_DIR,"ENA_filereport-PRJNA523380-CCLE.tsv")
# gencode_annot_file = file.path(RAW_DIR,"GENCODE","gencode.v44.annotation.gtf.gz")
# genome_annot_file = file.path(RAW_DIR,"ENSEMBL","Homo_sapiens.Gh38.110.sqlite")
# figs_dir = file.path(RESULTS_DIR,"figures","model_selection")

##### FUNCTIONS #####
load_ontologies = function(msigdb_dir, protein_impact_file, cosmic_genes_file){
    ontologies = list(
        "reactome" = read.gmt(file.path(msigdb_dir,"c2.cp.reactome.v7.4.symbols.gmt")),
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
    enrichments[["cosmic"]] = enricher(genes, TERM2GENE=ontologies[["cosmic"]], universe=universe[["genes"]])
    enrichments[["reactome"]] = enricher(genes, TERM2GENE=ontologies[["reactome"]], universe=universe[["genes"]])
    
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


thresh_eval_corr = function(models, rnai, spldep, ctl_neg, ctl_pos){
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


compute_auc = function(x, y){
    x[is.na(x)] = 0
    y[is.na(y)] = 0
    
    x = sort(x)
    y = sort(y)
    
    diffs.x = x[-1] - x[-length(x)]
    means.vert = (y[-1] + y[-length(y)])/2
    auc = sum(means.vert * diffs.x)

    return(auc)
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
    
    set.seed(RANDOM_SEED)
    X = models %>%
        distinct(EVENT, GENE, lr_pvalue, is_known_driver) %>%
        pivot_longer(-c(EVENT, GENE, lr_pvalue), names_to="class_type", values_to="in_class") %>%
        bind_rows(
            models %>%
            distinct(GENE, lr_pvalue, is_pan_essential, in_cosmic) %>%
            pivot_longer(-c(GENE,lr_pvalue), names_to="class_type", values_to="in_class") %>%
            group_by(GENE, class_type) %>%
            slice_sample(n=1) %>%
            ungroup()
        ) %>%
        mutate(class_type=factor(class_type, levels=c("is_known_driver","is_pan_essential","in_cosmic")))
    counts = X %>% group_by(class_type) %>% count(in_class) %>% mutate(label=sprintf("n=%s",n))
    plts[["model_sel-gt_lists_vs_lr_pvalue-random"]] = X %>%
        ggviolin(x="in_class", y="lr_pvalue", fill="in_class", color=NA, palette=PAL_DUAL, trim=TRUE) +
        geom_boxplot(width=0.1, outlier.size=0.1) +
        geom_text(aes(label=label, y=-0.05), counts, family=FONT_FAMILY, size=FONT_SIZE)+
        stat_compare_means(method="wilcox.test", family=FONT_FAMILY, size=FONT_SIZE) +
        facet_wrap(~class_type, scales="free") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(fill="none") + 
        labs(x="", y="LR Test p-value")

    set.seed(RANDOM_SEED)
    X = models %>%
        distinct(EVENT, GENE, lr_pvalue, is_known_driver) %>%
        pivot_longer(-c(EVENT, GENE, lr_pvalue), names_to="class_type", values_to="in_class") %>%
        bind_rows(
            models %>%
            distinct(GENE, lr_pvalue, is_pan_essential, in_cosmic) %>%
            pivot_longer(-c(GENE,lr_pvalue), names_to="class_type", values_to="in_class") %>%
            group_by(GENE, class_type) %>%
            slice_min(lr_pvalue, n=1) %>%
            ungroup()
        ) %>%
        mutate(class_type=factor(class_type, levels=c("is_known_driver","is_pan_essential","in_cosmic")))
    counts = X %>% group_by(class_type) %>% count(in_class) %>% mutate(label=sprintf("n=%s",n))
    plts[["model_sel-gt_lists_vs_lr_pvalue-lowest"]] = X %>%
        ggviolin(x="in_class", y="lr_pvalue", fill="in_class", color=NA, palette=PAL_DUAL, trim=TRUE) +
        geom_boxplot(width=0.1, outlier.size=0.1) +
        geom_text(aes(label=label, y=-0.05), counts, family=FONT_FAMILY, size=FONT_SIZE)+
        stat_compare_means(method="wilcox.test", family=FONT_FAMILY, size=FONT_SIZE) +
        facet_wrap(~class_type, scales="free") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(fill="none") + 
        labs(x="", y="LR Test p-value")
    
    ## ROC analysis - LR p-value
    ### to detect cancer-driver exons
    roc_cancer_driver = models %>%
        mutate(is_ctl_pos = is_known_driver) %>%
        roc(response=is_ctl_pos, predictor=lr_pvalue, 
            levels=c(TRUE,FALSE), direction="<") %>%
        coords(transpose=FALSE, re="all") %>% 
        mutate(
            fpr = 1 - specificity,
            auc_tf = compute_auc(fpr, sensitivity),
            auc_pr = compute_auc(recall, precision),
            ground_truth = "is_known_driver"
        )
    
    ### to detect common essentials
    roc_essentials = models %>%
        distinct(lr_pvalue, is_pan_essential, GENE) %>%
        group_by(GENE) %>%
        slice_min(lr_pvalue, n=1) %>%
        ungroup() %>%
        roc(response=is_pan_essential, predictor=lr_pvalue, 
            levels=c(TRUE,FALSE), direction="<") %>%
        coords(transpose=FALSE, re="all") %>% 
        mutate(
            fpr = 1 - specificity,
            auc_tf = compute_auc(fpr, sensitivity),
            auc_pr = compute_auc(recall, precision),
            ground_truth = "is_pan_essential"
        )
    
    ### to detect COSMIC
    roc_cosmic = models %>%
        distinct(lr_pvalue, in_cosmic, GENE) %>%
        group_by(GENE) %>%
        slice_min(lr_pvalue, n=1) %>%
        ungroup() %>%
        roc(response=in_cosmic, predictor=lr_pvalue, 
            levels=c(TRUE,FALSE), direction="<") %>%
        coords(transpose=FALSE, re="all") %>% 
        mutate(
            fpr = 1 - specificity,
            auc_tf = compute_auc(fpr, sensitivity),
            auc_pr = compute_auc(recall, precision),
            ground_truth = "in_cosmic"
        )
    
    roc_analysis = roc_cancer_driver %>%
        bind_rows(roc_essentials) %>%
        bind_rows(roc_cosmic)
    
    plts[["model_sel-roc_analysis-fpr_vs_tpr-lr_pvalue-line"]] = roc_analysis %>%
        mutate(auc_lab = sprintf("AUC=%s",round(auc_tf,2))) %>%
        ggplot(aes(x=fpr, y=sensitivity, color=ground_truth), size=1) +
        geom_line() +
        geom_abline(intercept = 0, slope = 1, linetype="dashed", size=LINE_SIZE) +
        color_palette(PAL_DRIVER_COMP) +
        theme_pubr() +
        theme(aspect.ratio=1) +
        geom_text(
            aes(x=x, y=y, 
            label=auc_lab),
            . %>% distinct(ground_truth, auc_lab) %>%
                mutate(x=0.75, y=c(0.12, 0.06, 0)),
            hjust=0, size=FONT_SIZE, family=FONT_FAMILY
        ) +
        labs(
            x="False Positive Rate (1 - Specificity)", 
            y="True Positive Rate (Sensitivity)"
        )
    
    plts[["model_sel-roc_analysis-recall_vs_precision-lr_pvalue-line"]] = roc_analysis %>%
        mutate(
            auc_lab = sprintf("AUC=%s",round(auc_pr,2)),
            recall = replace_na(recall, 0),
            precision = replace_na(precision, 0),
        ) %>%
        ggplot(aes(x=recall, y=precision, color=ground_truth), size=1) +
        geom_line() +
        color_palette(PAL_DRIVER_COMP) +
        theme_pubr() +
        theme(aspect.ratio=1) +
        geom_text(
            aes(x=x, y=y, 
            label=auc_lab),
            . %>% distinct(ground_truth, auc_lab) %>%
                mutate(x=0.75, y=c(0.8, 0.83, 0.86)),
            hjust=0, size=FONT_SIZE, family=FONT_FAMILY
        ) +
        labs(
            x="Recall", 
            y="Precision"
        )
    
    ## ROC analysis - Pearson
    ### to detect cancer-driver exons
    roc_cancer_driver = models %>%
        filter(lr_pvalue < THRESH_LR_PVALUE) %>%
        mutate(is_ctl_pos = is_known_driver) %>%
        roc(response=is_ctl_pos, predictor=pearson_correlation_mean, 
            levels=c(TRUE,FALSE), direction=">") %>%
        coords(transpose=FALSE, re="all") %>% 
        mutate(
            fpr = 1 - specificity,
            auc_tf = compute_auc(fpr, sensitivity),
            auc_pr = compute_auc(recall, precision),
            ground_truth = "is_known_driver"
        )
    
    ### to detect common essentials
    roc_essentials = models %>%
        filter(lr_pvalue < THRESH_LR_PVALUE) %>%
        distinct(pearson_correlation_mean, is_pan_essential, GENE) %>%
        group_by(GENE) %>%
        slice_max(pearson_correlation_mean, n=1) %>%
        ungroup() %>%
        roc(response=is_pan_essential, predictor=pearson_correlation_mean, 
            levels=c(TRUE,FALSE), direction=">") %>%
        coords(transpose=FALSE, re="all") %>% 
        mutate(
            fpr = 1 - specificity,
            auc_tf = compute_auc(fpr, sensitivity),
            auc_pr = compute_auc(recall, precision),
            ground_truth = "is_pan_essential"
        )
    
    ### to detect COSMIC
    roc_cosmic = models %>%
        filter(lr_pvalue < THRESH_LR_PVALUE) %>%
        distinct(pearson_correlation_mean, in_cosmic, GENE) %>%
        group_by(GENE) %>%
        slice_max(pearson_correlation_mean, n=1) %>%
        ungroup() %>%
        roc(response=in_cosmic, predictor=pearson_correlation_mean, 
            levels=c(TRUE,FALSE), direction=">") %>%
        coords(transpose=FALSE, re="all") %>% 
        mutate(
            fpr = 1 - specificity,
            auc_tf = compute_auc(fpr, sensitivity),
            auc_pr = compute_auc(recall, precision),
            ground_truth = "in_cosmic"
        )
    
    roc_analysis = roc_cancer_driver %>%
        bind_rows(roc_essentials) %>%
        bind_rows(roc_cosmic)
    
    plts[["model_sel-roc_analysis-fpr_vs_tpr-pearson-line"]] = roc_analysis %>%
        mutate(auc_lab = sprintf("AUC=%s",round(auc_tf,2))) %>%
        ggplot(aes(x=fpr, y=sensitivity, color=ground_truth), size=1) +
        geom_line() +
        geom_abline(intercept = 0, slope = 1, linetype="dashed", size=LINE_SIZE) +
        color_palette(PAL_DRIVER_COMP) +
        theme_pubr() +
        theme(aspect.ratio=1) +
        geom_text(
            aes(x=x, y=y, 
            label=auc_lab),
            . %>% distinct(ground_truth, auc_lab) %>%
                mutate(x=0.75, y=c(0.12, 0.06, 0)),
            hjust=0, size=FONT_SIZE, family=FONT_FAMILY
        ) +
        labs(
            x="False Positive Rate (1 - Specificity)", 
            y="True Positive Rate (Sensitivity)"
        )
    
    plts[["model_sel-roc_analysis-recall_vs_precision-pearson-line"]] = roc_analysis %>%
        mutate(
            auc_lab = sprintf("AUC=%s",round(auc_pr,2)),
            recall = replace_na(recall, 0),
            precision = replace_na(precision, 0),
        ) %>%
        ggplot(aes(x=recall, y=precision, color=ground_truth), size=1) +
        geom_line() +
        color_palette(PAL_DRIVER_COMP) +
        theme_pubr() +
        theme(aspect.ratio=1) +
        geom_text(
            aes(x=x, y=y, 
            label=auc_lab),
            . %>% distinct(ground_truth, auc_lab) %>%
                mutate(x=0.75, y=c(0.8, 0.83, 0.86)),
            hjust=0, size=FONT_SIZE, family=FONT_FAMILY
        ) +
        labs(
            x="Recall", 
            y="Precision"
        )    
    
    ### overlap ctl pos and common essentials
    plts[["model_sel-set_overlaps-all-venn"]] = models %>%
        group_by(GENE) %>%
        mutate(is_known_driver_gene = any(is_known_driver)) %>%
        ungroup() %>%
        distinct(GENE, is_known_driver_gene, is_pan_essential, in_cosmic) %>%
        ggplot(aes(A=in_cosmic, B=is_known_driver_gene, C=is_pan_essential)) +
        geom_venn(stroke_color=NA, 
                  fill_color=PAL_DRIVER_COMP,
                  set_name_size = FONT_SIZE+0.5, text_size = FONT_SIZE) +
        coord_fixed() +
        theme_void()
    
    plts[["model_sel-set_overlaps-selected-venn"]] = models %>%
        group_by(GENE) %>%
        mutate(is_known_driver_gene = any(is_known_driver)) %>%
        ungroup() %>%
        filter(is_selected) %>%
        distinct(GENE, is_known_driver_gene, is_pan_essential, in_cosmic) %>%
        ggplot(aes(A=in_cosmic, B=is_known_driver_gene, C=is_pan_essential)) +
        geom_venn(stroke_color=NA, 
                  fill_color=PAL_DRIVER_COMP,
                  set_name_size = FONT_SIZE+0.5, text_size = FONT_SIZE) +
        coord_fixed() +
        theme_void()
    
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
    
    # LR p-values and Pearson correlations are independent of number of observations
    plts[["model_sel-n_obs_vs_pvalue-scatter"]] = X %>%
        ggplot(aes(x=n_obs, y=lr_pvalue)) +
        geom_scattermore(pixels = c(1000,1000), pointsize=5, alpha=0.5, color=PAL_SINGLE_NEUTRAL) +
        geom_hline(yintercept=THRESH_LR_PVALUE, linetype="dashed", color="black", size=LINE_SIZE, alpha=0.2 ) +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        theme_pubr() +
        theme(aspect.ratio=1) +
        labs(x="N Observations", y="LR p-value")
    
    plts[["model_sel-n_obs_vs_pearson-scatter"]] = X %>%
        filter(lr_pvalue < THRESH_LR_PVALUE) %>%
        ggplot(aes(x=n_obs, y=pearson_correlation_mean)) +
        geom_scattermore(pixels = c(1000,1000), pointsize=5, alpha=0.5, color=PAL_SINGLE_NEUTRAL) +
        geom_hline(yintercept=THRESH_CORR, linetype="dashed", color="black", size=LINE_SIZE, alpha=0.2 ) +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        theme_pubr() +
        theme(aspect.ratio=1) +
        labs(x="N Observations", y="Pearson Correlation Mean")
    
    return(plts)
}


plot_model_properties = function(models, enrichment,
                                 spldep_stats, harm_stats){
    
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
    
    return(plts)
}


plot_model_validation = function(models, gene_mut_freq){
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


plot_cosmic_comparison = function(cosmic_comparison){
    plts = list()
    
    X = cosmic_comparison
    
    x = sapply(c("is_cancer_driver", "has_cancer_driver_exon"), function(col){
            genes_oi = X %>%
                distinct(gene, is_cancer_driver, has_cancer_driver_exon) %>%
                filter(get(col)) %>%
                pull(gene)
            return(genes_oi)
    }, simplify=FALSE)
    plts[["cosmic_comparison-overlap-venn"]] = x %>%
        ggvenn(
            fill_color = as.character(PAL_COSMIC_COMP),
            stroke_color = NA,
            set_name_size = FONT_SIZE+0.5,
            text_size = FONT_SIZE
        )
    
    
    x = X %>%
        pivot_longer(
            c(has_cancer_driver_exon, is_cancer_driver), 
            names_to="gene_set", values_to="in_gene_set"
        ) %>%
        filter(in_gene_set) %>%
        drop_na(term) %>%
        group_by(gene_set) %>%
        count(term) %>%
        mutate(perc = n / sum(n)) %>%
        ungroup()
    terms_oi = x %>%
        group_by(gene_set) %>%
        slice_max(perc, n=10) %>%
        ungroup() %>%
        pull(term) %>%
        unique()
    terms_order = x %>%
        filter(term %in% terms_oi) %>%
        pivot_wider(id_cols=term, names_from=gene_set, values_from=perc) %>%
        mutate(perc_diff = has_cancer_driver_exon - is_cancer_driver) %>%
        arrange(-has_cancer_driver_exon) %>%
        pull(term)
        
    plts[["cosmic_comparison-reactome-bar"]] = x %>%
        filter(term %in% terms_oi) %>%
        mutate(term = factor(term, levels=terms_order)) %>%
        ggplot(aes(x=term, y=perc, group=gene_set)) +
        geom_col(aes(fill=gene_set), color=NA, position="dodge") + 
        geom_text(
            aes(label=n), position=position_dodge(width=0.9),
            size=FONT_SIZE, family=FONT_FAMILY, hjust=0, vjust=0, angle=45
        ) +
        fill_palette(PAL_COSMIC_COMP) +
        theme_pubr(x.text.angle = 70) +
        labs(x="Reactome Pathway", y="% Overlap", fill="Gene Set")
    
    return(plts)
}


plot_rnai_qc = function(shrna_mapping_ensembl, event_info, models){
    plts = list()

    # add VastDB annotations to mapped shRNAs
    # event_coords = event_info %>%
    #     filter(str_detect(EVENT,"EX")) %>%
    #     mutate(
    #         seqnames = gsub("chr","",gsub(":.*","",COORD_o)),
    #         start_exon = gsub("-.*","",gsub(".*:","",COORD_o)) %>% as.integer(),
    #         end_exon = gsub(".*-","",gsub(".*:","",COORD_o)) %>% as.integer(),
    #         strand = gsub(".*:","",REF_CO)
    #     ) %>%
    #     distinct(EVENT, seqnames, start_exon, end_exon, strand, GENE) %>%
    #     drop_na()
    
    shrna_mapping = shrna_seqs %>%
        left_join(
            shrna_mapping_vastdb, 
            by=c("Barcode Sequence"="shrna_barcode", "Gene Symbol"="GENE")
        ) %>%
        dplyr::rename(barcode_sequence=`Barcode Sequence`, gene_name_prevmap=`Gene Symbol`) %>%
        distinct(barcode_sequence, gene_name_prevmap, EVENT) %>%
        left_join(
            shrna_mapping_ensembl %>%
                mutate(
                    exon_id_coords = sprintf("%s:%s-%s:%s",seqnames,start_exon,end_exon,strand),
                    gene_id_clean = gsub("\\..*","",gene_id)
                ) %>% distinct(barcode_sequence, gene_name_prevmap, exon_id_coords),
            by=c("barcode_sequence", "gene_name_prevmap")
        ) %>%
        mutate(
            in_models = gene_name_prevmap %in% models[["GENE"]]
        ) %>%
        left_join(models %>% distinct(EVENT, GENE, ENSEMBL, event_median, event_std, lr_pvalue), by="EVENT")
    
    # how many different exons are targeted by shRNAs in each gene
    n_targeted_exons_per_gene = shrna_mapping %>%
        distinct(gene_name_prevmap, exon_id_coords) %>%
        drop_na() %>%
        count(gene_name_prevmap, name="n_targeted_exons_per_gene") %>%
        mutate(mapping_type="ensembl") %>%
        bind_rows(
            shrna_mapping %>%
            distinct(gene_name_prevmap, EVENT) %>%
            drop_na() %>%
            count(gene_name_prevmap, name="n_targeted_exons_per_gene") %>%
            mutate(mapping_type="vastdb")
        ) %>%
        group_by(gene_name_prevmap) %>%
        slice_max(n_targeted_exons_per_gene, n=1) %>%
        ungroup() %>%
        distinct(gene_name_prevmap, n_targeted_exons_per_gene)
    
    shrna_mapping = shrna_mapping %>%
        left_join(n_targeted_exons_per_gene, by="gene_name_prevmap")
    
    # how many exons are targeted for each gene?
    plts[["rnai_qc-n_targeted_exons_per_gene-hist"]] = shrna_mapping %>%
        # filter out those shRNAs mapped to the wrong genes
        # all genes in models should have been mapped!!
        filter(in_models) %>%
        distinct(gene_name_prevmap, n_targeted_exons_per_gene) %>%
        count(n_targeted_exons_per_gene) %>%
        arrange(n_targeted_exons_per_gene) %>%
        ggbarplot(
            x="n_targeted_exons_per_gene", y="n", fill=PAL_SINGLE_NEUTRAL, color=NA,
            label=TRUE, lab.family=FONT_FAMILY, lab.size=FONT_SIZE+2
        ) +
        labs(x="N Exons targeted by shRNA per Gene", y="Count")
    
    # are targeted exons constitutive?
    singletons = shrna_mapping %>% 
        filter(in_models) %>%
        distinct(EVENT, gene_name_prevmap, n_targeted_exons_per_gene) %>%
        filter(n_targeted_exons_per_gene==1) %>%
        count(EVENT, gene_name_prevmap) %>% 
        drop_na() %>%
        pull(EVENT)
    X = models %>%
        left_join(n_targeted_exons_per_gene, by=c("GENE"="gene_name_prevmap")) %>%
        mutate(
            is_targeted = ifelse(
                EVENT %in% shrna_mapping[["EVENT"]], "Targeted", "Not Targeted"
            ),
            is_singleton = EVENT %in% singletons
        ) %>%
        group_by(is_targeted) %>%
        mutate(is_targeted_lab = sprintf("%s (n=%s)",is_targeted,n())) %>%
        ungroup()
    
    plts[["rnai_qc-exon_psi_median_vs_std_vs_shrna_targeted-scatter"]] = X %>%
        ggplot(aes(x=event_median, y=event_std)) +
        geom_scattermore(aes(color=is_singleton), . %>% filter(!is_singleton), pixels = c(1000,1000), pointsize=5, alpha=0.5) +
        geom_scattermore(aes(color=is_singleton), . %>% filter(is_singleton), pixels = c(1000,1000), pointsize=8, alpha=0.5) +
        geom_vline(xintercept=seq(0,100,10), linetype="dashed", size=LINE_SIZE, color="black") +
        color_palette(c(PAL_SINGLE_NEUTRAL, PAL_SINGLE_LIGHT)) +
        theme_pubr() +
        theme(aspect.ratio=1) +
        facet_wrap(~is_targeted_lab) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Exon PSI Median", y="Exon PSI Std", color="Is Singleton")

    plts[["rnai_qc-exon_psi_median_vs_shrna_targeted-bar"]] = X %>% 
        mutate(event_psi_bins = cut(event_median, breaks=seq(0,100,10), include.lowest=TRUE)) %>%
        count(is_targeted_lab, is_singleton, event_psi_bins) %>%
        ggbarplot(x="event_psi_bins", y="n", color=NA, fill="is_singleton", 
                  palette=c(PAL_SINGLE_NEUTRAL, PAL_SINGLE_LIGHT)) +    
        facet_wrap(~is_targeted_lab+is_singleton, scales="free_y") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(fill="none") +
        labs(x="Exon PSI Median Bin", y="N Exons")
    
    # impact of number of exons targeted per gene on our analysis
    plts[["rnai_qc-n_targeted_exons_per_gene_vs_lr_pvalue-box"]] = X %>%
        ggboxplot(x="n_targeted_exons_per_gene", y="lr_pvalue", fill=NA, color=PAL_SINGLE_NEUTRAL, width=0.5, outlier.size=0.1) +
        labs(x="N Exons Targeted by shRNA per Gene", y="LR p-value")
    
    plts[["rnai_qc-n_targeted_exons_per_gene_vs_abs_event_coef-box"]] = X %>% 
        mutate(abs_coef = abs(event_coefficient_mean)) %>%
        ggboxplot(x="n_targeted_exons_per_gene", y="abs_coef", fill=NA, color=PAL_SINGLE_NEUTRAL, width=0.5, outlier.size=0.1) +
        yscale("log10", .format=TRUE) +
        labs(x="N Exons Targeted by shRNA per Gene", y="abs(Event Coeff.)")

    return(plts)
}


plot_rnaseq_qc = function(metadata, rnai){
    plts = list()

    X = metadata %>%
        filter(DepMap_ID %in% colnames(rnai))
    
    plts[["rnaseq_qc-read_count_vs_n_detected_events-scatter"]] = X %>% 
        ggscatter(x="read_count", y="n_detected_events", alpha=0.5, size=1) + 
        geom_smooth(linetype="dashed", color="black", size=LINE_SIZE, alpha=0.2) +
        theme(aspect.ratio=1) +
        labs(x="Read Count", y="N Different Exons Detected")
    
    return(plts)
}


make_plots = function(models, rnai_stats, cancer_events, 
                      eval_pvalue, eval_corr, 
                      enrichment, spldep_stats, harm_stats, 
                      gene_mut_freq, 
                      rnai, spldep, splicing, genexpr, metadata, 
                      cosmic_comparison, shrna_mapping){
    plts = list(
        plot_model_selection(models, rnai_stats, cancer_events, eval_pvalue, eval_corr),
        plot_model_properties(models, enrichment,spldep_stats, harm_stats),
        plot_model_validation(models, gene_mut_freq),
        plot_events_oi(models, cancer_events, rnai, spldep, splicing, genexpr, metadata),
        plot_cosmic_comparison(cosmic_comparison),
        plot_rnai_qc(models, shrna_mapping, rnai_stats),
        plot_rnaseq_qc(metadata, rnai)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(models, rnai_stats, cancer_events, 
                      eval_pvalue, eval_corr, 
                      enrichment, spldep_stats, harm_stats,
                      gene_mut_freq, 
                      rnai, spldep, splicing, genexpr, metadata, 
                      cosmic_comparison, shrna_mapping){
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
        "qc" = list(
            "shrna_mapping" = shrna_mapping,
        ),
        "model_selection" = list(
            "model_summaries"= models,
            "rnai_stats" = rnai_stats,
            "cancer_events" = cancer_events,
            "evaluation_pvalue" = eval_pvalue,
            "evaluation_correlation" = eval_corr
        ),
        "model_properties" = list(
            "gsoa_selected" = df_enrichs,
            "splicing_dependecy_stats" = spldep_stats,
            "cosmic_comparison" = cosmic_comparison
        ),
        "model_validation" = list(
            "gene_mutation_frequency" = gene_mut_freq
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
    save_plt(plts, "model_sel-gt_lists_vs_lr_pvalue-random", ".pdf", figs_dir, width=8.5, height=5)
    save_plt(plts, "model_sel-gt_lists_vs_lr_pvalue-lowest", ".pdf", figs_dir, width=8.5, height=5)
    save_plt(plts, "model_sel-lr_pvalue", ".pdf", figs_dir, width=4, height=4)
    save_plt(plts, "model_sel-tpr_vs_fpr", ".pdf", figs_dir, width=5, height=5)
    ## pearson selection
    save_plt(plts, "model_sel-pearson_corr", ".pdf", figs_dir, width=5, height=5)
    save_plt(plts, "model_sel-pearson_corr_vs_spearman", ".pdf", figs_dir, width=5, height=5)
    ## roc analysis
    save_plt(plts, "model_sel-roc_analysis-fpr_vs_tpr-lr_pvalue-line", ".pdf", figs_dir, width=5, height=6)
    save_plt(plts, "model_sel-roc_analysis-recall_vs_precision-lr_pvalue-line", ".pdf", figs_dir, width=5, height=6)
    save_plt(plts, "model_sel-roc_analysis-fpr_vs_tpr-pearson-line", ".pdf", figs_dir, width=5, height=6)
    save_plt(plts, "model_sel-roc_analysis-recall_vs_precision-pearson-line", ".pdf", figs_dir, width=5, height=6)
    ## set overlaps: genes with known cancer-driver exons, cosmic genes, pan essentials
    save_plt(plts, "model_sel-set_overlaps-all-venn", ".pdf", figs_dir, width=5, height=5)
    save_plt(plts, "model_sel-set_overlaps-selected-venn", ".pdf", figs_dir, width=5, height=5)
    ## selection variables vs n observations
    save_plt(plts, "model_sel-n_obs_vs_pvalue-scatter", ".pdf", figs_dir, width=5, height=5)
    save_plt(plts, "model_sel-n_obs_vs_pearson-scatter", ".pdf", figs_dir, width=5, height=5)
    
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
    ## tumorigenesis
    save_plt(plts, "model_prop-tumorigenesis-scatter", ".pdf", figs_dir, width=6, height=6)
    save_plt(plts, "model_prop-tumorigenesis_vs_prot_imp", ".pdf", figs_dir, width=5, height=5)
    save_plt(plts, "model_prop-harm_scores-scatter", ".pdf", figs_dir, width=8, height=7)
    save_plt(plts, "model_prop-harm_scores_vs_prot_imp", ".pdf", figs_dir, width=8, height=7)

    # model validation
    save_plt(plts, "model_val-mutation_gene_count", ".pdf", figs_dir, width=8, height=8)
    save_plt(plts, "model_val-mutation_gene_frequency", ".pdf", figs_dir, width=8, height=8)
    save_plt(plts, "model_val-mutation_gene_frequency_silent_norm", ".pdf", figs_dir, width=8, height=8)
    
    # events oi
    genes_oi = c("KRAS","SMNDC1","NUMB","NUMB_lung")
    lapply(genes_oi, function(gene_oi){
        lapply(c("spldep_vs_rnai","genexpr_vs_rnai","genexpr_vs_splicing","splicing_vs_rnai"), 
               function(comparison){
            plt_name = sprintf("events_oi-%s-%s",gene_oi,comparison)
            save_plt(plts, plt_name, ".pdf", figs_dir, width=4, height=4)            
        })
    })
    
    # comparison cosmic
    save_plt(plts, "cosmic_comparison-overlap-venn", ".pdf", figs_dir, width=4, height=4)
    save_plt(plts, "cosmic_comparison-reactome-bar", ".pdf", figs_dir, width=8, height=18)
    
    # RNAi (shRNA) QC
    save_plt(plts, "rnai_qc-n_targeted_exons_per_gene-hist", ".pdf", figs_dir, width=5, height=2.5)
    save_plt(plts, "rnai_qc-exon_psi_median_vs_std_vs_shrna_targeted-scatter", ".pdf", figs_dir, width=7, height=7)
    save_plt(plts, "rnai_qc-exon_psi_median_vs_shrna_targeted-bar", ".pdf", figs_dir, width=7, height=3)

    # RNAseq QC
    save_plt(plts, "rnaseq_qc-read_count_vs_n_detected_events-scatter", ".pdf", figs_dir, width=5, height=5)
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
        make_option("--gene_mut_freq_file", type="character"),
        make_option("--msigdb_dir", type="character"),
        make_option("--protein_impact_file", type="character"),
        make_option("--cosmic_genes_file", type="character"),
        make_option("--spldep_file", type="character"),
        make_option("--rnai_file", type="character"),
        make_option("--genexpr_file", type="character"),
        make_option("--splicing_file", type="character"),
        make_option("--cancer_events_file", type="character"),
        make_option("--ascanceratlas_file", type="character"),
        make_option("--event_info_file", type="character"),
        make_option("--metadata_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}


main = function(){
    args = parseargs()
    
    models_file = args[["models_file"]]
    ccle_stats_file = args[["ccle_stats_file"]]
    gene_mut_freq_file = args[["gene_mut_freq_file"]]
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
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    models = read_tsv(models_file)
    ccle_stats = read_tsv(ccle_stats_file)
    gene_mut_freq = read_tsv(gene_mut_freq_file) %>% dplyr::rename(GENE=Hugo_Symbol)
    ontologies = load_ontologies(msigdb_dir, protein_impact_file, cosmic_genes_file)
    spldep = read_tsv(spldep_file)
    rnai = read_tsv(rnai_file)
    genexpr = read_tsv(genexpr_file)
    splicing = read_tsv(splicing_file)
    cancer_events = read_tsv(cancer_events_file)
    ascanceratlas = read_tsv(ascanceratlas_file)
    event_info = read_tsv(event_info_file)
    metadata = read_tsv(metadata_file)
    shrna_seqs = read_csv(shrna_seqs_file)
    shrna_mapping_ensembl = read_tsv(shrna_mapping_ensembl_file)
    shrna_mapping_vastdb = read_tsv(shrna_mapping_vastdb_file)
    #genome_annot = ensembldb::EnsDb(genome_annot_file)
    models_achilles = read_tsv(models_achilles_file)
    crispr_essentials = read_csv(crispr_essentials_file)
    ccle_info = read_tsv(ccle_info_file)
    
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
    
    cols_oi = c('EVENT','ENSEMBL','GENE','length','source_authors','source_year',
                'source_doi','description','event_id','source')
    cancer_events = cancer_events %>%
        mutate(source="HandCurated") %>%
        bind_rows(ascanceratlas) %>%
        dplyr::select(all_of(cols_oi)) %>%
        distinct()
    
    # prep pan essential
    crispr_essentials = crispr_essentials %>%
        separate(Essentials, into=c("GENE","ENTREZ"), sep=" ") %>%
        pull(GENE)
    
    rnai_essentials = rnai %>%
        pivot_longer(-index, names_to="DepMap_ID", values_to="demeter2") %>%
        drop_na() %>%
        group_by(DepMap_ID) %>%
        mutate(
            ranking = rank(demeter2),
            perc = ranking / n()
        ) %>%
        ungroup() %>%
        group_by(index) %>%
        summarize(
            percentile = quantile(perc, probs=0.9)
        ) %>%
        ungroup()

    density_est = rnai_essentials %>% pull(percentile) %>% density()
    minima = density_est$x[which(diff(sign(diff(density_est$y))) > 0) + 1]
    rnai_essentials = rnai_essentials %>%
        filter(percentile < min(minima)) %>%
        pull(index)
    
    pan_essentials = union(crispr_essentials, rnai_essentials)
    
    # log normalize gene expression
    genexpr = genexpr %>% mutate_at(vars(-("ID")), function(x){ log2(x+1) })

    # prep ccle info
    metadata = metadata %>%
        left_join(
            ccle_info %>% 
            filter(library_strategy=="RNA-Seq") %>%
            distinct(run_accession, read_count),
            by="run_accession"
        ) %>%
        left_join(
            splicing %>%
            summarize_if(is.numeric, ~ sum(is.finite(.))) %>%
            pivot_longer(everything(), names_to="DepMap_ID", values_to="n_detected_events"),
            by="DepMap_ID"
        )
        
    # add info to models
    models = models %>% 
        mutate(event_gene = paste0(EVENT,"_",GENE),
               event_type = gsub("Hsa","",gsub("[^a-zA-Z]", "",EVENT)),
               is_selected = pearson_correlation_mean > THRESH_CORR & 
                             lr_pvalue < THRESH_LR_PVALUE,
               is_selected = replace_na(is_selected, FALSE),
               is_known_driver = EVENT %in% (cancer_events %>% pull(EVENT)),
               is_pan_essential = GENE %in% pan_essentials,
               in_cosmic = GENE %in% ontologies[["cosmic"]][["gene"]]
        ) %>% 
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
    eval_corr = thresh_eval_corr(models, rnai, spldep, ctl_neg, ctl_pos)
    
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
    cosmic_comparison = data.frame(gene = available_genes) %>%
        mutate(
            has_cancer_driver_exon = gene %in% selected_genes,
            is_cancer_driver = gene %in% ontologies[["cosmic"]][["gene"]]
        ) %>% 
        left_join(ontologies[["reactome"]], by="gene")
    cosmic_test = cosmic_comparison %>%
        distinct(gene, has_cancer_driver_exon, is_cancer_driver) %>%
        with(table(has_cancer_driver_exon, is_cancer_driver)) %>%
        fisher.test()
    # 74 cancer driver genes have cancer-driver exons
    
    # add VastDB IDs to ENSEMBL exons
    #     event_coords = data.frame(
    #         EVENT = event_info[["EVENT"]],
    #         vastdb_chr = event_info[["COORD_o"]] %>% gsub(":.*","",.) %>% gsub("chr","",.),
    #         vastdb_start = event_info[["COORD_o"]] %>% gsub(".*:","",.) %>% gsub("-.*","",.) %>% as.integer(),
    #         vastdb_end = event_info[["COORD_o"]] %>% gsub(".*:","",.) %>% gsub(".*-","",.) %>% as.integer(),
    #         vastdb_strand = event_info[["REF_CO"]] %>% gsub(".*:","",.)
    #     ) %>%
    #     filter(str_detect(EVENT,"EX"))
    #     shrna_mapping = shrna_mapping %>%
    #         left_join(
    #             event_coords, by=c(
    #                 "start_exon_genomic"="vastdb_start", 
    #                 "end_exon_genomic"="vastdb_end", 
    #                 "strand"="vastdb_strand")
    #         ) %>%
    #         mutate(is_selected = EVENT %in% selected_events)
    
    # plot
    plts = make_plots(
        models, rnai_stats, cancer_events, 
        eval_pvalue, eval_corr, 
        enrichment, spldep_stats, harm_stats, 
        gene_mut_freq,
        rnai, spldep, splicing, genexpr, metadata,
        cosmic_comparison, shrna_mapping
    )

    # make figdata
    figdata = make_figdata(
        models, rnai_stats, cancer_events, 
        eval_pvalue, eval_corr, 
        enrichment, spldep_stats, harm_stats,
        gene_mut_freq,
        rnai, spldep, splicing, genexpr, metadata,
        cosmic_comparison, shrna_mapping
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