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

PAL_SINGLE_ACCENT = "orange"
PAL_SINGLE_LIGHT = "#6AC2BF"
PAL_SINGLE_DARK = "#716454"
PAL_DUAL = c(PAL_SINGLE_DARK, PAL_SINGLE_LIGHT)
LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

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


plot_eda_associations = function(models, drug_screen){
    top_n = 25
    X = models %>% 
        left_join(drug_screen %>% distinct(DRUG_ID,DRUG_NAME) %>% drop_na(), by="DRUG_ID")
    
    plts = list()
    
    # what are the distributions of p-values and FDR?
    plts[["associations-lr_pvalues"]] = X %>% 
        gghistogram(x="lr_pvalue", fill=PAL_SINGLE_LIGHT, color=NA, bins=100) + 
        geom_vline(xintercept=median(X[["lr_pvalue"]], na.rm=TRUE), 
                   linetype="dashed", size=LINE_SIZE) +
        labs(x="LR Test p-value")
    
    plts[["associations-lr_fdr"]] = X %>% 
        gghistogram(x="lr_padj", fill=PAL_SINGLE_LIGHT, color=NA, bins=100) + 
        geom_vline(xintercept=median(X[["lr_padj"]], na.rm=TRUE), 
                   linetype="dashed", size=LINE_SIZE) +
        labs(x="LR Test FDR")
    
    # how do number of observations relate to association coefficient?
    plts[["associations-nobs_vs_coef"]] = X %>% 
        filter(lr_padj < THRESH_FDR) %>% 
        ggplot(aes(x=n_obs, y=spldep_coefficient)) + 
        geom_scattermore(pixels=c(1000,1000), pointsize=8, alpha=0.2, color=PAL_SINGLE_DARK) + 
        theme_pubr() + 
        geom_vline(xintercept=THRESH_NOBS, linetype="dashed", size=LINE_SIZE) + 
        labs(x="N. Observations", y="Assoc. Coefficient")
    
    # how many drugs and events are significantly associated?
    plts[["associations-drug_counts"]] = X %>% 
        filter(lr_padj < THRESH_FDR & n_obs > THRESH_NOBS) %>% 
        count(drug_screen, DRUG_NAME) %>% 
        group_by(DRUG_NAME) %>%
        mutate(total=sum(n)) %>%
        arrange(total) %>%
        ggbarplot(x="DRUG_NAME", y="n", fill="drug_screen", palette=PAL_DUAL, color=NA) + 
        labs(x="Drug", y="N. Significant Associations", fill="Drug Screen")
    
    plts[["associations-top_drug_counts"]] = X %>% 
        filter(lr_padj < THRESH_FDR & n_obs > THRESH_NOBS) %>% 
        count(drug_screen, DRUG_NAME) %>% 
        group_by(DRUG_NAME) %>%
        mutate(total=sum(n)) %>%
        arrange(-total) %>%
        ungroup() %>%
        mutate(index=as.numeric(factor(DRUG_NAME, levels=unique(DRUG_NAME)))) %>%
        filter(index <= top_n) %>%
        ggbarplot(x="DRUG_NAME", y="n", fill="drug_screen", palette=PAL_DUAL, color=NA) + 
        labs(x="Drug", y="N. Significant Associations", fill="Drug Screen",
             title=sprintf("Top %s", top_n)) +
        coord_flip()
    
    plts[["associations-spldep_counts"]] = X %>% 
        filter(lr_padj < THRESH_FDR & n_obs > THRESH_NOBS) %>% 
        count(drug_screen, event_gene) %>% 
        group_by(event_gene) %>%
        mutate(total=sum(n)) %>%
        arrange(total) %>%
        ggbarplot(x="event_gene", y="n", fill="drug_screen", palette=PAL_DUAL, color=NA) + 
        labs(x="Event & Gene", y="N. Significant Associations", fill="Drug Screen")
    
    plts[["associations-top_spldep_counts"]] = X %>% 
        filter(lr_padj < THRESH_FDR & n_obs > THRESH_NOBS) %>% 
        count(drug_screen, event_gene) %>% 
        group_by(event_gene) %>%
        mutate(total=sum(n)) %>%
        arrange(-total) %>%
        ungroup() %>%
        mutate(index=as.numeric(factor(event_gene, levels=unique(event_gene)))) %>%
        filter(index <= top_n) %>%
        ggbarplot(x="event_gene", y="n", fill="drug_screen", palette=PAL_DUAL, color=NA) + 
        labs(x="Event & Gene", y="N. Significant Associations", fill="Drug Screen", 
             title=sprintf("Top %s", top_n)) +
        coord_flip()
    
    # - overlaps with drugs in GDSC1 and GDSC2
    samples_oi = list(
        "GDSC1" = drug_screen %>% filter(DATASET=="GDSC1") %>% pull(DRUG_ID) %>% unique() %>% na.omit(),
        "GDSC2" = drug_screen %>% filter(DATASET=="GDSC2") %>% pull(DRUG_ID) %>% unique() %>% na.omit()
    )
    m = samples_oi %>% 
        list_to_matrix() %>% 
        make_comb_mat()
    plts[["associations-agreement_between-upset"]] = m %>%
        UpSet(comb_order = order(comb_size(m)), 
              comb_col = PAL_SINGLE_DARK,
              top_annotation = upset_top_annotation(m, gp = gpar(fill = PAL_SINGLE_DARK, col=NA)),
              right_annotation = upset_right_annotation(m, gp = gpar(fill = PAL_SINGLE_DARK, col=NA))) %>%
            draw() %>%
            grid.grabExpr() %>%
            as.ggplot()
    
    # - agreement between GDSC1 and GDSC2
    drugs_oi = X %>% 
        distinct(DRUG_ID,drug_screen) %>% 
        count(DRUG_ID) %>% 
        filter(n>1) %>% 
        pull(DRUG_ID)
            
    correls = X %>% 
        filter(DRUG_ID %in% drugs_oi) %>% 
        pivot_wider(id_cols=c("DRUG_ID","DRUG_NAME","event_gene"), 
                    names_from="drug_screen", 
                    values_from="spldep_coefficient") %>% 
        group_by(DRUG_ID) %>% 
        summarize(correl=cor(GDSC1, GDSC2, method="pearson", use="pairwise.complete.obs"))
    
    plts[["associations-agreement_between-correlation"]] = correls %>% 
        gghistogram(x="correl", fill=PAL_SINGLE_DARK, bins=12, color=NA) + 
        labs(x=sprintf("Pearson Correlation (n=%s)", nrow(correls)),
             y="Count")
    
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
        UpSet(comb_order = order(comb_size(m)), 
              comb_col = PAL_SINGLE_DARK,
              top_annotation = upset_top_annotation(m, gp = gpar(fill = PAL_SINGLE_DARK, col=NA)),
              right_annotation = upset_right_annotation(m, gp = gpar(fill = PAL_SINGLE_DARK, col=NA))) %>%
            draw() %>%
            grid.grabExpr() %>%
            as.ggplot()
    
    return(plts)
}


plot_targets = function(models, drug_targets){
    X = models %>%
        left_join(
            drug_targets %>% distinct(DRUG_ID, TARGET) %>% mutate(is_target=TRUE),
            by = c("DRUG_ID","GENE"="TARGET")
        ) %>%
        mutate(is_target=replace_na(is_target, FALSE))
    
    plts = list()
    
    # - p-value and betas of drugs with target (sign. associated or not)
    plts[["targets-violin-lr_pvalue"]] = X %>% 
        ggviolin(x="is_target", y="lr_pvalue", trim=TRUE, fill="is_target", color=NA, palette=PAL_DUAL) + 
        geom_boxplot(width=0.1, outlier.size=0.1) + 
        guides(fill="none") + 
        stat_compare_means(method="wilcox", size=FONT_SIZE, family=FONT_FAMILY) + 
        labs(x="Is Target", y="LR p-value")
    
    plts[["targets-violin-lr_padj"]] = X %>% 
        ggviolin(x="is_target", y="lr_padj", trim=TRUE, fill="is_target", color=NA, palette=PAL_DUAL) + 
        geom_boxplot(width=0.1, outlier.size=0.1) + 
        guides(fill="none") + 
        stat_compare_means(method="wilcox", size=FONT_SIZE, family=FONT_FAMILY) + 
        labs(x="Is Target", y="LR FDR")
    
    plts[["targets-violin-spldep_coef"]] = X %>% 
        ggviolin(x="is_target", y="spldep_coefficient", trim=TRUE, color=NA, 
                 fill="is_target", palette=PAL_DUAL) + 
        geom_boxplot(width=0.1, outlier.size=0.1) + 
        guides(fill="none") + 
        stat_compare_means(method="wilcox", size=FONT_SIZE, family=FONT_FAMILY) + 
        labs(x="Is Target", y="Spl. Dep. Coefficient")
    
    plts[["targets-violin-spldep_coef_filt"]] = X %>% 
        filter(n_obs > THRESH_NOBS) %>%
        ggviolin(x="is_target", y="spldep_coefficient", trim=TRUE, color=NA, 
                 fill="is_target", palette=PAL_DUAL) + 
        geom_boxplot(width=0.1, outlier.size=0.1) + 
        guides(fill="none") + 
        stat_compare_means(method="wilcox", size=FONT_SIZE, family=FONT_FAMILY) + 
        labs(x="Is Target", y="Spl. Dep. Coefficient")
    
    # - from drugs with annotated target with a cancer-driver exon, how many sign. assocs? not signif?
    plts[["targets-bar-n_signifs"]] = X %>%
        mutate(is_sel = lr_padj < THRESH_FDR & n_obs > THRESH_NOBS,
               is_sel = replace_na(is_sel, FALSE)) %>%
        count(is_target, is_sel) %>%
        mutate(lab = paste0(is_target,"\n&\n",is_sel)) %>%
        ggbarplot(x="lab", y="n", fill="lab", color=NA, 
                  palette=get_palette(PAL_DUAL, 4), 
                  label=TRUE, lab.size=FONT_SIZE, lab.family=FONT_FAMILY) +
        guides(fill="none") +
        yscale("log10", .format=TRUE) +
        labs(x="Is Target & Selected", y="Count")
    
    return(plts)
}


plot_ppi = function(models, shortest_paths){
    X = shortest_paths %>%
        left_join(models, by = c("DRUG_ID","source"="GENE")) %>%
        mutate(is_sel = lr_padj < THRESH_FDR & n_obs > THRESH_NOBS,
               path_lab = as.character(shortest_path_length),
               path_lab = ifelse(shortest_path_length==1e6,"N.R.",path_lab),
               path_lab = ifelse(shortest_path_length>=5,">=5",path_lab),
               path_lab = factor(path_lab, levels=c("0","1","2","3","4",">=5","N.R.")))
    
    plts = list()
    
    # - Do significant associations tend to be closer to the drug target?
    plts[["ppi-bar-counts"]] = X %>% 
        count(is_sel, path_lab) %>% 
        ggbarplot(x="path_lab", y="n", fill="is_sel", color=NA,
                  palette=PAL_DUAL, position=position_dodge(0.7),
                  label = TRUE, lab.size=FONT_SIZE, lab.family=FONT_FAMILY) + 
        labs(x="Shortest Path Length to Drug Target(s)", y="Count", fill="Selected")
    
    plts[["ppi-bar-prop"]] = X %>% 
        count(is_sel, path_lab) %>% 
        group_by(is_sel) %>%
        mutate(prop = round(n / sum(n), 3)) %>%
        ungroup() %>%
        ggbarplot(x="path_lab", y="prop", fill="is_sel", color=NA,
                  palette=PAL_DUAL, position=position_dodge(0.7),
                  label = TRUE, lab.size=FONT_SIZE, lab.family=FONT_FAMILY) + 
        labs(x="Shortest Path Length to Drug Target(s)", y="Proportion", fill="Selected")
    
    plts[["ppi-bar-rel_prop"]] = X %>% 
        count(is_sel, path_lab) %>% 
        group_by(is_sel) %>%
        mutate(prop = n / sum(n)) %>%
        ungroup() %>%
        group_by(path_lab) %>%
        mutate(rel_prop = prop / sum(prop)) %>%
        ungroup() %>%
        ggbarplot(x="path_lab", y="rel_prop", fill="is_sel", color=NA,
                  palette=PAL_DUAL) + 
        labs(x="Shortest Path Length to Drug Target(s)", y="Rel. Proportion", fill="Selected")
    
    # - Are FDRs or effect sizes of significant associations informative of closeness to drug target?
    plts[["ppi-box-lr_padj"]] = X %>% 
        ggboxplot(x="path_lab", y="lr_padj", fill="is_sel", 
                  palette=PAL_DUAL, outlier.size=0.1) + 
        labs(x="Shortest Path Length to Drug Target(s)", y="LR FDR", fill="Selected") +
        facet_wrap(~is_sel, ncol=2, scales="free_y") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    plts[["ppi-box-spldep_coef"]] = X %>%
        ggboxplot(x="path_lab", y="spldep_coefficient", fill="is_sel", 
                  palette=PAL_DUAL, outlier.size=0.1) + 
        labs(x="Shortest Path Length to Drug Target(s)", y="Spl. Dep. Coef.", fill="Selected") +
        facet_wrap(~is_sel, ncol=2, scales="free_y") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    return(plts)
}


plot_mediators = function(spldep_models, models, shortest_paths_simple, drug_targets){
    # to get best candidates biologically relevant, we will only consider
    # exons whose PSI coefficient is larger than the TPM coefficient
    # in selected linear models
    events_oi = spldep_models %>% 
        filter(abs(event_coefficient_mean) > abs(gene_coefficient_mean)) %>% 
        pull(EVENT)
    
    X = models %>%
        left_join(shortest_paths_simple, by=c("DRUG_ID","GENE"="source")) %>%
        left_join(drug_targets %>% distinct(DRUG_ID,DRUG_NAME), by="DRUG_ID") %>%
        mutate(is_large = EVENT %in% events_oi,
               is_sel = lr_padj < THRESH_FDR & n_obs > THRESH_NOBS,
               is_sel = replace_na(is_sel, FALSE),
               path_lab = as.character(shortest_path_length),
               path_lab = ifelse(shortest_path_length==1e6,"N.R.",path_lab),
               path_lab = ifelse(shortest_path_length>=5,">=5",path_lab),
               path_lab = replace_na(shortest_path_length,"NA"),
               is_target = shortest_path_length=="0",
               is_target = replace_na(is_target,FALSE))
    
    
    plts = list()
    
    # - how many events are we considering?
    plts[["mediators-bar-n_selected"]] = X %>% 
        count(is_large, is_sel) %>% 
        mutate(lab=paste0(is_large,"\n&\n",is_sel)) %>%
        ggbarplot(x="lab", y="n", fill="lab", color=NA, 
                  palette=get_palette(PAL_DUAL, 4), 
                  label=TRUE, lab.size=FONT_SIZE, lab.family=FONT_FAMILY) +
        guides(fill="none") +
        yscale("log10", .format=TRUE) +
        labs(x="Larger Splicing Effect & Selected", y="Count")
    
    # what is the best significant association with larger PSI effect?    
    # - mediator on target
    assocs_oi = X %>%
        filter(is_sel & is_large & is_target) %>%
        distinct(EVENT,DRUG_ID)
    
    x = assocs_oi %>% 
        left_join(spldep_models %>% 
        distinct(EVENT, event_coefficient_mean), by="EVENT") %>% 
        left_join(X, by=c("EVENT","DRUG_ID")) %>%
        mutate(event_gene = as.factor(event_gene),
               name = reorder_within(event_gene, spldep_coefficient, DRUG_NAME))
    
    plts[["mediators-on_target-spldep_coef"]] = x %>% 
        mutate(lab = round(spldep_coefficient, 3)) %>%
        ggbarplot(x="name", y="spldep_coefficient", 
              color="drug_screen", fill=PAL_SINGLE_LIGHT, 
              position=position_dodge(0.7)) +
        color_palette(c("black","white")) + 
        geom_text(aes(y=spldep_coefficient+0.05, label=lab, group=drug_screen), 
                  position=position_dodge(0.7),
                  size=FONT_SIZE, family=FONT_FAMILY) +
        geom_hline(yintercept=0, linetype="dashed", size=LINE_SIZE) +
        scale_x_reordered() +
        labs(x="Event & Gene" , y="Spl. Dep. Coefficient", color="Drug Screen") +
        coord_flip()
    
    plts[["mediators-on_target-event_coef"]] = x %>% 
        distinct(name, event_coefficient_mean) %>%
        mutate(lab = round(event_coefficient_mean, 3)) %>%
        ggbarplot(x="name", y="event_coefficient_mean", 
                  color=NA, fill=PAL_SINGLE_LIGHT,
                  position=position_dodge(0.7)) +
        geom_text(aes(y=sign(event_coefficient_mean)*(abs(event_coefficient_mean)+0.015), label=lab),
                  position=position_dodge(0.7),
                  size=FONT_SIZE, family=FONT_FAMILY) + 
        geom_hline(yintercept=0, linetype="dashed", size=LINE_SIZE) +
        scale_x_reordered() +
        labs(x="Event & Gene" , y="PSI Coefficient") +
        coord_flip()
    
    # - mediator up/downstream (off target)    
    plts[["mediators-on_target-top_assocs_spldep"]] = assocs_oi %>% 
        left_join(X, by="DRUG_ID") %>% 
        filter(is_sel & is_large) %>% 
        group_by(DRUG_NAME) %>% 
        slice_max(spldep_coefficient, n=10) %>% 
        ungroup() %>%
        mutate(DRUG_NAME = as.factor(DRUG_NAME),
               name = reorder_within(event_gene, spldep_coefficient, DRUG_NAME)) %>%
        ggbarplot(x="name", y="spldep_coefficient", 
                  color="drug_screen", fill="is_target", palette=PAL_DUAL, 
                  position=position_dodge(0.7)) + 
        geom_text(aes(y=spldep_coefficient+0.08, label=path_lab, group=drug_screen), 
                  position=position_dodge(0.7),
                  size=FONT_SIZE, family=FONT_FAMILY) + 
        color_palette(c("black","white")) + 
        facet_wrap(~DRUG_NAME, scales="free") + 
        scale_x_reordered() +
        coord_flip() +
        labs(x="Event & Gene", y="Spl. Dep. Coef.",, fill="Is Target", color="Drug Screen") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    plts[["mediators-on_target-top_assocs_padj"]] = assocs_oi %>% 
        left_join(X, by="DRUG_ID") %>% 
        filter(is_sel & is_large) %>% 
        group_by(DRUG_NAME) %>% 
        slice_min(lr_padj, n=10) %>% 
        ungroup() %>%
        mutate(DRUG_NAME = as.factor(DRUG_NAME),
               log_padj = -log10(lr_padj),
               name = reorder_within(event_gene, log_padj, DRUG_NAME)) %>%
        ggbarplot(x="name", y="log_padj", 
                  color="drug_screen", fill="is_target", palette=PAL_DUAL, 
                  position=position_dodge(0.7)) + 
        geom_text(aes(y=log_padj+0.5, label=path_lab, group=drug_screen), 
                  position=position_dodge(0.7),
                  size=FONT_SIZE, family=FONT_FAMILY) + 
        color_palette(c("black","white")) + 
        facet_wrap(~DRUG_NAME, scales="free") + 
        scale_x_reordered() +
        labs(x="Event & Gene", y="LR -log10(FDR)", fill="Is Target", color="Drug Screen") +
        coord_flip() +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    # - examples exon-drug synergy
    events_oi = c("HsaEX0034998", # KRAS (reported in literature)
                  "HsaEX0038399", # MDM2
                  "HsaEX1013695") # EZH2
    
    x = X %>% 
        drop_na(spldep_coefficient) %>%
        filter(is_sel) %>%
        group_by(DRUG_ID, drug_screen) %>%
        # need to normalize at drug level
        mutate(norm_coef = order(spldep_coefficient, decreasing=TRUE) / n()) %>%
        ungroup()
        
    plts[["mediators-examples-synergy"]] = x %>% 
        filter(EVENT %in% events_oi) %>% 
        group_by(event_gene) %>% 
        slice_min(norm_coef, n=10) %>% 
        ungroup() %>%
        mutate(event_gene = as.factor(event_gene),
               name = reorder_within(DRUG_NAME, norm_coef, event_gene)) %>%
        ggbarplot(x="name", y="norm_coef", 
                  color="drug_screen", fill="is_target", palette=PAL_DUAL, 
                  position=position_dodge(0.7)) + 
        geom_text(aes(y=norm_coef+0.1, label=path_lab, group=drug_screen), 
                  position=position_dodge(0.7),
                  size=FONT_SIZE, family=FONT_FAMILY) + 
        color_palette(c("black","white")) + 
        facet_wrap(~event_gene, scales="free") + 
        scale_x_reordered() +
        labs(x="Drug", y="Norm. Spl. Dep. Coef.", fill="Is Target", color="Drug Screen") +
        coord_flip() +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
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
        gghistogram(x="n", fill=PAL_SINGLE_DARK, color=NA) +
        facet_wrap(~drug_screen+in_demeter2) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    plts[["drug_rec-correlations"]] = corrs %>% 
        ggplot(aes(x=drug_screen, y=correlation, 
                   group=interaction(drug_screen,in_demeter2))) +
        geom_violin(aes(fill=in_demeter2), color=NA) +
        geom_boxplot(width=0.1, outlier.size=0.1, position=position_dodge(0.9)) +
        geom_hline(yintercept=0, linetype="dashed") +
        fill_palette(PAL_DUAL) + 
        labs(x="Drug Screen", y="Pearson Correlation", fill="In Demeter2") +
        geom_text(aes(y=1.00, label=n), 
                  corrs %>% count(in_demeter2, drug_screen), 
                  position=position_dodge(0.9),
                  size=FONT_SIZE, family=FONT_FAMILY) +
        theme_pubr()
    
    plts[["drug_rec-npos_vs_corrs"]] = corrs %>% 
        ggscatter(x="n_pos", y="correlation", alpha=0.5, size=1,
                  color="in_demeter2", palette=PAL_DUAL) + 
        facet_wrap(~drug_screen+in_demeter2, scales="free_x") +
        labs(x="N. Associations per Sample", y="Pearson Correlation", color="In Demeter2") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))

    # best and worse correlations
    samples_oi = corrs %>% 
        group_by(drug_screen) %>%
        drop_na(correlation) %>%
        arrange(correlation) %>% 
        filter(row_number()==1 | row_number()==n()) %>%
        left_join(X, by=c("drug_screen","sample"))
    
    plts[["drug_rec-best_worse"]] = samples_oi %>% 
        ggscatter(x="real_ic50", y="predicted_ic50", size=0.5, color=PAL_SINGLE_DARK) + 
        facet_wrap(drug_screen~sample, scales="free") +
        stat_cor(method="pearson", size=FONT_SIZE) + 
        labs(x="Real log(IC50)", y="Predicted log(IC50)") + 
        theme_pubr(border=TRUE) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    return(plts)
}


make_plots = function(models, drug_screen,
                      drug_targets, shortest_paths, shortest_paths_simple,
                      spldep_models, estimated_drug_response){
    plts = list(
        plot_eda_associations(models, drug_screen),
        plot_targets(models, drug_targets),
        plot_ppi(models, shortest_paths),
        plot_mediators(spldep_models, models, shortest_paths_simple, drug_targets),
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
                    font.x=8, font.y=8, font.legend=6,
                    font.tickslab=6, font.family=FONT_FAMILY)    
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
    save_plt(plts, "associations-agreement_between-upset", ".pdf", figs_dir, width=7, height=6)
    save_plt(plts, "associations-agreement_between-correlation", ".pdf", figs_dir, width=5, height=5)
    save_plt(plts, "associations-overlaps_screens-upset", ".pdf", figs_dir, width=8, height=7)
    
    # targets
    save_plt(plts, "targets-violin-lr_pvalue", ".pdf", figs_dir, width=5, height=5)
    save_plt(plts, "targets-violin-lr_padj", ".pdf", figs_dir, width=5, height=5)
    save_plt(plts, "targets-bar-n_signifs", ".pdf", figs_dir, width=5, height=5)
    
    # ppi
    save_plt(plts, "ppi-bar-counts", ".pdf", figs_dir, width=5, height=6)
    save_plt(plts, "ppi-bar-prop", ".pdf", figs_dir, width=5, height=6)
    save_plt(plts, "ppi-bar-rel_prop", ".pdf", figs_dir, width=5, height=6)
    save_plt(plts, "ppi-box-lr_padj", ".pdf", figs_dir, width=8, height=6)
    save_plt(plts, "ppi-box-spldep_coef", ".pdf", figs_dir, width=8, height=6)
    
    # mediators
    save_plt(plts, "mediators-bar-n_selected", ".pdf", figs_dir, width=5, height=5)
    save_plt(plts, "mediators-on_target-spldep_coef", ".pdf", figs_dir, width=6, height=5)
    save_plt(plts, "mediators-on_target-event_coef", ".pdf", figs_dir, width=6, height=5)
    save_plt(plts, "mediators-on_target-top_assocs_spldep", ".pdf", figs_dir, width=10, height=7)
    save_plt(plts, "mediators-on_target-top_assocs_padj", ".pdf", figs_dir, width=10, height=7)
    save_plt(plts, "mediators-examples-synergy", ".pdf", figs_dir, width=12, height=7)
    
    # drug recommendations
    save_plt(plts, "drug_rec-n_screens", ".pdf", figs_dir, width=6, height=7)
    save_plt(plts, "drug_rec-correlations", ".pdf", figs_dir, width=5, height=5)
    save_plt(plts, "drug_rec-npos_vs_corrs", ".pdf", figs_dir, width=6, height=9)
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
    rnai = read_tsv(rnai_file)
    estimated_response = read_tsv(estimated_response_file)
    drug_screen = load_drug_screens(drug_screens_dir)
    spldep_ccle = read_tsv(spldep_ccle_file) %>%
        filter(index %in% unique(models[["EVENT"]]))
    paths_real = read_tsv(paths_real_file)
    spldep_models = read_tsv(spldep_models_file) %>%
        filter(EVENT %in% unique(models[["EVENT"]]))
    
    # prep inputs
    ## IC50
    estimated_response = estimated_response %>% mutate(in_demeter2 = sample %in% colnames(rnai))
    ## PPI
    shortest_paths = paths_real %>% mutate(shortest_path_length=replace_na(shortest_path_length,1e6))
    shortest_paths_simple = shortest_paths %>%
        distinct(source, shortest_path_length, DRUG_ID) %>%
        group_by(DRUG_ID, source) %>%
        slice_min(shortest_path_length, n=1) %>%
        ungroup()
    
    # make plots
    plts = make_plots(models, drug_screen,
                      drug_targets, shortest_paths, shortest_paths_simple,
                      spldep_models, estimated_drug_response)
    
    # make figdata
    figdata = make_figdata(models, drug_screen,
                      drug_targets, shortest_paths, shortest_paths_simple,
                      spldep_models, estimated_drug_response)
    
    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}