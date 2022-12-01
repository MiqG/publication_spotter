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

require(optparse)
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
require(ggtext)
require(ggnewscale)

require(ggraph)
require(igraph)
require(graphlayouts)

require(statebins)

ROOT = here::here()
source(file.path(ROOT,"src","R","utils.R"))

# variables
THRESH_FDR = 0.1
THRESH_PVALUE = 0.05
THRESH_NOBS = 20
RANDOM_SEED = 1234

# formatting
PAL_SINGLE_ACCENT = "orange"
PAL_SINGLE_LIGHT = "#6AC2BF"
PAL_SINGLE_DARK = "#716454"
PAL_DUAL = c(PAL_SINGLE_DARK, PAL_SINGLE_LIGHT)
PAL_IS_TARGET = get_palette(PAL_DUAL, 4)[3:4]
PAL_IN_DEMETER2 = setNames(c(PAL_SINGLE_LIGHT,"orange"), c(FALSE,TRUE))
LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

PAL_PROT_IMP = setNames(
    get_palette("jco", 10),
    c("NonCoding", "5\' UTR", "ORF disruption (exclusion)", "Alternative protein isoforms", 
      "ORF disruption (inclusion)", "In the CDS (uncertain)", "3\' UTR", "NA", 
      "ORF disruption when splice site is used", "Protein isoform when splice site is used"
    )
)

# Development
# -----------
# RAW_DIR = file.path(ROOT,"data","raw")
# PREP_DIR = file.path(ROOT,"data","prep")
# RESULTS_DIR = file.path(ROOT,"results","exon_drug_interactions")
# MODELS_DIR = file.path(ROOT,"results","model_splicing_dependency")

# metadata_file = file.path(PREP_DIR,"metadata","CCLE.tsv.gz")
# models_file = file.path(RESULTS_DIR,"files","model_summaries_drug_response-EX.tsv.gz")
# drug_targets_file = file.path(PREP_DIR,'drug_screens','drug_targets.tsv.gz')
# rnai_file = file.path(PREP_DIR,"demeter2","CCLE.tsv.gz")
# estimated_response_file = file.path(RESULTS_DIR,"files","estimated_drug_response_by_drug-EX.tsv.gz")
# drug_screens_dir = file.path(PREP_DIR,'drug_screens')
# spldep_ccle_file = file.path(MODELS_DIR,"files","splicing_dependency-EX","mean.tsv.gz")
# paths_real_file = file.path(RESULTS_DIR,"files","ppi","shortest_path_lengths_to_drug_targets-EX.tsv.gz")
# spldep_models_file = file.path(MODELS_DIR,"files","models_gene_dependency-EX","model_summaries.tsv.gz")
# msigdb_dir = file.path(ROOT,"data","raw","MSigDB","msigdb_v7.4","msigdb_v7.4_files_to_download_locally","msigdb_v7.4_GMTs")
# protein_impact_file = file.path(ROOT,"data","raw","VastDB","PROT_IMPACT-hg38-v3.tab.gz")
# splicing_file = file.path(PREP_DIR,"event_psi","CCLE-EX.tsv.gz")
# genexpr_file = file.path(PREP_DIR,"genexpr_tpm","CCLE.tsv.gz")
# snv_file = file.path(RAW_DIR,'DepMap','achilles_ccle','CCLE_mutations.csv')
# ppi_file = file.path(PREP_DIR,'ppi','STRINGDB.tsv.gz')
# event_info_file = file.path(RAW_DIR,"VastDB","EVENT_INFO-hg38_noseqs.tsv")
# gene_info_file = file.path(RAW_DIR,"ENSEMBL","gene_annotation-hg38.tsv.gz")
# figs_dir = file.path(RESULTS_DIR,"figures","model_drug_screens")

# paths_random_file = file.path(RESULTS_DIR,"files","ppi","shortest_path_lengths_to_drug_targets-random.tsv.gz")

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


load_ontologies = function(msigdb_dir, protein_impact_file){
    ontologies = list(
        "reactome" = read.gmt(file.path(msigdb_dir,"c2.cp.reactome.v7.4.symbols.gmt")),
        "hallmarks" = read.gmt(file.path(msigdb_dir,"h.all.v7.4.symbols.gmt")),
        "protein_impact" = read_tsv(protein_impact_file) %>%
            dplyr::rename(EVENT=EventID, term=ONTO) %>%
            dplyr::select(term,EVENT) %>%
            mutate(term_clean=gsub(" \\(.*","",term),
                   term_clean=gsub("ORF disruption upon sequence exclusion",
                                   "ORF disruption (exclusion)",term_clean),
                   term_clean=gsub("ORF disruption upon sequence inclusion",
                                   "ORF disruption (inclusion)",term_clean),
                   term_clean=gsub("In the CDS, with uncertain impact",
                                   "In the CDS (uncertain)",term_clean))
    )
    return(ontologies)
}


evaluate_reactome = function(models, drug_targets, ontology){
    # - are drug target(s) and indirect hits in the same pathway?
    X = models %>%
        left_join(
            drug_targets %>% distinct(DRUG_ID, TARGET) %>% mutate(is_target=TRUE),
            by = c("DRUG_ID","GENE"="TARGET")
        ) %>%
        mutate(is_target=replace_na(is_target, FALSE),
               is_sel = lr_padj<THRESH_FDR & n_obs>THRESH_NOBS) %>%
        filter(!is_target)
    
    ontology = ontology %>%
        group_by(term) %>%
        mutate(gene_set_size=n()) %>%
        ungroup()
    
    ## select pathways containing a drug target
    drug_pathways = drug_targets %>% 
        distinct(DRUG_ID, TARGET) %>% 
        left_join(ontology, by = c("TARGET"="gene")) %>%
        distinct(DRUG_ID, term, gene_set_size) %>%
        drop_na()
    
    ## add hits of every drug
    x_real = drug_pathways %>%
        left_join(X %>%
            filter(is_sel & !is_target) %>%
            distinct(DRUG_ID, GENE, is_target), 
            by="DRUG_ID"
        ) %>%
        drop_na() %>%
        left_join(
            ontology %>% mutate(in_pathway=TRUE), 
            by=c("term","gene_set_size","GENE"="gene")
        ) %>%
        mutate(in_pathway=replace_na(in_pathway, FALSE))
    
    set.seed(RANDOM_SEED)
    x_null = drug_pathways %>%
        left_join(X %>%
            group_by(DRUG_ID) %>%
            mutate(is_sel = sample(is_sel)) %>%
            ungroup() %>%
            filter(is_sel & !is_target) %>%
            distinct(DRUG_ID, GENE, is_target), 
            by="DRUG_ID"
        ) %>%
        drop_na() %>%
        left_join(
            ontology %>% mutate(in_pathway=TRUE), 
            by=c("term","gene_set_size","GENE"="gene")
        ) %>%
        mutate(in_pathway=replace_na(in_pathway, FALSE))
    
    ## for every drug, how many hits are in the pathway?
    ## for each drug, consider only the pathway with the maximum number of hits
    pathways_oi_real = x_real %>% 
        count(DRUG_ID, term, in_pathway) %>% 
        complete(DRUG_ID, term, in_pathway, fill = list(n=0)) %>% 
        filter(in_pathway) %>%
        left_join(
            ontology %>% distinct(term, gene_set_size), by="term"
        ) %>% 
        group_by(DRUG_ID) %>% 
        slice_max(n, n=1) %>%
        slice_max(gene_set_size, n=1) %>%
        ungroup() %>%
        mutate(dataset="Real")
    
    ## create random background of finding this number of hits per drug and pathway
    pathways_oi_null = x_null %>% 
        count(DRUG_ID, term, in_pathway) %>% 
        complete(DRUG_ID, term, in_pathway, fill = list(n=0)) %>% 
        filter(in_pathway) %>%
        left_join(
            ontology %>% distinct(term, gene_set_size), by="term"
        ) 
    
    pathways_oi_null = pathways_oi_real %>% 
        distinct(DRUG_ID, term) %>% 
        left_join(pathways_oi_null, by=c("DRUG_ID","term")) %>%
        mutate(dataset="Random")
    
    eval_reactome = pathways_oi_null %>%
        bind_rows(pathways_oi_real) %>%
        left_join(
            drug_targets %>% count(DRUG_ID, name="n_targets"),
            by="DRUG_ID"
        )
    
    return(eval_reactome)
}


plot_eda_associations = function(models, drug_screen){
    top_n = 25
    X = models %>% 
        left_join(drug_screen %>% distinct(DRUG_ID, DRUG_NAME) %>% drop_na(), by="DRUG_ID")
    
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
        mutate(is_target=replace_na(is_target, FALSE),
               is_sel = lr_padj<THRESH_FDR & n_obs>THRESH_NOBS)
    
    plts = list()
    
    # - p-value and betas of drugs with target (sign. associated or not)
    plts[["targets-violin-lr_pvalue"]] = X %>% 
        ggviolin(x="is_target", y="lr_pvalue", trim=TRUE, fill="is_target", color=NA, palette=PAL_IS_TARGET) + 
        geom_boxplot(width=0.1, outlier.size=0.1) + 
        guides(fill="none") + 
        stat_compare_means(method="wilcox", size=FONT_SIZE, family=FONT_FAMILY) + 
        labs(x="Is Target", y="LR p-value")
    
    plts[["targets-violin-lr_padj"]] = X %>% 
        ggviolin(x="is_target", y="lr_padj", trim=TRUE, fill="is_target", color=NA, palette=PAL_IS_TARGET) + 
        geom_boxplot(width=0.1, outlier.size=0.1) + 
        guides(fill="none") + 
        stat_compare_means(method="wilcox", size=FONT_SIZE, family=FONT_FAMILY) + 
        labs(x="Is Target", y="LR FDR")
    
    plts[["targets-violin-spldep_coef"]] = X %>% 
        ggviolin(x="is_target", y="spldep_coefficient", trim=TRUE, color=NA, 
                 fill="is_target", palette=PAL_IS_TARGET) + 
        geom_boxplot(width=0.1, outlier.size=0.1) + 
        guides(fill="none") + 
        stat_compare_means(method="wilcox", size=FONT_SIZE, family=FONT_FAMILY) + 
        labs(x="Is Target", y="Spl. Dep. Coefficient")
    
    plts[["targets-violin-spldep_coef_filt"]] = X %>% 
        filter(n_obs > THRESH_NOBS) %>%
        ggviolin(x="is_target", y="spldep_coefficient", trim=TRUE, color=NA, 
                 fill="is_target", palette=PAL_IS_TARGET) + 
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
               path_lab = ifelse(shortest_path_length>=5,">=5",path_lab),
               path_lab = ifelse(shortest_path_length==1e6,"N.R.",path_lab),
               path_lab = factor(path_lab, levels=c("0","1","2","3","4",">=5","N.R."))) %>%
        filter(path_lab != "N.R.") %>% # drop those for which we cannot measure a shortest path
        drop_na(is_sel)
    
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


plot_net_drug_target = function(x, drug_oi){
    
    # make networks
    events_sel = x %>% 
        group_by(DRUG_NAME) %>% 
        slice_max(spldep_coefficient, n=10) %>% 
        ungroup() %>%
        mutate(node=GENE) %>%
        filter(DRUG_NAME == drug_oi)
    genes_sel = events_sel %>% pull(GENE) %>% unique()
    targets_sel = events_sel %>% filter(is_target) %>% pull(GENE) %>% unique()
    
    # make big network
    ppi_net = graph_from_data_frame(ppi, directed=FALSE)
    
    # prepare edges
    ## find possible shortest paths
    notfound = setdiff(genes_sel, names(V(ppi_net)))
    print(sprintf("Genes not found in stringdb: %s", paste0(notfound, collapse=", ")))
    
    r = lapply(intersect(genes_sel, names(V(ppi_net))), function(from){
        r = lapply(targets_sel, function(to){
            path_edges = shortest_paths(ppi_net, from, to, output="both")[["epath"]][[1]]
            edges_oi = subgraph.edges(ppi_net, path_edges) %>% as_data_frame()
            return(edges_oi)
        })
        r = do.call(rbind, r)
        return(r)
    })
    edges = do.call(rbind, r) %>% distinct()
    
    ## prepare nodes
    nodes = data.frame(
            node = edges[,c("from","to")] %>% unlist() %>% unique()
        ) %>%
        left_join(events_sel %>% distinct(GENE, is_target), by=c("node"="GENE")) %>%
        mutate(is_target = replace_na(is_target, FALSE))
    
    net = graph_from_data_frame(edges, directed=FALSE, nodes)
    
    pal_colors = PAL_IS_TARGET
    set.seed(1234)
    lay = net %>% create_layout(layout = "stress")
    
    plt = lay %>% 
        ggraph() + 
        geom_edge_link(color="darkgray", alpha=0.5, width=0.2) + 
        geom_node_point(aes(color=is_target, shape=is_target), alpha=0.8) + 
        geom_node_text(aes(label=name), repel=TRUE, family="Arial", size=2, 
                       segment.size=0.1, box.padding=0.1, force=1) +
        scale_color_manual(values=pal_colors) +
        theme_void() +
        labs(color="Is Target", shape="Is Target", title=drug_oi) +
        theme(
            aspect.ratio = 1,
            plot.title = element_text(hjust = 0.5)
        )
    
    return(plt)
}


plot_gene_structure = function(gene, events, event_info, gene_info, prot_imp){
    coords_gene = gene_info %>%
        filter(Gene %in% gene)
    coords_events = event_info %>%
        filter(EVENT %in% events) %>%
        mutate(
            chr = gsub(":.*", "", COORD_o),
            coords = gsub(".*:", "", COORD_o),
            Start = as.integer(gsub("-.*", "", coords)),
            End = as.integer(gsub(".*-","", coords)),
            event_gene = paste0(EVENT,"_",GENE)
        ) %>%
        left_join(prot_imp %>% mutate(term_clean=replace_na(term_clean,"NA")), by="EVENT")
    chromosome = coords_events %>% pull(chr) %>% unique()
    #impacts = coords_events %>% pull(term_clean) %>% unique()
    
    plt = coords_gene %>%
        ggplot() +
        ggchicklet:::geom_rrect(
            aes(xmin=Start, xmax=End, ymin=-0.3, ymax=0.3),
            fill="gray90",
            color="black",
            radius=unit(1,"pt"),
            size=0.2
        ) +
        geom_rect(
            aes(xmin=Start, xmax=End, ymin=-1, ymax=1, fill=term_clean, color=term_clean),
            coords_events,
        ) +
        geom_text_repel(
            aes(x=Start, y=1, label=EVENT, color=term_clean),
            coords_events,
            ylim=c(2,10),
            segment.size=0.1,
            force=100,
            arrow=arrow(length=unit(1, "mm"), type="closed", ),
            size=FONT_SIZE,
            family=FONT_FAMILY
        ) +
        fill_palette(palette=PAL_PROT_IMP) +
        color_palette(palette=PAL_PROT_IMP) +
        ylim(-6,6) +
        scale_x_continuous(labels=scales::label_comma()) +
        guides(color="none") +
        labs(x=chromosome, title=gene, fill="Protein Impact") +
        theme_void() +
        theme(
            aspect.ratio = 0.3,
            plot.title = element_text(hjust=0.5, family=FONT_FAMILY, size=8),
            legend.position = "top",
            legend.text = element_text(family=FONT_FAMILY, size=6),
            legend.title = element_text(family=FONT_FAMILY, size=6),
            axis.line.x = element_line(size=0.5),
            axis.ticks.x = element_line(),
            axis.ticks.length.x = unit(0.1, "cm"),
            axis.text.x = element_text(angle=45, hjust=1, vjust=1, family=FONT_FAMILY, size=6),
            axis.title.x = element_text(family=FONT_FAMILY, size=8)
        )
    
    return(plt)
}


plot_mediators = function(spldep_models, models, shortest_paths_simple, drug_targets){
    
    X = models %>%
        left_join(shortest_paths_simple, by=c("DRUG_ID","GENE"="source")) %>%
        left_join(drug_targets %>% distinct(DRUG_ID,DRUG_NAME), by="DRUG_ID") %>%
        mutate(is_sel = lr_padj < THRESH_FDR & n_obs > THRESH_NOBS,
               is_sel = replace_na(is_sel, FALSE),
               path_lab = as.character(shortest_path_length),
               path_lab = ifelse(shortest_path_length>=5,">=5",path_lab),
               path_lab = ifelse(shortest_path_length==1e6,"N.R.",path_lab),
               path_lab = replace_na(path_lab,"NA"),
               path_lab = factor(path_lab, levels=c("0","1","2","3","4",">=5","N.R.","NA")),
               is_target = shortest_path_length=="0",
               is_target = replace_na(is_target,FALSE),
               DRUG_NAME = ifelse(DRUG_NAME=="AFATINIB", paste0(DRUG_NAME,"_",DRUG_ID), DRUG_NAME))
    
    plts = list()
    
    # - how many events are we considering?
    plts[["mediators-bar-n_selected"]] = X %>% 
        count(is_sel) %>% 
        mutate(lab= is_sel ) %>%
        ggbarplot(x="lab", y="n", fill="lab", color=NA, 
                  palette=get_palette(PAL_DUAL, 2), 
                  label=TRUE, lab.size=FONT_SIZE, lab.family=FONT_FAMILY) +
        guides(fill="none") +
        yscale("log10", .format=TRUE) +
        labs(x="Larger Splicing Effect & Selected", y="Count")
    
    # what is the best significant association with larger PSI effect?    
    # - mediator on target
    assocs_oi = X %>%
        filter(is_sel & is_target) %>%
        distinct(EVENT,DRUG_ID)
    
    x = assocs_oi %>% 
        left_join(spldep_models %>% 
            distinct(EVENT, event_coefficient_mean), by="EVENT") %>% 
        left_join(X, by=c("EVENT","DRUG_ID")) %>%
        mutate(event_gene = as.factor(event_gene),
               DRUG_ID = paste0(DRUG_ID),
               name = reorder_within(event_gene, spldep_coefficient, DRUG_ID))
    
    # color axis labels based on protein impact
    plts[["mediators-on_target-spldep_coef"]] = x %>% 
        mutate(lab = sprintf("%s | %s", round(spldep_coefficient, 3), DRUG_NAME)) %>%
        ggbarplot(x="name", y="spldep_coefficient", 
              color="drug_screen", fill="is_target", palette=PAL_SINGLE_LIGHT, 
              position=position_dodge(0.8)) +
        color_palette(c("black","white")) + 
        geom_text(aes(y=1.5, label=lab, group=drug_screen), hjust=0.2,
                  position=position_dodge(0.8),
                  size=FONT_SIZE, family=FONT_FAMILY) +
        geom_hline(yintercept=0, linetype="dashed", size=LINE_SIZE) +
        coord_flip(clip="off", ylim=c(-0.55, 1.5)) +
        scale_x_reordered() +
        labs(x="Event & Gene" , y="Spl. Dep. Coefficient", 
             color="Drug Screen", fill="Is Target") +
        new_scale_color() +
        geom_text(aes(y=-1.8,label=event_gene, group=drug_screen, color=term_clean), 
                  position=position_dodge(0), alpha=1,
                  size=FONT_SIZE, family=FONT_FAMILY) +
        color_palette(name="Protein Impact", palette=PAL_PROT_IMP) +
        guides(color=guide_legend(nrow=2)) +
        theme(
            axis.text.y = element_blank(),
            axis.title.y = element_blank()
        )
    
    plts[["mediators-on_target-event_coef"]] = x %>% 
        distinct(event_gene, event_coefficient_mean, is_target, term_clean) %>%
        mutate(lab = round(event_coefficient_mean, 3)) %>%
        ggbarplot(x="event_gene", y="event_coefficient_mean", 
                  color=NA, fill="is_target", palette=PAL_SINGLE_LIGHT,
                  position=position_dodge(0.7)) +
        geom_text(aes(y=sign(event_coefficient_mean)*(abs(event_coefficient_mean)+0.015), 
                      label=lab),
                  position=position_dodge(0.7),
                  size=FONT_SIZE, family=FONT_FAMILY) + 
        geom_hline(yintercept=0, linetype="dashed", size=LINE_SIZE) +
        scale_x_reordered() +
        labs(x="Event & Gene" , y="PSI Coefficient", fill="Is Target") +
        coord_flip(clip="off", ylim=c(-0.17, 0.07)) +
        new_scale_color() +
        geom_text(aes(y=-0.32, label=event_gene, color=term_clean), 
                  size=FONT_SIZE, family=FONT_FAMILY) +
        color_palette(name="Protein Impact", palette=PAL_PROT_IMP) +
        guides(color=guide_legend(nrow=2)) +
        theme(
            axis.text.y = element_blank(),
            axis.title.y = element_blank()
        )
    
    # - mediator up/downstream (off target)
    x = assocs_oi %>%
        distinct(DRUG_ID) %>%
        left_join(X, by="DRUG_ID") %>% 
        left_join(spldep_models %>% 
            distinct(EVENT, event_coefficient_mean, pearson_correlation_mean), by="EVENT") %>%
        filter(is_sel) %>% 
        distinct(DRUG_NAME, EVENT, GENE, event_gene, lr_padj, spldep_coefficient, pearson_correlation_mean,
                 event_coefficient_mean, drug_screen, is_target, path_lab, term_clean)
    
    plts[["mediators-on_target-top_assocs_spldep_all"]] = x %>% 
        group_by(DRUG_NAME) %>% 
        slice_max(abs(spldep_coefficient), n=10) %>% 
        ungroup() %>%
        mutate(DRUG_NAME = as.factor(DRUG_NAME),
               name = reorder_within(event_gene, spldep_coefficient, DRUG_NAME)) %>%
        ggbarplot(x="name", y="spldep_coefficient", 
                  color="drug_screen", fill="is_target", palette=PAL_IS_TARGET, 
                  position=position_dodge(0.7)) + 
        geom_text(aes(y=(abs(spldep_coefficient)+0.08)*sign(spldep_coefficient),
                      label=path_lab, group=drug_screen), 
                  position=position_dodge(0.7),
                  size=FONT_SIZE, family=FONT_FAMILY) + 
        color_palette(c("black","white")) + 
        facet_wrap(~DRUG_NAME, scales="free") + 
        scale_x_reordered() +
        coord_flip() +
        labs(x="Event & Gene", y="Spl. Dep. Coef.", fill="Is Target", color="Drug Screen") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    plts[["mediators-on_target-top_assocs_event_coef_all"]] = x %>% 
        group_by(DRUG_NAME) %>% 
        slice_max(abs(spldep_coefficient), n=10) %>% 
        ungroup() %>%
        mutate(DRUG_NAME = as.factor(DRUG_NAME),
               name = reorder_within(event_gene, spldep_coefficient, DRUG_NAME)) %>%
        ggbarplot(x="name", y="event_coefficient_mean", 
                  color="drug_screen", fill="is_target", palette=PAL_IS_TARGET, 
                  position=position_dodge(0.7)) + 
        geom_text(aes(y=(abs(event_coefficient_mean)+0.01)*sign(event_coefficient_mean), 
                      label=path_lab, group=drug_screen), 
                  position=position_dodge(0.7),
                  size=FONT_SIZE, family=FONT_FAMILY) + 
        color_palette(c("black","white")) + 
        facet_wrap(~DRUG_NAME, scales="free") + 
        scale_x_reordered() +
        coord_flip() +
        labs(x="Event & Gene", y="PSI Coefficient", fill="Is Target", color="Drug Screen") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    plts[["mediators-on_target-top_assocs_spldep_pos"]] = x %>% 
        group_by(DRUG_NAME) %>% 
        slice_max(spldep_coefficient, n=10) %>% 
        ungroup() %>%
        mutate(DRUG_NAME = as.factor(DRUG_NAME),
               name = reorder_within(event_gene, spldep_coefficient, DRUG_NAME)) %>%
        ggbarplot(x="name", y="spldep_coefficient", 
                  color="drug_screen", fill="is_target", palette=PAL_IS_TARGET, 
                  position=position_dodge(0.7)) + 
        geom_text(aes(y=(abs(spldep_coefficient)+0.08)*sign(spldep_coefficient),
                      label=path_lab, group=drug_screen), 
                  position=position_dodge(0.7),
                  size=FONT_SIZE, family=FONT_FAMILY) + 
        color_palette(c("black","white")) + 
        facet_wrap(~DRUG_NAME, scales="free") + 
        scale_x_reordered() +
        coord_flip() +
        labs(x="Event & Gene", y="Spl. Dep. Coef.",, fill="Is Target", color="Drug Screen") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    plts[["mediators-on_target-top_assocs_event_coef_pos"]] = x %>% 
        group_by(DRUG_NAME) %>% 
        slice_max(spldep_coefficient, n=10) %>% 
        ungroup() %>%
        mutate(DRUG_NAME = as.factor(DRUG_NAME),
               name = reorder_within(event_gene, spldep_coefficient, DRUG_NAME)) %>%
        ggbarplot(x="name", y="event_coefficient_mean", 
                  color="drug_screen", fill="is_target", palette=PAL_IS_TARGET, 
                  position=position_dodge(0.7)) + 
        geom_text(aes(y=(abs(event_coefficient_mean)+0.01)*sign(event_coefficient_mean), 
                      label=path_lab, group=drug_screen), 
                  position=position_dodge(0.7),
                  size=FONT_SIZE, family=FONT_FAMILY) + 
        color_palette(c("black","white")) + 
        facet_wrap(~DRUG_NAME, scales="free") + 
        scale_x_reordered() +
        coord_flip() +
        labs(x="Event & Gene", y="PSI Coefficient", fill="Is Target", color="Drug Screen") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    plts[["mediators-on_target-top_assocs_spldep_neg"]] = x %>% 
        group_by(DRUG_NAME) %>% 
        slice_min(spldep_coefficient, n=10) %>% 
        ungroup() %>%
        mutate(DRUG_NAME = as.factor(DRUG_NAME),
               name = reorder_within(event_gene, spldep_coefficient, DRUG_NAME)) %>%
        ggbarplot(x="name", y="spldep_coefficient", 
                  color="drug_screen", fill="is_target", palette=PAL_IS_TARGET, 
                  position=position_dodge(0.7)) + 
        geom_text(aes(y=(abs(spldep_coefficient)+0.08)*sign(spldep_coefficient),
                      label=path_lab, group=drug_screen), 
                  position=position_dodge(0.7),
                  size=FONT_SIZE, family=FONT_FAMILY) + 
        color_palette(c("black","white")) + 
        facet_wrap(~DRUG_NAME, scales="free") + 
        scale_x_reordered() +
        coord_flip() +
        labs(x="Event & Gene", y="Spl. Dep. Coef.",, fill="Is Target", color="Drug Screen") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    plts[["mediators-on_target-top_assocs_event_coef_neg"]] = x %>% 
        group_by(DRUG_NAME) %>% 
        slice_min(spldep_coefficient, n=10) %>% 
        ungroup() %>%
        mutate(DRUG_NAME = as.factor(DRUG_NAME),
               name = reorder_within(event_gene, spldep_coefficient, DRUG_NAME)) %>%
        ggbarplot(x="name", y="event_coefficient_mean", 
                  color="drug_screen", fill="is_target", palette=PAL_IS_TARGET, 
                  position=position_dodge(0.7)) + 
        geom_text(aes(y=(abs(event_coefficient_mean)+0.01)*sign(event_coefficient_mean), 
                      label=path_lab, group=drug_screen), 
                  position=position_dodge(0.7),
                  size=FONT_SIZE, family=FONT_FAMILY) + 
        color_palette(c("black","white")) + 
        facet_wrap(~DRUG_NAME, scales="free") + 
        scale_x_reordered() +
        coord_flip() +
        labs(x="Event & Gene", y="PSI Coefficient", fill="Is Target", color="Drug Screen") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))

    ## PPI MOA
    plts[["mediators-drugs_oi-top_assocs_ppi_pos"]] = list(
            "NUTLIN-3A" = plot_net_drug_target(x, "NUTLIN-3A (-)"),
            "AZD4547" = plot_net_drug_target(x, "AZD4547")
        ) %>% ggarrange(plotlist = ., common.legend=TRUE)
    
    ## gene structures
    genes_sel = x %>% 
        group_by(DRUG_NAME) %>% 
        slice_max(spldep_coefficient, n=10) %>% 
        ungroup() %>%
        filter(DRUG_NAME=="NUTLIN-3A (-)") %>%
        pull(GENE) %>% unique()
    plts[["mediators-drugs_oi-gene_structs"]] = lapply(genes_sel, function(gene_sel){
        print(gene_sel)
        events_sel = x %>% 
            group_by(DRUG_NAME) %>% 
            slice_max(spldep_coefficient, n=10) %>% 
            ungroup() %>%
            filter(DRUG_NAME=="NUTLIN-3A (-)" & GENE==gene_sel) %>%
            pull(EVENT) %>% unique()
        plt = plot_gene_structure(gene_sel, events_sel, event_info, gene_info, ontologies[["protein_impact"]])
        return(plt)
    }) %>% ggarrange(plotlist=., common.legend=TRUE)
    
    
    ## co-splicing
    drugs_oi = c("NUTLIN-3A (-)","AZD4547")
    events_oi = x %>% 
        filter(DRUG_NAME %in% drugs_oi) %>%
        group_by(DRUG_NAME) %>% 
        slice_max(spldep_coefficient, n=10) %>% 
        ungroup() %>%
        distinct(EVENT, event_gene, DRUG_NAME, is_target, spldep_coefficient)
    
    corrs = sapply(drugs_oi, function(drug_oi){
        corr = events_oi %>% 
            filter(DRUG_NAME %in% drug_oi) %>%
            left_join(splicing, by="EVENT") %>%
            select(-one_of(c("EVENT","DRUG_NAME"))) %>%
            column_to_rownames("event_gene") %>%
            as.matrix() %>%
            t() %>%
            cor(method="spearman",use="pairwise.complete.obs")
        return(corr)
    }, simplify=FALSE)
    
    mat = corrs[["NUTLIN-3A (-)"]]
    plts[["mediators-drugs_oi-cosplicing-nutlin"]] = Heatmap(
        mat, name="Spear. Corr.", 
        row_names_gp = gpar(fontsize=6, fontfamily=FONT_FAMILY),
        column_names_gp = gpar(fontsize=6, fontfamily=FONT_FAMILY),
        heatmap_legend_param = list(legend_gp = gpar(fontsize = 6, fontfamily=FONT_FAMILY)),
        cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 6))
        }) %>% 
        draw() %>%
        grid.grabExpr() %>%
        as.ggplot()
    
    mat = corrs[["AZD4547"]]
    plts[["mediators-drugs_oi-cosplicing-AZD4547"]] = Heatmap(
        mat, name="Spear. Corr.", 
        row_names_gp = gpar(fontsize=6, fontfamily=FONT_FAMILY),
        column_names_gp = gpar(fontsize=6, fontfamily=FONT_FAMILY),
        heatmap_legend_param = list(legend_gp = gpar(fontsize = 6, fontfamily=FONT_FAMILY)),
        cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 6))
        }) %>% 
        draw() %>%
        grid.grabExpr() %>%
        as.ggplot()
    
    
    plts[["mediators-drugs_oi-box-splicing"]] = events_oi %>% 
        filter(DRUG_NAME %in% drugs_oi) %>%
        left_join(splicing, by="EVENT") %>%
        pivot_longer(-one_of(c("EVENT","DRUG_NAME","event_gene","is_target","spldep_coefficient")), 
                     names_to = "sampleID", values_to="psi") %>%
        drop_na(psi) %>%
        mutate(DRUG_NAME = as.factor(DRUG_NAME),
               name = reorder_within(event_gene, spldep_coefficient, DRUG_NAME)) %>%
        ggplot(aes(x=name, y=psi)) +
        geom_boxplot(aes(color=is_target), width=0.8, outlier.size=0.1) +
        color_palette(palette=PAL_IS_TARGET) +
        facet_wrap(~DRUG_NAME, ncol=1, scales="free_y") +
        labs(x="Event & Gene", y="PSI", color="Is Target") +
        scale_x_reordered() +
        theme_pubr() +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        coord_flip()
    
    # - examples exon-drug synergy
    x = X %>%
        drop_na(spldep_coefficient) %>%
        filter(is_sel)
    
    ## associations with drugs
    drugs_oi = c("VENETOCLAX")
    
    plts[["mediators-drugs_oi-synergy"]] = x %>% 
        filter(DRUG_NAME %in% drugs_oi) %>% 
        group_by() %>%
        slice_max(abs(spldep_coefficient), n=10) %>%
        ungroup() %>%
        mutate(DRUG_NAME = as.factor(DRUG_NAME),
               name = reorder_within(event_gene, spldep_coefficient, DRUG_NAME)) %>%
        ggbarplot(x="name", y="spldep_coefficient", 
                  color="drug_screen", fill="is_target", palette=PAL_DUAL, 
                  position=position_dodge(0.7)) + 
        geom_text(aes(y=sign(spldep_coefficient)*(abs(spldep_coefficient)+0.1), label=path_lab, group=drug_screen), 
                  position=position_dodge(0.7),
                  size=FONT_SIZE, family=FONT_FAMILY) + 
        color_palette(c("black","white")) + 
        facet_wrap(~DRUG_NAME, scales="free") + 
        scale_x_reordered() +
        labs(x="Drug", y="Spl. Dep. Coefficient", fill="Is Target", color="Drug Screen") +
        coord_flip() +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    ## associations with events
    events_oi = c('HsaEX0034998','HsaEX0038400') # KRAS
    # MDM2
    
    plts[["mediators-events_oi-synergy"]] = x %>% 
        filter(EVENT %in% events_oi) %>% 
        mutate(event_gene = as.factor(event_gene),
               name = reorder_within(DRUG_NAME, spldep_coefficient, event_gene)) %>%
        ggbarplot(x="name", y="spldep_coefficient", 
                  color="drug_screen", fill="is_target", palette=PAL_DUAL, 
                  position=position_dodge(0.7)) + 
        geom_text(aes(y=sign(spldep_coefficient)*(abs(spldep_coefficient)+0.1), label=path_lab, group=drug_screen), 
                  position=position_dodge(0.7),
                  size=FONT_SIZE, family=FONT_FAMILY) + 
        color_palette(c("black","white")) + 
        facet_wrap(~event_gene, scales="free") + 
        scale_x_reordered() +
        labs(x="Drug", y="Spl. Dep. Coefficient", fill="Is Target", color="Drug Screen") +
        coord_flip() +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    return(plts)
}


plot_drug_rec = function(estimated_response, drug_screen, models){
    
    # all models
    models_oi = models %>%
        filter(lr_padj<THRESH_FDR  & n_obs>THRESH_NOBS) %>% # we don't filter in the estimation
        distinct(ID, pearson_correlation)
    
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
        fill_palette(palette=PAL_IN_DEMETER2) + 
        labs(x="Drug Screen", y="Pearson Correlation", fill="In Demeter2") +
        geom_text(aes(y=1.00, label=n), 
                  corrs %>% count(in_demeter2, drug_screen), 
                  position=position_dodge(0.9),
                  size=FONT_SIZE, family=FONT_FAMILY) +
        theme_pubr()
    
    plts[["drug_rec-npos_vs_corrs"]] = corrs %>% 
        ggscatter(x="n_pos", y="correlation", alpha=0.5, size=1,
                  color="in_demeter2", palette=PAL_IN_DEMETER2) + 
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
        stat_cor(method="pearson", size=FONT_SIZE, family=FONT_FAMILY) + 
        labs(x="Real log(IC50)", y="Predicted log(IC50)") + 
        theme_pubr(border=TRUE) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    return(plts)
}


plot_reactome = function(eval_reactome){
    
    X = eval_reactome
    
    plts = list()
    
    plts[["reactome-violin-hits_in_pathway"]] = X %>% 
        ggviolin(x="dataset", y="n", trim=TRUE, 
                 fill="dataset", color=NA, palette=c("gray","#6CA29B")) + 
        geom_boxplot(width=0.1, outlier.size=0.1) +
        stat_compare_means(method="wilcox.test", family=FONT_FAMILY, size=FONT_SIZE) +
        guides(fill="none") + 
        labs(x="Drug-Exon Association", y="N. Associations in Drug Pathway")
    
    return(plts)
}

prep_examples = function(splicing, genexpr, snv, ontology, gene_info){
## get splicing of HsaEX0038414_MDM4 across samples
    splicing_oi = splicing %>%
        filter(EVENT == "HsaEX0038414") %>%
        pivot_longer(-EVENT, names_to="DepMap_ID", values_to="psi")
    
    ## log-transform TPMs and get genes in P53 pathway only
    tpm = genexpr %>% 
        mutate_at(vars(-ID), function(x){log2(x+1)})
    
    ## get samples with TP53 mutations
    muts_oi = snv %>%
        filter(Hugo_Symbol=="TP53" & Variant_Classification!="Silent") %>%
        group_by(DepMap_ID) %>%
        summarize(mutated_tp53 = TRUE) %>%
        ungroup()
    
    ## combine info
    genes_oi = ontology %>% filter(term == "HALLMARK_P53_PATHWAY") %>% pull(gene)
    genexpr_pathway = gene_info %>% 
        filter(Gene %in% genes_oi) %>% 
        distinct(EnsID,Gene) %>%
        left_join(tpm, by=c("EnsID"="ID")) %>%
        pivot_longer(-c(Gene,EnsID), names_to="DepMap_ID", values_to="tpm") %>%
        left_join(splicing_oi, by="DepMap_ID") %>%
        mutate(psi_bin = cut(psi, breaks=c(0,25,50,75,100))) %>%
        drop_na(psi_bin) %>%
        left_join(muts_oi, by="DepMap_ID") %>% 
        mutate(
            mutated_tp53 = replace_na(mutated_tp53, FALSE),
            mutated_tp53 = ifelse(mutated_tp53, "TP53mut", "WT")
        )
    
    ctl = genexpr_pathway %>%
        filter(psi_bin == "(75,100]") %>%
        group_by(Gene) %>%
        summarize(tpm_ctl = median(tpm)) %>%
        ungroup()
    
    X = genexpr_pathway %>%
        left_join(ctl) %>%
        mutate(tpm_fc = tpm - tpm_ctl)
    
    return(X)
}


plot_examples = function(examples){
    plts = list()
    
    # Does the activation of the P53 pathway change with the splicing of HsaEX0038414_MDM4?
    X = examples
    
    ## distribution of event inclusion in each sample
    plts[["examples-mdm4-psi-distr"]] = X %>%
        distinct(DepMap_ID, psi) %>%
        gghistogram(x="psi", fill=PAL_IS_TARGET[1], color=NA) +
        geom_vline(xintercept=c(25,50,75), linetype="dashed", size=LINE_SIZE) +
        labs(x="PSI HsaEX0038414_MDM4", y="Count")
    
    plts[["examples-mdm4-psi-bins"]] = X %>%
        distinct(DepMap_ID, psi_bin) %>%
        count(psi_bin) %>%
        drop_na(psi_bin) %>%
        ggbarplot(x="psi_bin", y="n", fill=PAL_IS_TARGET[1], color=NA,
                  label=TRUE, lab.size=FONT_SIZE, lab.family=FONT_FAMILY) +
        labs(x="PSI HsaEX0038414_MDM4", y="Count") +
        theme(aspect.ratio=1)
    
    ## distributions of gene expression
#     comparisons = list(
#         c("(0,25]","(75,100]"),
#         c("(25,50]","(75,100]"),
#         c("(50,75]","(75,100]")
#     )
    plts[["examples-mdm4-bins_vs_genexpr_fc_tp53"]] = X %>%
        filter(Gene=="TP53") %>%
        ggviolin(x="psi_bin", y="tpm", fill=PAL_SINGLE_DARK, color=NA, trim=TRUE) +
        geom_boxplot(width=0.1, outlier.size=0.1) + 
        geom_text(aes(label=n, y=-0.5), . %>% count(psi_bin,mutated_tp53) %>% mutate(n=sprintf("n=%s",n)), 
                  size=FONT_SIZE, family=FONT_FAMILY) +
        stat_compare_means(method="wilcox.test", ref.group="(0,25]", label="p.signif",
                           size=FONT_SIZE, family=FONT_FAMILY) +
        facet_wrap(~mutated_tp53) +
        labs(x="PSI HsaEX0038414_MDM4", y="log2(TPM+1) TP53") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY), aspect.ratio=1)
    
    return(plts)
}


make_plots = function(
    models, drug_screen,
    drug_targets, shortest_paths, 
    spldep_models, shortest_paths_simple,
    estimated_response, 
    eval_reactome, 
    examples
){
    plts = list(
        plot_eda_associations(models, drug_screen),
        plot_targets(models, drug_targets),
        plot_ppi(models, shortest_paths),
        plot_mediators(spldep_models, models, shortest_paths_simple, drug_targets),
        plot_drug_rec(estimated_response, drug_screen, models),
        plot_reactome(eval_reactome),
        plot_examples(examples)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(
    models, drug_screen,
    drug_targets, shortest_paths, 
    spldep_models, shortest_paths_simple,
    estimated_response, 
    eval_reactome, 
    examples
){
    
    figdata = list(
        "drug_event_assoc" = list(
            "drug_sensitivities" = drug_screen %>% mutate(in_training_set=in_demeter2),
            "drug_targets" = drug_targets,
            "model_summaries" = models,
            "evaluation_reactome" = eval_reactome,
            "shortest_paths" = shortest_paths,
            "shortest_paths_simple" = shortest_paths_simple,
            "pred_ic50_by_drug" = estimated_response,
            "example_mdm4" = examples %>% filter(Gene=="TP53")
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
    save_plt(plts, "associations-top_drug_counts", ".pdf", figs_dir, width=6, height=8)
    save_plt(plts, "associations-spldep_counts", ".pdf", figs_dir, width=5, height=5)
    save_plt(plts, "associations-top_spldep_counts", ".pdf", figs_dir, width=6, height=8)
    save_plt(plts, "associations-agreement_between-upset", ".pdf", figs_dir, width=7, height=6)
    save_plt(plts, "associations-agreement_between-correlation", ".pdf", figs_dir, width=5, height=5)
    save_plt(plts, "associations-overlaps_screens-upset", ".pdf", figs_dir, width=8, height=7)
    
    # targets
    save_plt(plts, "targets-violin-lr_pvalue", ".pdf", figs_dir, width=3.5, height=3)
    save_plt(plts, "targets-violin-lr_padj", ".pdf", figs_dir, width=5, height=5)
    save_plt(plts, "targets-bar-n_signifs", ".pdf", figs_dir, width=3.5, height=3.5)
    
    # ppi
    save_plt(plts, "ppi-bar-counts", ".pdf", figs_dir, width=5, height=6)
    save_plt(plts, "ppi-bar-prop", ".pdf", figs_dir, width=5, height=6)
    save_plt(plts, "ppi-bar-rel_prop", ".pdf", figs_dir, width=4.5, height=5)
    save_plt(plts, "ppi-box-lr_padj", ".pdf", figs_dir, width=8, height=4)
    save_plt(plts, "ppi-box-spldep_coef", ".pdf", figs_dir, width=8, height=6)
    
    # reactome
    save_plt(plts, "reactome-violin-hits_in_pathway", ".pdf", figs_dir, width=3.5, height=3.5)
    
    # mediators
    save_plt(plts, "mediators-bar-n_selected", ".pdf", figs_dir, width=5, height=5)
    save_plt(plts, "mediators-on_target-spldep_coef", ".pdf", figs_dir, width=7, height=9)
    save_plt(plts, "mediators-on_target-event_coef", ".pdf", figs_dir, width=7, height=7)
    save_plt(plts, "mediators-on_target-top_assocs_spldep_all", ".pdf", figs_dir, width=15, height=13)
    save_plt(plts, "mediators-on_target-top_assocs_spldep_pos", ".pdf", figs_dir, width=15, height=13)
    save_plt(plts, "mediators-on_target-top_assocs_spldep_neg", ".pdf", figs_dir, width=15, height=13)
    save_plt(plts, "mediators-on_target-top_assocs_event_coef_all", ".pdf", figs_dir, width=15, height=13)
    save_plt(plts, "mediators-on_target-top_assocs_event_coef_pos", ".pdf", figs_dir, width=15, height=13)
    save_plt(plts, "mediators-on_target-top_assocs_event_coef_neg", ".pdf", figs_dir, width=15, height=13)
    save_plt(plts, "mediators-drugs_oi-gene_structs", ".pdf", figs_dir, width=10, height=10)
    save_plt(plts, "mediators-drugs_oi-top_assocs_ppi_pos", ".pdf", figs_dir, width=8, height=5)
    save_plt(plts, "mediators-drugs_oi-cosplicing-AZD4547", ".pdf", figs_dir, width=12, height=10)
    save_plt(plts, "mediators-drugs_oi-cosplicing-nutlin", ".pdf", figs_dir, width=12, height=10)
    save_plt(plts, "mediators-drugs_oi-box-splicing", ".pdf", figs_dir, width=5.5, height=9)
    save_plt(plts, "mediators-events_oi-synergy", ".pdf", figs_dir, width=10, height=5)
    save_plt(plts, "mediators-drugs_oi-synergy", ".pdf", figs_dir, width=5, height=5)
    
    # drug recommendations
    save_plt(plts, "drug_rec-n_screens", ".pdf", figs_dir, width=6, height=7)
    save_plt(plts, "drug_rec-correlations", ".pdf", figs_dir, width=5, height=5)
    save_plt(plts, "drug_rec-npos_vs_corrs", ".pdf", figs_dir, width=6, height=9)
    save_plt(plts, "drug_rec-best_worse", ".pdf", figs_dir, width=7, height=8)
    
    # examples
    save_plt(plts, "examples-mdm4-psi-distr", ".pdf", figs_dir, width=4.35, height=1.5)
    save_plt(plts, "examples-mdm4-psi-bins", ".pdf", figs_dir, width=4.5, height=4.5)
    save_plt(plts, "examples-mdm4-bins_vs_genexpr_fc_tp53", ".pdf", figs_dir, width=7.5, height=6)
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
        make_option("--metadata_file", type="character"),
        make_option("--models_file", type="character"),
        make_option("--drug_targets_file", type="character"),
        make_option("--rnai_file", type="character"),
        make_option("--estimated_response_file", type="character"),
        make_option("--drug_screens_dir", type="character"),
        make_option("--spldep_ccle_file", type="character"),
        make_option("--paths_real_file", type="character"),
        make_option("--spldep_models_file", type="character"),
        make_option("--msigdb_dir", type="character"),
        make_option("--protein_impact_file", type="character"),
        make_option("--splicing_file", type="character"),
        make_option("--genexpr_file", type="character"),
        make_option("--snv_file", type="character"),
        make_option("--ppi_file", type="character"),
        make_option("--event_info_file", type="character"),
        make_option("--gene_info_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    metadata_file = args[["metadata_file"]]
    models_file = args[["models_file"]]
    drug_targets_file = args[["drug_targets_file"]]
    rnai_file = args[["rnai_file"]]
    estimated_response_file = args[["estimated_response_file"]]
    drug_screens_dir = args[["drug_screens_dir"]]
    spldep_ccle_file = args[["spldep_ccle_file"]]
    paths_real_file = args[["paths_real_file"]]
    spldep_models_file = args[["spldep_models_file"]]
    msigdb_dir = args[["msigdb_dir"]]
    protein_impact_file = args[["protein_impact_file"]]
    splicing_file = args[["splicing_file"]]
    genexpr_file = args[["genexpr_file"]]
    snv_file = args[["snv_file"]]
    ppi_file = args[["ppi_file"]]
    event_info_file = args[["event_info_file"]]
    gene_info_file = args[["gene_info_file"]]
    figs_dir = args[["figs_dir"]]
    
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
    ontologies = load_ontologies(msigdb_dir, protein_impact_file)
    splicing = read_tsv(splicing_file)
    genexpr = read_tsv(genexpr_file)
    snv = read_csv(snv_file)
    ppi = read_tsv(ppi_file)
    event_info = read_tsv(event_info_file)
    gene_info = read_tsv(gene_info_file)

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
    ## ReactomeDB
    eval_reactome = evaluate_reactome(models, drug_targets, ontologies[["reactome"]])
    
    ## protein impact
    models = models %>%
        left_join(ontologies[["protein_impact"]], by="EVENT")
    
    ## mechanistic examples
    examples = prep_examples(splicing, genexpr, snv, ontologies[["hallmarks"]], gene_info)
    
    # make plots
    plts = make_plots(
            models, drug_screen,
            drug_targets, shortest_paths, 
            spldep_models, shortest_paths_simple,
            estimated_response, 
            eval_reactome, 
            examples
    )
    
    # make figdata
    figdata = make_figdata(
            models, drug_screen,
            drug_targets, shortest_paths, 
            spldep_models, shortest_paths_simple,
            estimated_response, 
            eval_reactome, 
            examples
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