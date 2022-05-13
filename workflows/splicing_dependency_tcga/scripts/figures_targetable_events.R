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
# Targetable events
# - find vulnerable exons that are differentially spliced in tumor samples
# - validate that low-variant exons in CCLE cannot be differentially spliced
# - validate how splicing dependencies vary between tissues in CCLE and TCGA
#
# Effect of tumor microenvironment
# - how ranges of splicing dependencies change in tumors compared with CCLE

require(tidyverse)
require(ggpubr)
require(cowplot)
require(scattermore)
require(ggrepel)
require(clusterProfiler)
require(ComplexHeatmap)
require(gridExtra)
require(ggplotify)
require(extrafont)
require(writexl)

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

# variables
THRESH_FDR = 0.05
THRESH_MEDIAN_DIFF = 5
MIN_SAMPLES = 10

# formatting
PAL_SINGLE_ACCENT = "orange"
PAL_SINGLE_LIGHT = "#6AC2BF"
PAL_SINGLE_DARK = "#716454"
PAL_DUAL = c(PAL_SINGLE_DARK, PAL_SINGLE_LIGHT)
LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

# Development
# -----------
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# RESULTS_DIR = file.path(ROOT,'results','splicing_dependency_tcga')
# diff_result_sample_file = file.path(RESULTS_DIR,'files','PANCAN','mannwhitneyu-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz')
# annotation_file = file.path(RAW_DIR,'VastDB','EVENT_INFO-hg38_noseqs.tsv')
# selected_events_file = file.path(ROOT,'results','model_splicing_dependency','files','selected_models-EX.txt')
# psi_ccle_file = file.path(PREP_DIR,'event_psi','CCLE-EX.tsv.gz')
# spldep_stats_file = file.path(RESULTS_DIR,'files','PANCAN','summary_splicing_dependency-EX.tsv.gz')
# spldep_ccle_file = file.path(ROOT,'results','model_splicing_dependency','files','splicing_dependency-EX','mean.tsv.gz')
# diff_result_subtypes_file = file.path(RESULTS_DIR,'files','PANCAN_subtypes','mannwhitneyu-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz')
# spldep_stats_subtypes_file = file.path(RESULTS_DIR,'files','PANCAN_subtypes','summary_splicing_dependency-EX.tsv.gz')
# spldep_stats_ccle_cancers_file = file.path(ROOT,'results','model_splicing_dependency','files','splicing_dependency_summaries-EX','primary_disease.tsv.gz')
# indices_file = file.path(RESULTS_DIR,'files','PANCAN','correlation_spldep_indices-EX.tsv.gz')
# figs_dir = file.path(RESULTS_DIR,'figures','targetable_events')
# ontology_file = file.path(ROOT,"results","model_splicing_dependency",'files','ENCORE','kd_gene_sets-EX.tsv.gz')
# cancer_events_file = file.path(ROOT,"support","cancer_events.tsv")

# MODELS_DIR = file.path(ROOT,"results","splicing_dependency_drugs")
# drug_models_file = file.path(MODELS_DIR,"files","model_summaries_drug_response-EX.tsv.gz")
# drug_targets_file = file.path(PREP_DIR,'drug_screens','drug_targets.tsv.gz')
# drug_treatments_file = file.path(PREP_DIR,'drug_treatments','PANCAN.tsv.gz')
# metadata_file = file.path(PREP_DIR,'metadata','PANCAN.tsv.gz')
# selected_drugs_file = file.path(MODELS_DIR,'files','selected_models-EX.txt')


##### FUNCTIONS #####
prep_diff_result = function(diff_result, spldep_stats){
    diff_result = diff_result %>%
        rename_all(recode, index = "EVENT") %>%
        mutate(event_type = gsub('Hsa','',gsub("[^a-zA-Z]", "",EVENT))) %>%
        left_join(spldep_stats, by=c('cancer_type','EVENT')) %>%
        drop_na(event_gene) %>%
        group_by(cancer_type) %>%
        mutate(psi__padj = p.adjust(psi__pvalue, 'fdr'),
               psi__log10_padj = -log10(psi__padj)) %>%
        mutate(psi__is_significant = psi__padj<THRESH_FDR & 
                                     abs(psi__median_diff)>THRESH_MEDIAN_DIFF)
    return(diff_result)
}


plot_top_candidates_sample_type = function(diff_result, cancer_events, patt=''){
    
    plts = list()

    # how many samples do we have for each cancer and sample type?
    X = diff_result %>% 
        mutate(PT = `psi__condition_a-n_present` + `psi__condition_a-n_nan`,
               STN = `psi__condition_b-n_present` + `psi__condition_b-n_nan`) %>% 
        distinct(cancer_type,PT,STN) %>%
        pivot_longer(c(PT,STN), names_to='sample_type',values_to='n')
    
    plts[['top_samples-sample_counts_cancer']] = X %>% 
        ggbarplot(x='cancer_type', y='n', fill='sample_type', color=NA,
                  position=position_dodge(0.7), label=TRUE, lab.size=FONT_SIZE, 
                  lab.family=FONT_FAMILY, palette='lancet') + 
        labs(x='Cancer Type', y='N. Samples') + 
        theme_pubr(x.text.angle = 45) + 
        geom_hline(yintercept=10, linetype='dashed', size=LINE_SIZE)
    
    # ranking with differential analyses
    cancers_oi = X %>% 
        group_by(cancer_type) %>% 
        summarize(to_keep=sum(n >= MIN_SAMPLES)==2) %>% 
        filter(to_keep) %>% 
        pull(cancer_type)
    
    X = diff_result %>%
        filter(cancer_type %in% cancers_oi) %>%
        mutate(sign_dpsi = sign(psi__median_diff),
               sign_spldep = sign(mean),
               event_gene = ifelse(EVENT%in%cancer_events[["EVENT"]], 
                                   paste0("*",event_gene), event_gene))
    
    plts[['top_samples-dpsi_vs_spldep-scatter']] = X %>%
        ggplot(aes(x=psi__median_diff, 
                   y=mean, 
                   color=psi__is_significant)) +
        geom_scattermore(
            pixels = c(1000,1000), 
            pointsize = 10) +
        facet_wrap(~cancer_type, ncol=4) +
        theme_pubr(border=TRUE) +
        guides(color="none") +
        labs(x='Delta PSI', y='median(Spl. Dep. in PT)') +
        ggpubr::color_palette(c("black","orange")) +
        geom_hline(yintercept=0, linetype='dashed', size=LINE_SIZE) +
        geom_vline(xintercept=0, linetype='dashed', size=LINE_SIZE) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    plts[['top_samples-dpsi_vs_spldep-selection']] = X %>%
        filter(psi__is_significant) %>%
        mutate(sign_combined = sprintf('sign_dpsi=%s & sign_spldep=%s',sign_dpsi,sign_spldep)) %>%
        count(sign_combined) %>% 
        drop_na() %>%
        ggbarplot(x='cancer_type', y='n', fill='sign_combined', palette='jco', 
                  color=FALSE, label=TRUE, lab.size=FONT_SIZE, lab.family=FONT_FAMILY,
                  position=position_dodge(0.9)) +
        theme_pubr(x.text.angle = 45, legend='right') +
        labs(x='Cancer Type', y='N. Significant Events')
   
    plts[['top_samples-dpsi_vs_spldep-candidates_spldep']] = X %>% 
        filter(cancer_type %in% cancers_oi) %>%
        filter(psi__is_significant & 
               ((sign_dpsi>0 & sign_spldep<0) | (sign_dpsi<0 & sign_spldep>0))) %>%
        arrange(-mean) %>%
        ggbarplot(x='event_gene', y='mean', 
                  fill='cancer_type', position=position_dodge(0.9), color=FALSE, 
                  palette=get_palette('Paired',length(unique(X[['cancer_type']])))) + 
        geom_hline(yintercept=0, linetype='dashed', size=LINE_SIZE) +
        labs(x='Event & Gene', y='mean(Spl. Dep. in PT)', fill='Cancer Type') +
        coord_flip()
    
    plts[['top_samples-dpsi_vs_spldep-candidates_dpsi']] = X %>% 
        filter(cancer_type %in% cancers_oi) %>%
        filter(psi__is_significant & 
               ((sign_dpsi>0 & sign_spldep<0) | (sign_dpsi<0 & sign_spldep>0))) %>%
        arrange(-mean) %>%
        ggbarplot(x='event_gene', y='psi__median_diff', 
                  fill='cancer_type', position=position_dodge(0.9), color=FALSE, 
                  palette=get_palette('Paired',length(unique(X[['cancer_type']])))) + 
        geom_hline(yintercept=c(-5,5), linetype='dashed', size=LINE_SIZE) +
        labs(x='Event & Gene', y='Delta PSI', fill='Cancer Type') +
        coord_flip()
    
    plts[['top_samples-dpsi_vs_spldep-candidates_harm']] = X %>% 
        mutate(harm_score = ifelse(mean<0,
                                   (-1) * mean * (0-`psi__condition_a-median`), # remove
                                   (-1) * mean * (100-`psi__condition_a-median`)) # include
               
               ) %>%
        filter(cancer_type %in% cancers_oi) %>%
        filter(psi__is_significant & 
               ((sign_dpsi>0 & sign_spldep<0) | (sign_dpsi<0 & sign_spldep>0))) %>%
        arrange(-mean) %>%
        ggbarplot(x='event_gene', y='harm_score', 
                  fill='cancer_type', position=position_dodge(0.9), color=FALSE, 
                  palette=get_palette('Paired',length(unique(X[['cancer_type']])))) + 
        geom_hline(yintercept=0, linetype='dashed', size=LINE_SIZE) +
        labs(x='Event & Gene', y='Harm Score', fill='Cancer Type') +
        coord_flip()
    
    x = X %>% 
        filter(cancer_type %in% cancers_oi) %>%
        filter(psi__is_significant & 
               ((sign_dpsi>0 & sign_spldep<0) | (sign_dpsi<0 & sign_spldep>0)))
    a = x %>% dplyr::select(c("EVENT","event_gene","psi__condition_a-median",
                              "psi__condition_a-mad","psi__condition_a"))
    colnames(a) = gsub("_a","",colnames(a))
    b = x %>% dplyr::select(c("EVENT","event_gene","psi__condition_b-median",
                              "psi__condition_b-mad","psi__condition_b"))
    colnames(b) = gsub("_b","",colnames(b))
    x = bind_rows(a,b)
    
    plts[['top_samples-dpsi_vs_spldep-candidates_psi']] = x %>% 
        mutate(ymin=`psi__condition-median` - `psi__condition-mad`,
               ymax=`psi__condition-median` + `psi__condition-mad`) %>%
        ggplot(aes(x=cancer_type, y=`psi__condition-median`)) +
        geom_pointrange(aes(ymin=ymin, ymax=ymax, color=psi__condition), 
                        fatten=1, position=position_dodge(0.5)) + 
        facet_wrap(~event_gene, ncol=5) +
        color_palette("npg") + 
        labs(x="Cancer Type", y="PSI", color="Sample Type") +
        theme_pubr() + 
        coord_flip() +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
        
    # enrichment ENCORE
    events_oi = x %>% pull(EVENT) %>% unique()
    universe = X %>% pull(EVENT) %>% unique()
    
    events_oi %>% length() %>% print()
    events_oi %in% ontology[["EVENT"]] %>% sum() %>% print()
    
    enrichment = enricher(
        events_oi, TERM2GENE=ontology, universe=universe, maxGSSize=Inf
        ) %>%
        as.data.frame() %>%
        rowwise() %>%
        mutate(gene_ratio = eval(parse(text=GeneRatio))) %>%
        ungroup()
    
    plts[["top_samples-enrichment-KD-dot"]] = enrichment %>%
        slice_max(gene_ratio, n=10) %>%
        arrange(gene_ratio) %>%
        ggscatter(x="Description", y="gene_ratio", 
                  size="Count", color="p.adjust") +
        gradient_color(c("blue", "red")) +
        scale_size(range=c(0.5,3)) +
        labs(x="Event Set", y="Event Ratio") +
        coord_flip()
    
    names(plts) = paste0(patt,names(plts))
    
    return(plts)
}


plot_ccle_vs_tcga = function(diff_result_sample_raw, psi_ccle){
    # are there any low-varying events in cell lines that could be differentially spliced?
    cancers_oi = c('BRCA','COAD','HNSC','KICH','KIRC','KIRP','LIHC','LUAD',
                   'LUSC','PRAD','READ','THCA','UCEC')
    common_events = intersect(diff_result_sample_raw[['index']], psi_ccle[['EVENT']])
    psi = psi_ccle %>% filter(EVENT %in% common_events) %>% column_to_rownames('EVENT')
    X = data.frame(ccle_std = apply(psi,1,sd,na.rm=TRUE)) %>% 
        rownames_to_column('index') %>% 
        left_join(diff_result_sample_raw, by='index') %>%
        filter(cancer_type %in% cancers_oi) %>%
        mutate(is_significant = psi__padj<THRESH_FDR & abs(psi__median_diff)>THRESH_MEDIAN_DIFF)
    
    plts = list()
    plts[['ccle_vs_tcga-std_ccle_vs_dpsi']] = X %>%
        ggplot(aes(x=ccle_std, y=abs(psi__median_diff))) +
        geom_scattermore(
            pixels = c(1000,1000), 
            pointsize = 5,
            color='#52b788') +
        facet_wrap(~cancer_type, ncol=4) +
        theme_pubr(border=TRUE) +
        labs(x='Event Std. in CCLE', y='|Delta PSI|') +
        geom_hline(yintercept=THRESH_MEDIAN_DIFF, linetype='dashed', size=LINE_SIZE) +
        geom_vline(xintercept=1, linetype='dashed', size=LINE_SIZE) +
        stat_cor(method='spearman', label.sep = '\n', size=2) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    return(plts)
}


plot_tumorigenesis = function(spldep_stats, spldep_stats_ccle, indices){
    X = spldep_stats
    ref = spldep_stats_ccle %>% column_to_rownames("event_gene")
    n_cancers = X %>% distinct(cancer_type) %>% nrow()
    
    plts = list()
    plts[['tumorigenesis-scatters']] = X %>%
        ggplot(aes(x=q25, y=q75)) +
        geom_scattermore(
            pixels = c(1000,1000), 
            pointsize = 8,
            color='black') +
        labs(x='0.25 Quantile', y='0.75 Quantile') +
        geom_hline(yintercept=0, linetype='dashed', size=LINE_SIZE) +
        geom_vline(xintercept=0, linetype='dashed', size=LINE_SIZE) +
        geom_text_repel(aes(label=event_gene), size=FONT_SIZE, segment.size=0.1, family="Arial",
                        X %>% group_by(cancer_type) %>% slice_max(order_by = q25*q75, n=5)) +
        geom_text_repel(aes(label=event_gene), size=FONT_SIZE, segment.size=0.1, family="Arial",
                        X %>% group_by(cancer_type) %>% filter(q75>1)) +
        geom_text_repel(aes(label=event_gene), size=FONT_SIZE, segment.size=0.1, family="Arial",
                        X %>% group_by(cancer_type) %>% filter(q75>0 & q25<(-0.9))) +
        facet_wrap(~cancer_type, ncol=4, scales="free") +
        theme_pubr(border=TRUE) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    idx = X %>% pull(event_gene)
    diffs = data.frame(
        event_gene = idx,
        GENE = ref[idx,'GENE'],
        median_diff = X[['median']] - ref[idx, 'med'],
        range_diff = abs(X[['q25']] - X[['q75']]) - abs(ref[idx,'q25'] - ref[idx,'q75']),
        cancer_type = X[['cancer_type']]
    )
    
    plts[['tumorigenesis-diff_medians']] = diffs %>% 
        ggviolin(x="cancer_type", y="median_diff", fill="cancer_type", color=NA, 
                 palette=get_palette("Paired",n_cancers)) + 
        geom_boxplot(width=0.1, outlier.size=0.1) + 
        guides(fill='none') + 
        geom_hline(yintercept=0, linetype='dashed', size=LINE_SIZE) + 
        geom_text_repel(aes(label=event_gene), 
                        diffs %>% 
                        filter(abs(median_diff)>0.5),
                        size=FONT_SIZE,
                        segment.size=0.1, family="Arial") +  
        theme_pubr(x.text.angle=70) + 
        labs(x="Cancer Type", 
             y="median(Spl. Dep. TCGA) - median(Spl. Dep. CCLE)")
    
    plts[['tumorigenesis-diff_ranges']] = diffs %>% 
        ggviolin(x="cancer_type", y="range_diff", fill="cancer_type", color=NA, 
                 palette=get_palette("Paired",n_cancers)) + 
        geom_boxplot(width=0.1, outlier.size=0.1) + 
        guides(fill='none') + 
        geom_hline(yintercept=0, linetype='dashed', size=LINE_SIZE) + 
        geom_text_repel(aes(label=event_gene), diffs %>% 
                        filter(abs(range_diff)>0.15), family="Arial",
                        size=FONT_SIZE, segment.size=0.1, max.overlaps=15) +  
        theme_pubr(x.text.angle=70) + 
        labs(x="Cancer Type", y="range(Spl. Dep. TCGA) - range(Spl. Dep. CCLE)")
    
    plts[['tumorigenesis-top_diff_ranges']] = diffs %>% 
        group_by(cancer_type) %>% 
        filter(range_diff>0.01) %>% # only higher than 0.05
        count(event_gene) %>%
        ungroup() %>%
        group_by(event_gene) %>%
        mutate(total=sum(n)) %>%
        arrange(total) %>%
        filter(total>5) %>%
        ggbarplot(x="event_gene", y="n", fill="cancer_type", 
                  color=NA, palette=get_palette("Paired",n_cancers)) +
        labs(x="Event & Gene", y="Count", fill="Cancer Type") +
        coord_flip()
    
    X = indices %>% 
        filter(index%in%selected_events) %>% 
        left_join(annot %>% distinct(GENE,EVENT), 
                  by=c("index"="EVENT")) %>% 
        mutate(event_gene=paste0(index,"_",GENE)) %>% 
        group_by(index_name,event_gene) %>% 
        summarize(med=median(correlation))
        
    plts[["tumorigenesis-median_corr_estimate"]] = X %>% 
        ggviolin(x="index_name", y="med", fill="bisque3", color=NA) + 
        geom_boxplot(width=0.1, outlier.size=0.1) + 
        geom_text_repel(aes(label=event_gene), X %>% slice_max(abs(med), n=10),
                        size=FONT_SIZE, family="Arial", segment.size=0.1) +
        labs(x="Index", y="median(Spearman Correlation)")
    
    return(plts)
}


plot_spldep_distrs = function(spldep_stats_ccle_cancers, spldep_stats){
    plts = list()
    
    X = spldep_stats_ccle_cancers %>% 
        pivot_wider(id_cols="EVENT", 
                    names_from="primary_disease", 
                    values_from="median") %>% 
        column_to_rownames("EVENT")
    X = X[,!(colSums(is.na(X))==nrow(X))]
    X = X[!(rowSums(is.na(X))>=(ncol(X)-1)),]
    X = t(apply(X, 1, function(x){replace_na(x, median(x, na.rm=TRUE))}))
    mat = t(scale(t(X)))
    plts[["spldep_distrs-ccle-heat"]] = mat %>% 
        Heatmap(show_row_names = FALSE, 
                column_names_gp=gpar(fontsize=6),
                show_row_dend = FALSE)
    plts[["spldep_distrs-ccle-heat"]] = as.ggplot(grid.grabExpr(draw(plts[["spldep_distrs-ccle-heat"]])))
    
    X = spldep_stats %>% 
        pivot_wider(id_cols="EVENT", 
                    names_from="cancer_type", 
                    values_from="median") %>% 
        column_to_rownames("EVENT")
    X = X[,!(colSums(is.na(X))==nrow(X))]
    X = X[!(rowSums(is.na(X))>=(ncol(X)-1)),]
    X = t(apply(X, 1, function(x){replace_na(x, median(x, na.rm=TRUE))}))
    mat = t(scale(t(X)))
    plts[["spldep_distrs-tcga-heat"]] = mat %>% 
        Heatmap(show_row_names=FALSE, 
                column_names_gp=gpar(fontsize=6),
                show_row_dend = FALSE)
    plts[["spldep_distrs-tcga-heat"]] = as.ggplot(grid.grabExpr(draw(plts[["spldep_distrs-tcga-heat"]])))
    
    return(plts)
}


make_plots = function(diff_result_sample, 
                      diff_result_subtypes, 
                      diff_result_sample_raw, psi_ccle, 
                      spldep_stats, spldep_stats_ccle, indices,
                      cancer_events){
    plts = list(
        plot_top_candidates_sample_type(diff_result_sample, cancer_events),
        plot_top_candidates_sample_type(diff_result_subtypes, cancer_events, 'subtypes-'),
        plot_ccle_vs_tcga(diff_result_sample_raw, psi_ccle),
        plot_tumorigenesis(spldep_stats, spldep_stats_ccle, indices),
        plot_spldep_distrs(spldep_stats_ccle_cancers, spldep_stats)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(diff_result_sample, 
                        diff_result_subtypes, 
                        spldep_stats){
    figdata = list(
        "targetable_events" = list(
            "differential_analysis-by_cancer_type" = diff_result_sample,
            "differential_analysis-by_cancer_subtype" = diff_result_subtypes,
            "splicing_dependency_stats" = spldep_stats
        )
    )
    return(figdata)
}


save_plt = function(plts, plt_name, extension='.pdf', 
                    directory='', dpi=350, format=TRUE,
                    width = par("din")[1], height = par("din")[2]){
    print(plt_name)
    plt = plts[[plt_name]]
    if (format){
        plt = ggpar(plt, font.title=8, font.subtitle=8, font.caption=8, 
                    font.x=8, font.y=8, font.legend=6,
                    font.tickslab=6, font.family=FONT_FAMILY, device=cairo_pdf)   
    }
    filename = file.path(directory,paste0(plt_name,extension))
    save_plot(filename, plt, base_width=width, base_height=height, dpi=dpi, units='cm', )
}


save_plots = function(plts, figs_dir){
    # top candidates sample type
    save_plt(plts, 'top_samples-sample_counts_cancer', '.pdf', figs_dir, width=8, height=5)
    save_plt(plts, 'top_samples-dpsi_vs_spldep-scatter', '.pdf', figs_dir, width=10, height=10)
    save_plt(plts, 'top_samples-dpsi_vs_spldep-selection', '.pdf', figs_dir, width=10, height=6)
    save_plt(plts, 'top_samples-dpsi_vs_spldep-candidates_spldep', '.pdf', figs_dir, width=5.5, height=12.5)
    save_plt(plts, 'top_samples-dpsi_vs_spldep-candidates_dpsi', '.pdf', figs_dir, width=5.5, height=12.5)
    save_plt(plts, 'top_samples-dpsi_vs_spldep-candidates_harm', '.pdf', figs_dir, width=5.5, height=12.5)
    save_plt(plts, 'top_samples-dpsi_vs_spldep-candidates_psi', '.pdf', figs_dir, width=18, height=30)
    save_plt(plts, 'top_samples-enrichment-KD-dot', '.pdf', figs_dir, width=5, height=7)
    
    # top candidates sample type (cancer subtypes)
    save_plt(plts, 'subtypes-top_samples-sample_counts_cancer', '.pdf', figs_dir, width=8, height=8)
    save_plt(plts, 'subtypes-top_samples-dpsi_vs_spldep-scatter', '.pdf', figs_dir, width=10, height=10)
    save_plt(plts, 'subtypes-top_samples-dpsi_vs_spldep-selection', '.pdf', figs_dir, width=10, height=6)
    save_plt(plts, 'subtypes-top_samples-dpsi_vs_spldep-candidates_spldep', '.pdf', figs_dir, width=5.5, height=8)
    save_plt(plts, 'subtypes-top_samples-dpsi_vs_spldep-candidates_dpsi', '.pdf', figs_dir, width=5.5, height=8)
    save_plt(plts, 'subtypes-top_samples-dpsi_vs_spldep-candidates_harm', '.pdf', figs_dir, width=5.5, height=8)
    save_plt(plts, 'subtypes-top_samples-dpsi_vs_spldep-candidates_psi', '.pdf', figs_dir, width=15, height=15)
    save_plt(plts, 'subtypes-top_samples-enrichment-KD-dot', '.pdf', figs_dir, width=5, height=7)

    # CCLE vs TCGA
    save_plt(plts, 'ccle_vs_tcga-std_ccle_vs_dpsi', '.pdf', figs_dir, width=10, height=10)
    
    # tumorigenensis
    save_plt(plts, 'tumorigenesis-scatters', '.pdf', figs_dir, width=15, height=15)
    save_plt(plts, 'tumorigenesis-diff_medians', '.pdf', figs_dir, width=8, height=8)
    save_plt(plts, 'tumorigenesis-diff_ranges', '.pdf', figs_dir, width=8, height=8)
    save_plt(plts, 'tumorigenesis-top_diff_ranges', '.pdf', figs_dir, width=6, height=8)
    save_plt(plts, 'tumorigenesis-median_corr_estimate', '.pdf', figs_dir, width=5, height=5)
    
    # splicing dependency distributions
    save_plt(plts, 'spldep_distrs-ccle-heat', '.pdf', figs_dir, width=8, height=8)
    save_plt(plts, 'spldep_distrs-tcga-heat', '.pdf', figs_dir, width=8, height=8)
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
    diff_result_sample_file = args$diff_result_sample_file
    diff_result_response_file = args$diff_result_response_file
    selected_events_file = args$selected_events_file
    spldep_file = args$spldep_file
    psi_ccle_file = args$psi_ccle_file
    spldep_ccle_file = args$spldep_ccle_file
    figs_dir = args$figs_dir
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    diff_result_sample = read_tsv(diff_result_sample_file)
    diff_result_subtypes = read_tsv(diff_result_subtypes_file)
    annot = read_tsv(annotation_file)
    psi_ccle = read_tsv(psi_ccle_file)
    selected_events = readLines(selected_events_file)
    spldep_stats = read_tsv(spldep_stats_file) %>% 
        filter(EVENT %in% selected_events)
    spldep_stats_subtypes = read_tsv(spldep_stats_subtypes_file) %>% 
        filter(EVENT %in% selected_events)
    spldep_ccle = read_tsv(spldep_ccle_file) %>% 
        filter(index %in% selected_events) 
    spldep_stats_ccle_cancers = read_tsv(spldep_stats_ccle_cancers_file) %>% 
        filter(EVENT %in% selected_events) 
    indices = read_tsv(indices_file) %>% filter(index_name=="ESTIMATE")
    ontology = read_tsv(ontology_file)
    cancer_events = read_tsv(cancer_events_file)
    
    # metadata = read_tsv(metadata_file)
    
    # explore possible drug validation
#     selected_drugs = readLines(selected_drugs_file)
#     drug_targets = read_tsv(drug_targets_file)
#     drug_models = read_tsv(drug_models_file) %>% 
#         filter(ID %in% selected_drugs & lr_padj<0.1) %>%
#         left_join(drug_targets %>% distinct(DRUG_ID,DRUG_NAME), by="DRUG_ID") %>%
#         group_by(ID) %>%
#         arrange(-spldep_coefficient) %>%
#         mutate(event_gene = paste0(EVENT,"_",GENE),
#                event_type = gsub("Hsa","",gsub("[^a-zA-Z]", "",EVENT)),
#                ranking = row_number()) %>%
#         ungroup()
#     drug_treatments = read_tsv(drug_treatments_file) %>% 
#         mutate(DRUG_NAME=toupper(DRUG_NAME))
    
#     x = 
#     x %>% 
#         ungroup() %>%
#         left_join(drug_models %>% 
#         distinct(event_gene, DRUG_ID,DRUG_NAME,ranking,spldep_coefficient), by="event_gene") %>%
#         drop_na(ranking) %>%
#         filter(DRUG_NAME %in% drugs_oi) %>%
#         arrange(mean) %>%
#         distinct(event_gene, DRUG_NAME,DRUG_ID,ranking, spldep_coefficient, cancer_type, mean) %>%
#         filter(ranking <= 10)
    
#     drugs_oi = drug_treatments %>% count(DRUG_NAME) %>% filter(n>50) %>% pull(DRUG_NAME)
#     drug_treatments %>% filter(DRUG_NAME=="TEMSIROLIMUS")

    # add event gene
    spldep_stats = spldep_stats %>%
        left_join(annot[,c('EVENT','GENE')], by='EVENT') %>%
        mutate(event_gene = paste0(EVENT,'_',GENE))
    spldep_stats_subtypes = spldep_stats_subtypes %>%
        left_join(annot[,c('EVENT','GENE')], by='EVENT') %>%
        mutate(event_gene = paste0(EVENT,'_',GENE))
    
    # prep results differential analyses
    ## cancer types
    diff_result_sample_raw = diff_result_sample
    diff_result_sample = prep_diff_result(diff_result_sample, spldep_stats)
    ## cancer subtypes
    diff_result_subtypes_raw = diff_result_subtypes
    diff_result_subtypes = prep_diff_result(
        diff_result_subtypes %>%
            mutate(cancer = cancer_type, 
                   cancer_type = paste0(cancer_type,'_',cancer_subtype)), 
        spldep_stats_subtypes %>%
            mutate(cancer = cancer_type, 
                   cancer_type = paste0(cancer_type,'_',cancer_subtype)) %>%
            dplyr::select(-c(cancer, cancer_subtype)))
    
    # prep summary stats CCLE
    X = spldep_ccle %>%
        column_to_rownames('index')
    spldep_stats_ccle = data.frame(
            med = apply(X,1,median, na.rm=TRUE),
            std = apply(X,1,sd, na.rm=TRUE),
            min = apply(X,1,min, na.rm=TRUE),
            max = apply(X,1,max, na.rm=TRUE),
            q75 = apply(X,1,quantile, na.rm=TRUE, probs=0.75),
            q25 = apply(X,1,quantile, na.rm=TRUE, probs=0.25),
            n_missing = rowSums(is.na(X))
        ) %>% rownames_to_column('EVENT') %>%
        left_join(annot[,c('EVENT','GENE')], by='EVENT') %>%
        mutate(event_gene = paste0(EVENT,'_',GENE))
    
    # plot
    plts = make_plots(diff_result_sample, 
                      diff_result_subtypes, 
                      diff_result_sample_raw, psi_ccle, 
                      spldep_stats, spldep_stats_ccle, indices,
                      cancer_events)
    
    # make figdata
    figdata = make_figdata(diff_result_sample, 
                           diff_result_subtypes, 
                           spldep_stats)

    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}