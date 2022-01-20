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

require(tidyverse)
require(tidytext)
require(ggpubr)
require(cowplot)
require(scattermore)
require(ggrepel)
require(clusterProfiler)
require(ComplexHeatmap)
require(impute)
require(gridExtra)
require(ggplotify)
require(extrafont)

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

# variables
THRESH_FDR = 0.05
THRESH_MEDIAN_DIFF = 5
MIN_SAMPLES = 10

# Development
# -----------
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# RESULTS_DIR = file.path(ROOT,'results','splicing_dependency_tcga')
# diff_result_sample_file = file.path(RESULTS_DIR,'files','PANCAN','mannwhitneyu-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz')
# diff_result_response_file = file.path(RESULTS_DIR,'files','PANCAN','mannwhitneyu-RESPONDER_vs_NONRESPONDER-EX.tsv.gz')
# annotation_file = file.path(RAW_DIR,'VastDB','EVENT_INFO-hg38_noseqs.tsv')
# gsea_result_file = file.path(RESULTS_DIR,'files','PANCAN','gsea-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz')
# figs_dir = file.path(RESULTS_DIR,'figures','targetable_events-EX')
# selected_events_file = file.path(ROOT,'results','model_splicing_dependency','files','selected_models-EX.txt')
# spldep_file = file.path(RESULTS_DIR,'files','BRCA','splicing_dependency-EX','mean.tsv.gz')
# spldep_luad_file = file.path(RESULTS_DIR,'files','LUAD','splicing_dependency_mean-EX.tsv.gz')
# spldep_lgg_file = file.path(RESULTS_DIR,'files','LGG','splicing_dependency_mean-EX.tsv.gz')
# metadata_response_file = file.path(PREP_DIR,'Moiso2021','drug_response.tsv.gz')
# metadata_file = file.path(PREP_DIR,'metadata','PANCAN.tsv.gz')
# psi_ccle_file = file.path(PREP_DIR,'event_psi','CCLE-EX.tsv.gz')
# spldep_stats_file = file.path(RESULTS_DIR,'files','PANCAN','summary_splicing_dependency-EX.tsv.gz')
# spldep_ccle_file = file.path(ROOT,'results','model_splicing_dependency','files','splicing_dependency-EX','mean.tsv.gz')
# MODELS_DIR = file.path(ROOT,'results','splicing_dependency_drugs')
# models_file = file.path(MODELS_DIR,'files','model_summaries_drug_response-EX.tsv.gz')
# drug_targets_file = file.path(RAW_DIR,'GDSC','screened_compunds_rel_8.2.csv')
# msigdb_dir = file.path(ROOT,'data','raw','MSigDB','msigdb_v7.4','msigdb_v7.4_files_to_download_locally','msigdb_v7.4_GMTs')
# diff_result_subtypes_file = file.path(RESULTS_DIR,'files','PANCAN_subtypes','mannwhitneyu-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz')
# spldep_stats_subtypes_file = file.path(RESULTS_DIR,'files','PANCAN_subtypes','summary_splicing_dependency-EX.tsv.gz')
# embedding_file = file.path(RESULTS_DIR,'files','PANCAN','embedded_splicing_dependency_mean-EX.tsv.gz')
# spldep_stats_ccle_cancers_file = file.path(ROOT,'results','model_splicing_dependency','files','splicing_dependency_summaries-EX','primary_disease.tsv.gz')

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


plot_top_candidates_sample_type = function(diff_result, models, drug_targets, rankings, patt=''){
    
    plts = list()

    # how many samples do we have for each cancer and sample type?
    X = diff_result %>% 
        mutate(PT = `psi__condition_a-n_present` + `psi__condition_a-n_nan`,
               STN = `psi__condition_b-n_present` + `psi__condition_b-n_nan`) %>% 
        distinct(cancer_type,PT,STN) %>%
        pivot_longer(c(PT,STN),names_to='sample_type',values_to='n')
    
    plts[['top_samples-sample_counts_cancer']] = X %>% 
        ggbarplot(x='cancer_type', y='n', fill='sample_type', color=NA,
                  position=position_dodge(0.7), label=TRUE, lab.size=1, palette='lancet') + 
        labs(x='Cancer Type', y='N. Samples') + 
        theme_pubr(x.text.angle = 45) + 
        geom_hline(yintercept=10, linetype='dashed')
    
    # ranking with differential analyses
    cancers_oi = X %>% 
        group_by(cancer_type) %>% 
        summarize(to_keep=sum(n >= MIN_SAMPLES)==2) %>% 
        filter(to_keep) %>% 
        pull(cancer_type)
    
    X = diff_result %>%
        filter(cancer_type %in% cancers_oi) %>%
        mutate(sign_dpsi = sign(psi__median_diff),
               sign_spldep = sign(mean))
    
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
        geom_hline(yintercept=0, linetype='dashed') +
        geom_vline(xintercept=0, linetype='dashed') +
        theme(strip.text.x = element_text(size=6, family='Arial'))
    
    plts[['top_samples-dpsi_vs_spldep-selection']] = X %>%
        filter(psi__is_significant) %>%
        mutate(sign_combined = sprintf('sign_dpsi=%s & sign_spldep=%s',sign_dpsi,sign_spldep)) %>%
        count(sign_combined) %>% 
        drop_na() %>%
        ggbarplot(x='cancer_type', y='n', fill='sign_combined', palette='jco', 
                  color=FALSE, label=TRUE, lab.size=1, position=position_dodge(0.9)) +
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
        geom_hline(yintercept=0, linetype='dashed') +
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
        geom_hline(yintercept=c(-5,5), linetype='dashed') +
        labs(x='Event & Gene', y='Delta PSI', fill='Cancer Type') +
        coord_flip()
    
    x = X %>% 
        filter(cancer_type %in% cancers_oi) %>%
        filter(psi__is_significant & 
               ((sign_dpsi>0 & sign_spldep<0) | (sign_dpsi<0 & sign_spldep>0)))
    a = x %>% dplyr::select(c("event_gene","psi__condition_a-median",
                              "psi__condition_a-mad","psi__condition_a"))
    colnames(a) = gsub("_a","",colnames(a))
    b = x %>% dplyr::select(c("event_gene","psi__condition_b-median",
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
        theme(strip.text.x = element_text(size=6, family='Arial'))
        
    
    # which drugs are associated with targetable exons?
    targetable_events = X %>% 
        filter(cancer_type %in% cancers_oi) %>%
        filter(psi__is_significant & 
               ((sign_dpsi>0 & sign_spldep<0) | (sign_dpsi<0 & sign_spldep>0))) %>%
        pull(EVENT) %>%
        unique()
    
    plts[['top_samples-candidates-drug_assocs']] = models %>% 
        left_join(drug_targets %>% distinct(DRUG_ID,DRUG_NAME), by="DRUG_ID") %>%
        left_join(drug_targets %>% distinct(DRUG_ID,TARGET) %>% mutate(is_target=TRUE), 
                  by=c("DRUG_ID", "GENE"="TARGET")) %>%
        mutate(is_target=!is.na(is_target)) %>%
        filter(lr_padj<0.1 & EVENT%in%targetable_events) %>%
        group_by(event_gene) %>% 
        slice_max(order_by=spldep_coefficient, n=15) %>% 
        mutate(DRUG_NAME=as.factor(DRUG_NAME),
               name=reorder_within(DRUG_NAME, spldep_coefficient, event_gene)) %>% 
        ggplot(aes(x=name, y=spldep_coefficient, fill)) +
        geom_col(aes(fill=is_target, color=drug_screen, group=drug_screen), 
                     position=position_dodge(1)) +
        color_palette(c("white","black")) + 
        fill_palette("Set2") +    
        facet_wrap(~event_gene, scales='free', ncol=4) +
        coord_flip() +
        scale_x_reordered() +     
        labs(x='Drug', y='Effect Size', color='Drug Screen', fill='Is Drug Target') + 
        theme_pubr(border=TRUE) +
        theme(strip.text.x = element_text(size=6, family='Arial'))
    
    X = rankings %>% 
        filter(EVENT %in% targetable_events) %>% 
        group_by(event_gene) %>% 
        slice_min(combined_ranking, n=1) %>% 
        ungroup() %>% 
        arrange(index) %>%
        left_join(drug_targets %>% 
                      distinct(DRUG_ID,TARGET) %>% 
                      mutate(is_target=TRUE), 
                  by=c("DRUG_ID", "GENE"="TARGET")) %>%
        mutate(is_target=!is.na(is_target))
    
    plts[['top_samples-candidates-drug_assocs_best']] = X %>%
        ggbarplot(x="event_gene", y="combined_ranking", 
                  fill="is_target", color="drug_screen") +
        color_palette(c("white","black")) + 
        fill_palette("Set2") + 
        geom_text(aes(label=index), X, size=1, family='Arial') +
        geom_text(aes(y=0.015, label=DRUG_NAME), X, 
                  size=1, family='Arial') +
        labs(x="Event & Gene", y="Ranking Ratio Sum", 
             fill="Is Drug Target", color="Drug Screen") +
        coord_flip()
        
    
    
    names(plts) = paste0(patt,names(plts))
    
    return(plts)
}


plot_spldeps_pt = function(diff_result, spldep, metadata){
    plts=list()
    
    # overall median splicing dependencies
    X = diff_result
    mat = X %>% 
        pivot_wider(id_cols = 'event_gene', names_from = 'cancer_type', 
                    values_from = 'spldep__condition_a-median') %>% 
        column_to_rownames('event_gene') %>% 
        as.matrix()
    imputed = impute.knn(mat)
    mat = imputed[['data']]
    pca = apply(mat,1,function(x){ (x-mean(x))/sd(x)}) %>%
        t() %>%
        prcomp()
    pcs = pca[['rotation']] %>% 
        as.data.frame() %>% 
        rownames_to_column('cancer_type')
    n = length(unique(pcs[['cancer_type']]))
    plts[['spldeps_pt-pancan-pca']] = pcs %>% 
        ggscatter(x='PC1', y='PC2', color='cancer_type') +
        geom_text_repel(aes(label=cancer_type, color=cancer_type), family="Arial") +
        ggpubr::color_palette(get_palette('Paired',n)) +
        guides(color='none')
    
    
    # BRCA subtypes
    samples_oi = metadata %>% 
        filter(sample_type=='Primary Tumor' & 
               PAM50Call_RNAseq!='Unknown') %>% 
        pull(sampleID) %>% 
        intersect(colnames(spldep))
    annot = metadata %>% 
        filter(sampleID %in% samples_oi) %>%
        arrange(match(sampleID, samples_oi))
    X = spldep[,c('event_gene',samples_oi)]
    mat = X %>% column_to_rownames('event_gene') %>% as.matrix()
    imputed = impute.knn(mat, colmax = 0.9)
    mat = imputed[['data']]
#     vals = unique(annot[['PAM50Call_RNAseq']]) %>% na.omit()
#     ha = HeatmapAnnotation(subtype=annot[['PAM50Call_RNAseq']], 
#                            col = list(subtype=setNames(get_palette('Dark2', k=length(vals)),vals)))
#     ht = mat %>% Heatmap(show_row_names = FALSE, bottom_annotation = ha)
#     plts[['spldeps_pt-brca_subtypes-heatmap']] = draw(ht)
    
    pca = apply(mat,1,function(x){ (x-mean(x))/sd(x)}) %>%
        t() %>%
        prcomp()
    pcs = pca[['rotation']] %>% 
        as.data.frame() %>% 
        rownames_to_column('sampleID') %>%
        left_join(annot, by='sampleID')
    plts[['spldeps_pt-brca_subtypes-counts']] = pcs %>% 
        count(PAM50Call_RNAseq) %>% 
        ggbarplot(x='PAM50Call_RNAseq', y='n', fill='PAM50Call_RNAseq', 
                  color=NA, label=TRUE, lab.size=1) + guides(fill='none') + 
        labs(y='N. Samples')
    plts[['spldeps_pt-brca_subtypes-pca']] = pcs %>% 
        ggscatter(x='PC1', y='PC2', color='PAM50Call_RNAseq', 
                  palette='Set1', alpha=0.5)
    X = spldep %>% 
        pivot_longer(-event_gene, names_to = 'sampleID', values_to='spldep') %>% 
        left_join(metadata[,c('sampleID','PAM50Call_RNAseq')], by='sampleID') %>% 
        drop_na() %>%
        filter(PAM50Call_RNAseq!='Unknown') %>%
        group_by(PAM50Call_RNAseq,event_gene) %>% 
        summarize(med_spldep=median(spldep)) %>% 
        arrange(med_spldep)
    
    plts[['spldeps_pt-brca_subtypes-bar']] = X %>%
        ggbarplot(x='event_gene', y='med_spldep', fill='PAM50Call_RNAseq', 
                  color=NA, position=position_dodge(0.9)) + 
        labs(x='Event & Gene', y='median(Spl. Dep.)') + 
        coord_flip()
    
    X = X %>% 
        ungroup() %>% 
        group_by(event_gene) %>% 
        mutate(prop_spldep=abs(med_spldep)/sum(abs(med_spldep)), 
               prop_spldep_sign=sign(med_spldep)*prop_spldep) %>% 
        arrange(prop_spldep_sign) 
    plts[['spldeps_pt-brca_subtypes-specificity']] = X %>%
        ggscatter(x='prop_spldep', y='med_spldep', color='PAM50Call_RNAseq', alpha=0.7) +
        geom_text_repel(
            aes(label=event_gene),
            X %>%
            filter(prop_spldep>0.25 & abs(med_spldep)>0.25) %>%
            ungroup() %>%
            arrange(event_gene),
            max.overlaps=50,
            family="Arial"
        ) +
        labs(x='Proportion of Spl. Dep.', y='median(Spl. Dep.)')
    
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
        geom_hline(yintercept=THRESH_MEDIAN_DIFF, linetype='dashed') +
        geom_vline(xintercept=1, linetype='dashed') +
        stat_cor(method='spearman', label.sep = '\n', size=2) +
        theme(strip.text.x = element_text(size=6, family='Arial'))
    
    return(plts)
}


plot_tumorigenesis = function(spldep_stats, spldep_stats_ccle){
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
        labs(x='0.05 Quantile', y='0.95 Quantile') +
        geom_hline(yintercept=0, linetype='dashed') +
        geom_vline(xintercept=0, linetype='dashed') +
        geom_text_repel(aes(label=event_gene), size=1, segment.size=0.1, family="Arial",
                        X %>% group_by(cancer_type) %>% slice_max(order_by = q25*q75, n=5)) +
        geom_text_repel(aes(label=event_gene), size=1, segment.size=0.1, family="Arial",
                        X %>% group_by(cancer_type) %>% filter(q75>1)) +
        geom_text_repel(aes(label=event_gene), size=1, segment.size=0.1, family="Arial",
                        X %>% group_by(cancer_type) %>% filter(q75>0 & q25<(-0.9))) +
        facet_wrap(~cancer_type, ncol=4, scales="free") +
        theme_pubr(border=TRUE) +
        theme(strip.text.x = element_text(size=6, family='Arial'))
    
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
        geom_hline(yintercept=0, linetype='dashed') + 
        geom_text_repel(aes(label=event_gene), 
                        diffs %>% 
                        filter(abs(median_diff)>0.5),
                        size=1,
                        segment.size=0.1, family="Arial") +  
        theme_pubr(x.text.angle=70) + 
        labs(x="Cancer Type", 
             y="median(Spl. Dep. TCGA) - median(Spl. Dep. CCLE)")
    
    plts[['tumorigenesis-diff_ranges']] = diffs %>% 
        ggviolin(x="cancer_type", y="range_diff", fill="cancer_type", color=NA, 
                 palette=get_palette("Paired",n_cancers)) + 
        geom_boxplot(width=0.1, outlier.size=0.1) + 
        guides(fill='none') + geom_hline(yintercept=0, linetype='dashed') + 
        geom_text_repel(aes(label=event_gene), diffs %>% 
                        filter(abs(range_diff)>0.15), family="Arial",
                        size=1, segment.size=0.1, max.overlaps=15) +  
        theme_pubr(x.text.angle=70) + 
        labs(x="Cancer Type", y="range(Spl. Dep. TCGA) - range(Spl. Dep. CCLE)")
    
    plts[['tumorigenesis-top_diff_ranges']] = diffs %>% 
        group_by(cancer_type) %>% 
        filter(range_diff>0.05) %>% # only higher than 0.05
        count(event_gene) %>%
        ungroup() %>%
        group_by(event_gene) %>%
        mutate(total=sum(n)) %>%
        arrange(total) %>%
        filter(total>1) %>%
        ggbarplot(x="event_gene", y="n", fill="cancer_type", 
                  color=NA, palette=get_palette("Paired",n_cancers)) +
        labs(x="Event & Gene", y="Count", fill="Cancer Type") +
        coord_flip()
    
    return(plts)
}


plot_patient_clusters = function(embedding, metadata){
    X = embedding %>% 
        left_join(metadata %>% distinct(sampleID, cancer, sample_type, OS, OS.time), 
                  by=c("index"="sampleID")) %>%
        filter(sample_type=="Primary Tumor")
    n_cancers = length(unique(X[["cancer"]]))
    n_clusters = length(unique(X[["leiden_labels"]]))
    
    plts = list()
    plts[["pat_clust-umap-cancers"]] = X %>% 
        ggplot(aes(x=UMAP0, y=UMAP1)) + 
        geom_scattermore(aes(color=cancer), pointsize=5, 
                         pixels=c(1000,1000), alpha=0.5) + 
        color_palette(get_palette("Paired", n_cancers)) + 
        theme_pubr() + 
        guides(color = guide_legend(override.aes = list(alpha = 1))) + 
        labs(color="Cancer Type") + 
        theme(aspect.ratio=1)
    
    plts[["pat_clust-umap-clusters"]] = X %>% 
        ggplot(aes(x=UMAP0, y=UMAP1)) + 
        geom_scattermore(aes(color=leiden_labels), pointsize=5, 
                         pixels=c(1000,1000), alpha=0.5) + 
        color_palette(get_palette("default", n_clusters)) + 
        theme_pubr() + 
        guides(color = guide_legend(override.aes = list(alpha = 1))) + 
        labs(color="Cluster") + 
        theme(aspect.ratio=1)
    
    plts[["pat_clust-clusters_vs_cancers-balloon"]] = X %>% 
        count(leiden_labels, cancer) %>%
        group_by(leiden_labels) %>%
        mutate(perc = n / sum(n)) %>%
        ggballoonplot(x="cancer", y="leiden_labels", size="perc", 
                      fill="tomato3", color="white") +
        scale_size(range=c(1,3)) + 
        labs(x="Cancer Type", y="Cluster")
    
    plts[["pat_clust-clusters_vs_cancers-bar"]] = X %>% 
        count(leiden_labels, cancer) %>%
        ggbarplot(x="leiden_labels", y="n", fill="cancer", position=position_dodge(0.9),
                  color=NA, palette=get_palette("Paired", n_cancers)) +
        labs(x="Cluster", y="Count", fill="Cancer Type") +
        coord_flip()
    
        
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


make_plots = function(diff_result_sample, diff_result_subtypes, spldep, 
                      metadata, diff_result_sample_raw, psi_ccle, spldep_stats, 
                      spldep_stats_ccle, selected_events, models, drug_targets, 
                      rankings, embedding){
    plts = list(
        plot_top_candidates_sample_type(diff_result_sample, models, 
                                        drug_targets, rankings),
        plot_top_candidates_sample_type(diff_result_subtypes, models, 
                                        drug_targets, rankings, 'subtypes-'),
        plot_spldeps_pt(diff_result_sample, spldep, metadata),
        plot_ccle_vs_tcga(diff_result_sample_raw, psi_ccle),
        plot_tumorigenesis(spldep_stats, spldep_stats_ccle),
        plot_patient_clusters(embedding, metadata),
        plot_spldep_distrs(spldep_stats_ccle_cancers, spldep_stats)
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
    # top candidates sample type
    save_plt(plts, 'top_samples-sample_counts_cancer', '.pdf', figs_dir, width=8, height=5)
    save_plt(plts, 'top_samples-dpsi_vs_spldep-scatter', '.pdf', figs_dir, width=10, height=10)
    save_plt(plts, 'top_samples-dpsi_vs_spldep-selection', '.pdf', figs_dir, width=10, height=6)
    save_plt(plts, 'top_samples-dpsi_vs_spldep-candidates_spldep', '.pdf', figs_dir, width=8, height=11)
    save_plt(plts, 'top_samples-dpsi_vs_spldep-candidates_dpsi', '.pdf', figs_dir, width=6, height=11)
    save_plt(plts, 'top_samples-dpsi_vs_spldep-candidates_psi', '.pdf', figs_dir, width=18, height=30)
    save_plt(plts, 'top_samples-candidates-drug_assocs', '.pdf', figs_dir, width=25, height=40)
    save_plt(plts, 'top_samples-candidates-drug_assocs_best', '.pdf', figs_dir, width=7, height=10)
    
    # top candidates sample type (cancer subtypes)
    save_plt(plts, 'subtypes-top_samples-sample_counts_cancer', '.pdf', figs_dir, width=8, height=8)
    save_plt(plts, 'subtypes-top_samples-dpsi_vs_spldep-scatter', '.pdf', figs_dir, width=10, height=10)
    save_plt(plts, 'subtypes-top_samples-dpsi_vs_spldep-selection', '.pdf', figs_dir, width=10, height=6)
    save_plt(plts, 'subtypes-top_samples-dpsi_vs_spldep-candidates_spldep', '.pdf', figs_dir, width=8, height=9)
    save_plt(plts, 'subtypes-top_samples-dpsi_vs_spldep-candidates_dpsi', '.pdf', figs_dir, width=8, height=9)
    save_plt(plts, 'subtypes-top_samples-dpsi_vs_spldep-candidates_psi', '.pdf', figs_dir, width=15, height=15)
    save_plt(plts, 'subtypes-top_samples-candidates-drug_assocs', '.pdf', figs_dir, width=27, height=27)
    save_plt(plts, 'subtypes-top_samples-candidates-drug_assocs_best', '.pdf', figs_dir, width=8, height=8)
    
    # splicing dependencies in primary tumors
    ## PANCAN
    save_plt(plts, 'spldeps_pt-pancan-pca', '.pdf', figs_dir, width=6, height=6)
    ## BRCA
    save_plt(plts, 'spldeps_pt-brca_subtypes-counts', '.pdf', figs_dir, width=6, height=6)
    save_plt(plts, 'spldeps_pt-brca_subtypes-bar', '.pdf', figs_dir, width=10, height=15)
    save_plt(plts, 'spldeps_pt-brca_subtypes-pca', '.pdf', figs_dir, width=6, height=6)
    save_plt(plts, 'spldeps_pt-brca_subtypes-specificity', '.pdf', figs_dir, width=6, height=6)
    
    # CCLE vs TCGA
    save_plt(plts, 'ccle_vs_tcga-std_ccle_vs_dpsi', '.pdf', figs_dir, width=10, height=10)
    
    # tumorigenensis
    save_plt(plts, 'tumorigenesis-scatters', '.pdf', figs_dir, width=15, height=15)
    save_plt(plts, 'tumorigenesis-diff_medians', '.pdf', figs_dir, width=8, height=8)
    save_plt(plts, 'tumorigenesis-diff_ranges', '.pdf', figs_dir, width=8, height=8)
    save_plt(plts, 'tumorigenesis-top_diff_ranges', '.pdf', figs_dir, width=6, height=8)
    
    # patient clustering
    save_plt(plts, 'pat_clust-umap-cancers', '.pdf', figs_dir, width=9, height=9)
    save_plt(plts, 'pat_clust-umap-clusters', '.pdf', figs_dir, width=9, height=9)
    save_plt(plts, 'pat_clust-clusters_vs_cancers-balloon', '.pdf', figs_dir, width=9, height=8)
    save_plt(plts, 'pat_clust-clusters_vs_cancers-bar', '.pdf', figs_dir, width=5, height=12)
    
    # splicing dependency distributions
    save_plt(plts, 'spldep_distrs-ccle-heat', '.pdf', figs_dir, width=8, height=8)
    save_plt(plts, 'spldep_distrs-tcga-heat', '.pdf', figs_dir, width=8, height=8)
}


main = function(){
    args = getParsedArgs()
    diff_result_sample_file = args$diff_result_sample_file
    diff_result_response_file = args$diff_result_response_file
    selected_events_file = args$selected_events_file
    metadata_file = args$metadata_file
    spldep_file = args$spldep_file
    # gsea_result_file = args$gsea_result_file
    annotation_file = args$annotation_file
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
    embedding = read_tsv(embedding_file) %>%
        mutate(leiden_labels=as.factor(leiden_labels))
    spldep_stats_ccle_cancers = read_tsv(spldep_stats_ccle_cancers_file) %>% 
        filter(EVENT %in% selected_events) 
    
    # add event gene
    spldep_stats = spldep_stats %>%
        left_join(annot[,c('EVENT','GENE')], by='EVENT') %>%
        mutate(event_gene = paste0(EVENT,'_',GENE))
    spldep_stats_subtypes = spldep_stats_subtypes %>%
        left_join(annot[,c('EVENT','GENE')], by='EVENT') %>%
        mutate(event_gene = paste0(EVENT,'_',GENE))
    
    # prep results differential analyses
    diff_result_sample_raw = diff_result_sample
    diff_result_sample = prep_diff_result(diff_result_sample, spldep_stats)
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
    
    # BRCA
    metadata = read_tsv(metadata_file) %>%
        mutate(PAM50Call_RNAseq = ifelse(is.na(PAM50Call_RNAseq),'Unknown',PAM50Call_RNAseq))
    spldep = read_tsv(spldep_file) %>%
        filter(index %in% selected_events) %>%
        left_join(
            diff_result_sample[,c('EVENT','event_gene')] %>% distinct(),
            by = c('index'='EVENT')
        ) %>%
        dplyr::select(-index)
        
    # drug screen modeling
    models = read_tsv(models_file) %>% 
        mutate(event_gene = paste0(EVENT,'_',GENE),
               event_type = gsub('Hsa','',gsub("[^a-zA-Z]", "",EVENT)))
    drug_targets = read_csv(drug_targets_file) %>% 
        mutate(DRUG_NAME=toupper(DRUG_NAME)) %>%
        dplyr::select(DRUG_ID,DRUG_NAME,TARGET,TARGET_PATHWAY) %>%
        separate_rows(TARGET) %>%
        distinct()
    
    # drug-event rankings
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
    
    # plot
    plts = make_plots(diff_result_sample, diff_result_subtypes, spldep, 
                      metadata, diff_result_sample_raw, psi_ccle, spldep_stats, 
                      spldep_stats_ccle, selected_events, models, drug_targets, 
                      rankings, embedding)

    # save
    save_plots(plts, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}