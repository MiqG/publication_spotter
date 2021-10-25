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
require(ggpubr)
require(cowplot)
require(scattermore)
require(ggrepel)
require(clusterProfiler)
require(ComplexHeatmap)
require(impute)
require(gridExtra)
require(ggplotify)

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
# figs_dir = file.path(RESULTS_DIR,'figures','differential_analysis-EX')
# selected_events_file = file.path(ROOT,'results','model_splicing_dependency','files','selected_models-EX.txt')
# spldep_file = file.path(RESULTS_DIR,'files','BRCA','splicing_dependency_mean-EX.tsv.gz')
# spldep_luad_file = file.path(RESULTS_DIR,'files','LUAD','splicing_dependency_mean-EX.tsv.gz')
# spldep_lgg_file = file.path(RESULTS_DIR,'files','LGG','splicing_dependency_mean-EX.tsv.gz')
# metadata_response_file = file.path(PREP_DIR,'Moiso2021','drug_response.tsv.gz')
# metadata_file = file.path(PREP_DIR,'metadata','PANCAN.tsv.gz')
# psi_ccle_file = file.path(PREP_DIR,'event_psi','CCLE-EX.tsv.gz')

##### FUNCTIONS #####
prep_diff_result = function(diff_result, annot, selected_events){
    diff_result = diff_result %>%
        rename_all(recode, index = "EVENT") %>%
        mutate(event_type = gsub('Hsa','',gsub("[^a-zA-Z]", "",EVENT))) %>%
        left_join(annot[,c('EVENT','GENE')], by='EVENT') %>%
        mutate(event_gene = paste0(EVENT,'_',GENE),
               is_selected = EVENT %in% selected_events) %>%
        filter(is_selected) %>%
        group_by(cancer_type) %>%
        mutate(psi__padj = p.adjust(psi__pvalue, 'fdr'),
               spldep__padj = p.adjust(spldep__pvalue, 'fdr'),
               psi__log10_padj = -log10(psi__padj),
               spldep__log10_padj = -log10(spldep__padj)) %>%
        mutate(psi__is_significant = psi__padj<THRESH_FDR & 
                                     abs(psi__median_diff)>THRESH_MEDIAN_DIFF,
               spldep__is_significant = spldep__padj<THRESH_FDR &
                                        psi__is_significant)
    return(diff_result)
}


plot_top_candidates_response = function(diff_result){
    
    plts = list()
    
    # how many samples do we have for each cancer - drug and treatment cond.?
    X = diff_result %>% 
        mutate(responder = `psi__condition_a-n_present` + `psi__condition_a-n_nan`,
               non_responder = `psi__condition_b-n_present` + `psi__condition_b-n_nan`) %>% 
        distinct(cancer_type,treatment,responder,non_responder) %>%
        pivot_longer(c(responder,non_responder),names_to='response',values_to='n')
    plts[['top_response-sample_counts_cancer']] = X %>% 
        ggbarplot(x='cancer_type', y='n', fill='response', color=NA, 
                  position=position_dodge(0.7), label=TRUE, palette='lancet') + 
        labs(x='Cancer Type', y='No. Samples') + 
        theme_pubr(x.text.angle = 45) + 
        geom_hline(yintercept=10, linetype='dashed')
    
    plts[['top_response-sample_counts_treatment']] = X %>% 
        ggbarplot(x='treatment', y='n', fill='cancer_type', color=NA, 
                  facet.by='response', position=position_dodge(0.7), label=TRUE,
                  palette=get_palette('Paired', length(unique(X[['cancer_type']])))) + 
        labs(x='Treatment', y='No. Samples') + 
        geom_hline(yintercept=10, linetype='dashed') +
        coord_flip()
    
    # differential analyses
    cancers_oi = c('LGG','LUAD','COAD') # with at least 10 samples per condition
    X = diff_result
    
    plts[['top_response-diff_psi-volcanos']] = X %>%
        filter(cancer_type %in% cancers_oi) %>%
        ggplot(aes(x=psi__median_diff, 
                   y=psi__log10_padj, 
                   color=psi__is_significant)) +
        geom_point(size=1) +
        facet_wrap(~cancer_type, ncol=4) +
        theme_pubr(border=TRUE) +
        guides(color="none") +
        labs(x='Delta PSI', y='-log10(FDR)') +
        ggpubr::color_palette(c("black","orange")) +
        geom_text(
            aes(x=x,y=y,label=n),
            X %>% 
            filter(psi__is_significant & cancer_type%in%cancers_oi) %>% 
            mutate(is_pos = psi__median_diff>0) %>%
            count(cancer_type, psi__is_significant, is_pos) %>%
            mutate(x=ifelse(is_pos, 15, -15),
                   y=2.5),
            color='black'
        )
    
    plts[['top_response-diff_spldep-volcanos']] = X %>%
        filter(cancer_type %in% cancers_oi) %>%
        ggplot(aes(x=spldep__median_diff, 
                   y=spldep__log10_padj, 
                   color=spldep__is_significant)) +
        geom_point(size=1) +
        facet_wrap(~cancer_type, ncol=4) +
        theme_pubr(border=TRUE) +
        guides(color="none") +
        labs(x='Delta Spl. Dep.', y='-log10(FDR)') +
        ggpubr::color_palette(c("black","orange")) +
        geom_text(
            aes(x=x,y=y,label=n),
            X %>% 
            filter(spldep__is_significant & cancer_type%in%cancers_oi) %>% 
            mutate(is_pos = spldep__median_diff>0) %>%
            count(cancer_type, spldep__is_significant, is_pos) %>%
            mutate(x=ifelse(is_pos, 0.45, -0.45),
                   y=2.5),
            color='black'
        )
    
    plts[['top_response-diff_spldep_vs_psi-scatters']] = X %>%
        filter(cancer_type %in% cancers_oi) %>%
        ggplot(aes(x=spldep__median_diff,
                   y=psi__median_diff)) +
        geom_point(size=1) +
        facet_wrap(~cancer_type, ncol=4, scales='free') +
        theme_pubr(border=TRUE) +
        stat_cor(method='spearman') +
        labs(x='Delta Splicing Dep.', y='Delta PSI') +
        geom_hline(yintercept=0, linetype='dashed') +
        geom_vline(xintercept=0, linetype='dashed')
    
    
    plts[['top_response-candidates']] = X %>%
        filter(cancer_type%in%cancers_oi & spldep__is_significant) %>%
        ggbarplot(x='event_gene', y='spldep__median_diff', fill='cancer_type', 
                  color=NA, palette='Paired', position=position_dodge(0.7)) +
        labs(x='Event & Gene', y='Delta Spl. Dep.', fill='Cancer Type') + 
        coord_flip()
    
    return(plts)
}


plot_top_candidates_sample_type = function(diff_result){
    plts = list()

    # how many samples do we have for each cancer and sample type?
    X = diff_result %>% 
        mutate(PT = `psi__condition_a-n_present` + `psi__condition_a-n_nan`,
               STN = `psi__condition_b-n_present` + `psi__condition_b-n_nan`) %>% 
        distinct(cancer_type,PT,STN) %>%
        pivot_longer(c(PT,STN),names_to='sample_type',values_to='n')
    plts[['top_samples-sample_counts_cancer']] = X %>% 
        ggbarplot(x='cancer_type', y='n', fill='sample_type', color=NA, lab.size=3,
                  position=position_dodge(0.7), label=TRUE, palette='lancet') + 
        labs(x='Cancer Type', y='No. Samples') + 
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
               sign_spldep = sign(`spldep__condition_a-median`))
    
    plts[['top_samples-dpsi_vs_spldep-scatter']] = X %>%
        ggplot(aes(x=psi__median_diff, 
                   y=`spldep__condition_a-median`, 
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
        geom_vline(xintercept=0, linetype='dashed')
    
    plts[['top_samples-dpsi_vs_spldep-selection']] = X %>%
        filter(psi__is_significant) %>%
        mutate(sign_combined = sprintf('sign_dpsi=%s & sign_spldep=%s',sign_dpsi,sign_spldep)) %>%
        count(sign_combined) %>% 
        drop_na() %>%
        ggbarplot(x='cancer_type', y='n', fill='sign_combined', palette='jco', 
                  color=FALSE, label=TRUE, position=position_dodge(0.9)) +
        theme_pubr(x.text.angle = 45, legend='right') +
        labs(x='Cancer Type', y='No. Significant Events')
   
    plts[['top_samples-dpsi_vs_spldep-candidates']] = X %>% 
        filter(cancer_type %in% cancers_oi) %>%
        filter(psi__is_significant & 
               ((sign_dpsi>0 & sign_spldep<0) | (sign_dpsi<0 & sign_spldep>0))) %>%
        arrange(-`spldep__condition_a-median`) %>%
        ggbarplot(x='event_gene', y='spldep__condition_a-median', 
                  fill='cancer_type', position=position_dodge(0.9), color=FALSE, 
                  palette=get_palette('Paired',length(unique(X[['cancer_type']])))) + 
        geom_hline(yintercept=0, linetype='dashed') +
        labs(x='Event & Gene', y='median(Spl. Dep. in PT)', fill='Cancer Type') +
        coord_flip()
    
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
        geom_text_repel(aes(label=cancer_type, color=cancer_type)) +
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
                  color=NA, label=TRUE) + guides(fill='none') + 
        labs(y='No. Samples')
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
            max.overlaps=50
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
        stat_cor(method='spearman', label.sep = '\n')
    
    return(plts)
}

# plot_enrichment = function(result, 
#                             pattern='',
#                             palette='Paired', 
#                             legend_label='Cancer Type'){
#     res = new("compareClusterResult", compareClusterResult = result)
#     plt_title = sprintf('ontology=%s | dataset=%s',unique(result[['ontology']]),unique(result[['dataset']]))
    
#     # prepare palette
#     n = length(unique(result$Cluster))
#     palette = get_palette(palette, n)
    
#     # plot
#     plts = list()
#     plts[['enrichment-dotplot']] = dotplot(res) + 
#         labs(x=legend_label, title=plt_title)
#     plts[['enrichment-cnetplot']] = cnetplot(res) + 
#             scale_fill_manual(values=palette) + 
#             labs(fill=legend_label)
#     plts[['enrichment-barplot']] = result %>% 
#         count(Description,cancer_type) %>% 
#         ggbarplot(x='Description', y='n', fill='cancer_type', 
#                   color=NA, palette=palette) + 
#         theme_pubr(x.text.angle=45) + 
#         labs(x='Term', y='Count', title=plt_title, fill=legend_label) + 
#         guides(color='none')
    
#     names(plts) = paste0(pattern,'-',names(plts))

#     return(plts)
# }


# plot_enrichments = function(gsea_result){
#     ontologies = unique(gsea_result[['ontology']])
#     datasets = unique(gsea_result[['dataset']])
    
#     plts = lapply(datasets, function(dataset_oi){
#         tmp = lapply(ontologies, function(ontology_oi){
#             result = gsea_result %>% 
#                 filter(ontology==ontology_oi & dataset==dataset_oi)
#             pattern = sprintf('%s-%s',dataset_oi,ontology_oi)
#             plt = plot_enrichment(result, pattern=pattern)
#             return(plt)
#         })
#         tmp = do.call(c,tmp)
#         return(tmp)
#     })
#     plts = do.call(c,plts)
    
#     return(plts)
# }


make_plots = function(diff_result_sample, diff_result_response, spldep, metadata, diff_result_sample_raw, psi_ccle){
    plts = list(
        plot_top_candidates_sample_type(diff_result_sample),
        plot_spldeps_pt(diff_result_sample, spldep, metadata),
        plot_ccle_vs_tcga(diff_result_sample_raw, psi_ccle),
        plot_top_candidates_response(diff_result_response)
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
    # top candidates sample type
    save_plt(plts, 'top_samples-sample_counts_cancer', '.pdf', figs_dir, width=6, height=4)
    save_plt(plts, 'top_samples-dpsi_vs_spldep-scatter', '.pdf', figs_dir, width=10, height=10)
    save_plt(plts, 'top_samples-dpsi_vs_spldep-selection', '.pdf', figs_dir, width=10, height=6)
    save_plt(plts, 'top_samples-dpsi_vs_spldep-candidates', '.pdf', figs_dir, width=10, height=10)
    
    # splicing dependencies in primary tumors
    ## PANCAN
    save_plt(plts, 'spldeps_pt-pancan-pca', '.pdf', figs_dir, width=6, height=6)
    ## BRCA
    save_plt(plts, 'spldeps_pt-brca_subtypes-counts', '.pdf', figs_dir, width=6, height=6)
    save_plt(plts, 'spldeps_pt-brca_subtypes-bar', '.pdf', figs_dir, width=10, height=15)
    save_plt(plts, 'spldeps_pt-brca_subtypes-pca', '.pdf', figs_dir, width=6, height=6)
    save_plt(plts, 'spldeps_pt-brca_subtypes-specificity', '.pdf', figs_dir, width=6, height=6)
    
    # top candidates drug response
    save_plt(plts, 'top_response-sample_counts_cancer', '.pdf', figs_dir, width=6, height=4)
    save_plt(plts, 'top_response-sample_counts_treatment', '.pdf', figs_dir, width=8, height=4)
    save_plt(plts, 'top_response-diff_psi-volcanos', '.pdf', figs_dir, width=6, height=2.5)
    save_plt(plts, 'top_response-diff_spldep-volcanos', '.pdf', figs_dir, width=6, height=2.5)
    save_plt(plts, 'top_response-diff_spldep_vs_psi-scatters', '.pdf', figs_dir, width=6, height=2.5)
    save_plt(plts, 'top_response-candidates', '.pdf', figs_dir, width=6, height=6)

    # CCLE vs TCGA
    save_plt(plts, 'ccle_vs_tcga-std_ccle_vs_dpsi', '.pdf', figs_dir, width=10, height=10)

    
    # enrichments
    ## dotplots
#     plts_oi = grep('enrichment-dotplot',names(plts),value=TRUE)
#     lapply(plts_oi, function(plt_oi){
#         save_plt(plts, plt_oi, '.pdf', figs_dir, width=10, height=15)
#     })
#     ## barplots
#     plts_oi = grep('enrichment-barplot',names(plts),value=TRUE)
#     lapply(plts_oi, function(plt_oi){
#         save_plt(plts, plt_oi, '.pdf', figs_dir, width=10, height=10)
#     })
#     ## cnetplots
#     plts_oi = grep('enrichment-cnetplot',names(plts),value=TRUE)
#     lapply(plts_oi, function(plt_oi){
#         save_plt(plts, plt_oi, '.png', figs_dir, width=15, height=15)
#     })
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
    figs_dir = args$figs_dir
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    diff_result_sample = read_tsv(diff_result_sample_file)
    diff_result_response = read_tsv(diff_result_response_file)
    annot = read_tsv(annotation_file)
#     gsea_result = read_tsv(gsea_result_file) %>%
#         mutate(Cluster = as.factor(cancer_type))
    psi_ccle = read_tsv(psi_ccle_file)
    selected_events = readLines(selected_events_file)
    
    # prep results differential analyses
    diff_result_sample_raw = diff_result_sample
    diff_result_sample = prep_diff_result(diff_result_sample, annot, selected_events)
    diff_result_response = prep_diff_result(diff_result_response, annot, selected_events)
    
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
        
    # plot
    plts = make_plots(diff_result_sample, diff_result_response, spldep, metadata, diff_result_sample_raw, psi_ccle)

    # save
    save_plots(plts, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}