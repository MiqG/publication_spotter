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

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

# variables
THRESH_FDR = 0.05
THRESH_MEDIAN_DIFF = 5

# Development
# -----------
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# RESULTS_DIR = file.path(ROOT,'results','splicing_dependency_tcga')
# diff_result_file = file.path(RESULTS_DIR,'files','PANCAN','mannwhitneyu-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz')
# annotation_file = file.path(RAW_DIR,'VastDB','EVENT_INFO-hg38_noseqs.tsv')
# gsea_result_file = file.path(RESULTS_DIR,'files','PANCAN','gsea-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz')
# figs_dir = file.path(RESULTS_DIR,'figures','diffspldep_sample_type')


##### FUNCTIONS #####
plot_differential_analyses = function(diff_result){
    X = diff_result
    
    plts = list()
    plts[['diff-psi-volcanos']] = X %>%
        ggplot(aes(x=psi__median_diff, 
                   y=psi__log10_padj, 
                   color=psi__is_significant)) +
        geom_scattermore(
            pixels = c(1000,1000), 
            pointsize = 2) +
        facet_wrap(~cancer_type, ncol=4) +
        theme_pubr(border=TRUE) +
        guides(color="none") +
        labs(x='Delta PSI', y='-log10(FDR)') +
        ggpubr::color_palette(c("black","orange")) +
        geom_text(
            aes(x=x,y=y,label=n),
            X %>% 
            filter(psi__is_significant) %>% 
            mutate(is_pos = psi__median_diff>0) %>%
            count(cancer_type, psi__is_significant, is_pos) %>%
            mutate(x=ifelse(is_pos, -60, 60),
                   y=50),
            color='black'
        )
    
    plts[['diff-spldep-volcanos']] = X %>%
        ggplot(aes(x=spldep__median_diff, 
                   y=spldep__log10_padj, 
                   color=spldep__is_significant)) +
        geom_scattermore(
            pixels = c(1000,1000), 
            pointsize = 2) +
        facet_wrap(~cancer_type, ncol=4) +
        theme_pubr(border=TRUE) +
        guides(color="none") +
        labs(x='Delta Splicing Dep.', y='-log10(FDR)') +
        ggpubr::color_palette(c("black","orange")) +
        geom_text(
            aes(x=x,y=y,label=n),
            X %>% 
            filter(spldep__is_significant) %>% 
            mutate(is_pos = spldep__median_diff>0) %>%
            count(cancer_type, spldep__is_significant, is_pos) %>%
            mutate(x=ifelse(is_pos, -60, 60),
                   y=50),
            color='black'
        )
    
    plts[['diff-spldep_vs_psi-scatters']] = X %>%
        ggplot(aes(x=spldep__median_diff,
                   y=psi__median_diff)) +
        geom_scattermore(
            pixels = c(1000,1000), 
            pointsize = 2) +
        facet_wrap(~cancer_type, ncol=4, scales='free') +
        theme_pubr(border=TRUE) +
        stat_cor(method='spearman') +
        labs(x='Delta Splicing Dep.', y='Delta PSI')
    
    
    plts[['diff-spldep_vs_psi-corrs_hist']] = X %>%
        group_by(cancer_type) %>%
        summarize(corr=cor(spldep__median_diff, psi__median_diff, 
                        use='pairwise.complete.obs', method='spearman')) %>%
        gghistogram(x='corr', fill='lightblue', color=NA, bins=50) +
        xlim(-1,1) +
        labs(x='Spearman Correlation', y='Count')
    
    
    return(plts)
}


plot_top_candidates = function(diff_result){
    X = diff_result %>% filter(spldep__is_significant)
    
    plts = list()
    plts[['top-barplots']] = X %>% 
        group_by(cancer_type) %>% 
        slice_max(order_by=abs(spldep__median_diff), n=5) %>% 
        ggbarplot(x='event_gene', y='spldep__median_diff') + 
        facet_wrap(~cancer_type, scales='free', ncol=2) + 
        theme_pubr(x.text.angle=45, border=TRUE) + 
        labs(x='Event & Gene', y='Delta Splicing Dep.')
    return(plts)
}


plot_enrichment = function(result, 
                            pattern='',
                            palette='Paired', 
                            legend_label='Cancer Type'){
    res = new("compareClusterResult", compareClusterResult = result)
    plt_title = sprintf('ontology=%s | dataset=%s',unique(result[['ontology']]),unique(result[['dataset']]))
    
    # prepare palette
    n = length(unique(result$Cluster))
    palette = get_palette(palette, n)
    
    # plot
    plts = list()
    plts[['enrichment-dotplot']] = dotplot(res) + 
        labs(x=legend_label, title=plt_title)
    plts[['enrichment-cnetplot']] = cnetplot(res) + 
            scale_fill_manual(values=palette) + 
            labs(fill=legend_label)
    plts[['enrichment-barplot']] = result %>% 
        count(Description,cancer_type) %>% 
        ggbarplot(x='Description', y='n', fill='cancer_type', 
                  color=NA, palette=palette) + 
        theme_pubr(x.text.angle=45) + 
        labs(x='Term', y='Count', title=plt_title, fill=legend_label) + 
        guides(color='none')
    
    names(plts) = paste0(pattern,'-',names(plts))

    return(plts)
}


plot_enrichments = function(gsea_result){
    ontologies = unique(gsea_result[['ontology']])
    datasets = unique(gsea_result[['dataset']])
    
    plts = lapply(datasets, function(dataset_oi){
        tmp = lapply(ontologies, function(ontology_oi){
            result = gsea_result %>% 
                filter(ontology==ontology_oi & dataset==dataset_oi)
            pattern = sprintf('%s-%s',dataset_oi,ontology_oi)
            plt = plot_enrichment(result, pattern=pattern)
            return(plt)
        })
        tmp = do.call(c,tmp)
        return(tmp)
    })
    plts = do.call(c,plts)
    
    return(plts)
}


make_plots = function(diff_result, gsea_result){
    plts = list(
        plot_differential_analyses(diff_result),
        plot_top_candidates(diff_result),
        plot_enrichments(gsea_result)
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
    # differential analyses
    save_plt(plts, 'diff-psi-volcanos', '.pdf', figs_dir, width=10, height=10)
    save_plt(plts, 'diff-spldep-volcanos', '.pdf', figs_dir, width=10, height=10)
    save_plt(plts, 'diff-spldep_vs_psi-scatters', '.pdf', figs_dir, width=10, height=10)
    save_plt(plts, 'diff-spldep_vs_psi-corrs_hist', '.pdf', figs_dir, width=5, height=5)
    
    # top candidates
    save_plt(plts, 'top-barplots', '.pdf', figs_dir, width=10, height=15)
    
    # enrichments
    ## dotplots
    plts_oi = grep('enrichment-dotplot',names(plts),value=TRUE)
    lapply(plts_oi, function(plt_oi){
        save_plt(plts, plt_oi, '.pdf', figs_dir, width=10, height=15)
    })
    ## barplots
    plts_oi = grep('enrichment-barplot',names(plts),value=TRUE)
    lapply(plts_oi, function(plt_oi){
        save_plt(plts, plt_oi, '.pdf', figs_dir, width=10, height=10)
    })
    ## cnetplots
    plts_oi = grep('enrichment-cnetplot',names(plts),value=TRUE)
    lapply(plts_oi, function(plt_oi){
        save_plt(plts, plt_oi, '.png', figs_dir, width=15, height=15)
    })
}


main = function(){
    args = getParsedArgs()
    diff_result_file = args$diff_result_file
    gsea_result_file = args$gsea_result_file
    annotation_file = args$annotation_file
    figs_dir = args$figs_dir
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    diff_result = read_tsv(diff_result_file)
    annot = read_tsv(annotation_file)
    gsea_result = read_tsv(gsea_result_file) %>%
        mutate(Cluster = as.factor(cancer_type))
    
    # prep
    diff_result = diff_result %>%
        rename_all(recode, index = "EVENT") %>%
        mutate(event_type = gsub('Hsa','',gsub("[^a-zA-Z]", "",EVENT))) %>%
        left_join(annot[,c('EVENT','GENE')], by='EVENT') %>%
        mutate(event_gene = paste0(EVENT,'_',GENE)) %>%
        group_by(cancer_type) %>%
        mutate(psi__padj = p.adjust(psi__pvalue, 'fdr'),
               spldep__padj = p.adjust(spldep__pvalue, 'fdr'),
               psi__log10_padj = -log10(psi__padj),
               spldep__log10_padj = -log10(spldep__padj)) %>%
        mutate(psi__is_significant = psi__padj<THRESH_FDR & 
                                     abs(psi__median_diff)>THRESH_MEDIAN_DIFF,
               spldep__is_significant = spldep__padj<THRESH_FDR &
                                        psi__is_significant)
        
    
    # plot
    plts = make_plots(diff_result, gsea_result)

    # save
    save_plots(plts, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}