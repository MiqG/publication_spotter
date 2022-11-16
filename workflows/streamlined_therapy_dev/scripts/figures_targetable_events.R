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
require(extrafont)

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

# variables
THRESH_FDR = 0.05
THRESH_MEDIAN_DIFF = 5
MIN_SAMPLES = 10

CANCERS_OI = c(
    'BRCA','COAD','HNSC','KICH','KIRC','KIRP',
    'LIHC','LUAD','LUSC','PRAD','THCA','UCEC',
    'BRCA_brca_Basal','BRCA_brca_Her2','BRCA_brca_NotBasal',
    'UCEC_ucec_CN_high','UCEC_ucec_CN_low','UCEC_ucec_MSI','UCEC_ucec_POLE'
)

VALIDATED_EXONS = c(
    "HsaEX0070392_VLDLR",
    "HsaEX0052877_RCC1",
    "HsaEX0049558_PPP1R12A",
    "HsaEX0044398_NUP85",
    "HsaEX0071941_YAP1",
    "HsaEX0050345_PRPF18",
    "HsaEX0026116_FNBP1"
)

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
# RESULTS_DIR = file.path(ROOT,'results','streamlined_therapy_dev')
# annotation_file = file.path(RAW_DIR,'VastDB','EVENT_INFO-hg38_noseqs.tsv')
# selected_events_file = file.path(ROOT,'results','model_splicing_dependency','files','selected_models-EX.txt')
# spldep_stats_file = file.path(RESULTS_DIR,'files','PANCAN','summary_splicing_dependency-EX.tsv.gz')
# spldep_stats_subtypes_file = file.path(RESULTS_DIR,'files','PANCAN_subtypes','summary_splicing_dependency-EX.tsv.gz')
# diff_result_sample_file = file.path(RESULTS_DIR,'files','PANCAN','mannwhitneyu-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz')
# diff_result_subtypes_file = file.path(RESULTS_DIR,'files','PANCAN_subtypes','mannwhitneyu-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz')
# cancer_events_file = file.path(ROOT,"support","cancer_events.tsv")
# ascanceratlas_file = file.path(RAW_DIR,"ASCancerAtlas","CASE_all-VastDB_mapped.tsv.gz")

# figs_dir = file.path(RESULTS_DIR,'figures','targetable_events')

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


plot_top_candidates_sample_type = function(diff_result, patt=''){
    
    plts = list()

    # how many samples do we have for each cancer and sample type?
    X = diff_result %>% 
        mutate(PT = `psi__condition_a-n_present` + `psi__condition_a-n_nan`,
               STN = `psi__condition_b-n_present` + `psi__condition_b-n_nan`) %>% 
        distinct(cancer_type,PT,STN) %>%
        pivot_longer(c(PT,STN), names_to='sample_type',values_to='n')
    
    plts[['top_samples-sample_counts_cancer']] = X %>% 
        filter(cancer_type %in% CANCERS_OI) %>%
        ggbarplot(x='cancer_type', y='n', fill='sample_type', color=NA,
                  position=position_dodge(0.7), label=TRUE, lab.size=FONT_SIZE, 
                  lab.family=FONT_FAMILY, palette='lancet') + 
        labs(x='Cancer Type', y='N. Samples', fill="Sample Type") + 
        theme_pubr(x.text.angle = 45) + 
        geom_hline(yintercept=10, linetype='dashed', size=LINE_SIZE)
    
    # ranking with differential analyses
    X = diff_result %>%
        filter(cancer_type %in% CANCERS_OI) %>%
        mutate(sign_dpsi = sign(psi__median_diff),
               sign_spldep = sign(mean),
               event_gene = ifelse(event_gene %in% VALIDATED_EXONS, 
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
        filter(cancer_type %in% CANCERS_OI) %>%
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
        filter(cancer_type %in% CANCERS_OI) %>%
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
        filter(cancer_type %in% CANCERS_OI) %>%
        filter(psi__is_significant & 
               ((sign_dpsi>0 & sign_spldep<0) | (sign_dpsi<0 & sign_spldep>0))) %>%
        arrange(-mean) %>%
        ggbarplot(x='event_gene', y='harm_score', 
                  fill='cancer_type', position=position_dodge(0.9), color=FALSE, 
                  palette=get_palette('Paired',length(unique(X[['cancer_type']])))) + 
        geom_hline(yintercept=0, linetype='dashed', size=LINE_SIZE) +
        labs(x='Event & Gene', y='Max. Harm Score', fill='Cancer Type') +
        coord_flip()
    
    x = X %>% 
        filter(cancer_type %in% CANCERS_OI) %>%
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
        
    names(plts) = paste0(patt,names(plts))
    
    return(plts)
}


make_plots = function(diff_result_sample, diff_result_subtypes){
    plts = list(
        plot_top_candidates_sample_type(diff_result_sample),
        plot_top_candidates_sample_type(diff_result_subtypes, 'subtypes-')
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
    save_plt(plts, 'top_samples-dpsi_vs_spldep-candidates_spldep', '.pdf', figs_dir, width=5.5, height=19)
    save_plt(plts, 'top_samples-dpsi_vs_spldep-candidates_dpsi', '.pdf', figs_dir, width=5.5, height=19)
    save_plt(plts, 'top_samples-dpsi_vs_spldep-candidates_harm', '.pdf', figs_dir, width=5.5, height=19)
    save_plt(plts, 'top_samples-dpsi_vs_spldep-candidates_psi', '.pdf', figs_dir, width=18, height=30)
    
    # top candidates sample type (cancer subtypes)
    save_plt(plts, 'subtypes-top_samples-sample_counts_cancer', '.pdf', figs_dir, width=8, height=8)
    save_plt(plts, 'subtypes-top_samples-dpsi_vs_spldep-scatter', '.pdf', figs_dir, width=10, height=10)
    save_plt(plts, 'subtypes-top_samples-dpsi_vs_spldep-selection', '.pdf', figs_dir, width=10, height=6)
    save_plt(plts, 'subtypes-top_samples-dpsi_vs_spldep-candidates_spldep', '.pdf', figs_dir, width=5.5, height=8)
    save_plt(plts, 'subtypes-top_samples-dpsi_vs_spldep-candidates_dpsi', '.pdf', figs_dir, width=5.5, height=8)
    save_plt(plts, 'subtypes-top_samples-dpsi_vs_spldep-candidates_harm', '.pdf', figs_dir, width=5.5, height=8)
    save_plt(plts, 'subtypes-top_samples-dpsi_vs_spldep-candidates_psi', '.pdf', figs_dir, width=15, height=15)
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
    figs_dir = args$figs_dir
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    diff_result_sample = read_tsv(diff_result_sample_file)
    diff_result_subtypes = read_tsv(diff_result_subtypes_file)
    annot = read_tsv(annotation_file)
    selected_events = readLines(selected_events_file)
    spldep_stats = read_tsv(spldep_stats_file) %>% 
        filter(EVENT %in% selected_events)
    spldep_stats_subtypes = read_tsv(spldep_stats_subtypes_file) %>% 
        filter(EVENT %in% selected_events)
    cancer_events = read_tsv(cancer_events_file)
    ascanceratlas = read_tsv(ascanceratlas_file)
    
    # prep cancer events
    ascanceratlas = ascanceratlas %>%
        dplyr::rename(EVENT = EVENT_perf)
    
    cancer_events = cancer_events %>%
        bind_rows(ascanceratlas)
    
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
    
    # plot
    plts = make_plots(diff_result_sample, diff_result_subtypes)
    
    # make figdata
    figdata = make_figdata(diff_result_sample, diff_result_subtypes, spldep_stats)

    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}