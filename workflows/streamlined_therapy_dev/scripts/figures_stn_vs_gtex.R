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

require(optparse)
require(tidyverse)
require(ggpubr)
require(cowplot)
require(scattermore)
require(ggrepel)
require(extrafont)

# variables
THRESH_FDR = 0.05
THRESH_MEDIAN_DIFF = 5
MIN_SAMPLES = 10

CANCERS_OI = c(
    'BLCA','BRCA','COAD','HNSC','KICH','KIRC','KIRP',
    'LIHC','LUAD','LUSC','PRAD','THCA','UCEC',
    'BRCA_brca_Basal','BRCA_brca_Her2','BRCA_brca_NotBasal',
    'UCEC_ucec_CN_high','UCEC_ucec_CN_low','UCEC_ucec_MSI','UCEC_ucec_POLE'
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
ROOT = here::here()
RAW_DIR = file.path(ROOT,'data','raw')
PREP_DIR = file.path(ROOT,'data','prep')
RESULTS_DIR = file.path(ROOT,'results','streamlined_therapy_dev')
annotation_file = file.path(RAW_DIR,'VastDB','EVENT_INFO-hg38_noseqs.tsv')
selected_events_file = file.path(ROOT,'results','model_splicing_dependency','files','selected_models-EX.txt')
diff_result_sample_file = file.path(RESULTS_DIR,'files','PANCAN','mannwhitneyu-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz')
gtex_file = file.path(RAW_DIR,"inhouse","Sodaei","event_psi","breast-mammary_tissue-EX.tsv.gz")
figs_dir = file.path(RESULTS_DIR,'figures','stn_vs_gtex')

##### FUNCTIONS #####
plot_stn_tcga_vs_gtex = function(diff_result_sample, correlations){
    plts = list()
    
    X = correlations %>%
        pivot_longer(
            -c(cancer_type,n_events), 
            names_to="correlation_method", values_to="correlation_coef"
        )
    
    plts[["stn_tcga_vs_gtex-correlations-violin"]] = X %>%
        ggviolin(x="correlation_method", y="correlation_coef", fill=PAL_SINGLE_LIGHT, color=NA, trim=TRUE) +
        geom_boxplot(width=0.1, outlier.size=0.1) +
        geom_text(
            aes(label=label, y=1), 
            . %>% count(correlation_method) %>% mutate(label=sprintf("n=%s",n)), 
            family=FONT_FAMILY, size=FONT_SIZE
        ) +
        geom_text(aes(label=cancer_type), . %>% filter(cancer_type=="BRCA"), family=FONT_FAMILY, size=FONT_SIZE) +
        facet_wrap(~correlation_method, scales="free") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Correlation Method", y="Correlation Coef.\nSTN TCGA vs Breast GTEx")
    
    plts[["stn_tcga_vs_gtex-delta_psi-scatter"]] = diff_result_sample %>%
        filter(psi__is_significant) %>% 
        drop_na(psi__median_diff, gtex_median_diff) %>%
        ggplot(aes(x=psi__median_diff, y=gtex_median_diff)) +
        geom_scattermore(pixels = c(1000,1000), pointsize=5, alpha=0.5, color="brown") +
        stat_cor(method="spearman", size=FONT_SIZE, family=FONT_FAMILY) +
        theme_pubr() +
        facet_wrap(~cancer_type) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="DeltaPSI PT TCGA vs STN TCGA", y="DeltaPSI PT TCGA vs Breast GTEx")
    
    return(plts)
}

make_plots = function(diff_result_sample, correlations){
    plts = list(
        plot_stn_tcga_vs_gtex(diff_result_sample, correlations),
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(diff_result_sample, correlations){
    figdata = list(
        "stn_vs_gtex" = list(
            "diff_result_sample" = diff_result_sample %>%
                            filter(psi__is_significant) %>% 
                            drop_na(psi__median_diff, gtex_median_diff),
            "correlations" = correlations
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
    save_plt(plts, 'stn_tcga_vs_gtex-correlations-violin', '.pdf', figs_dir, width=8, height=5)
    save_plt(plts, 'stn_tcga_vs_gtex-delta_psi-scatter', '.pdf', figs_dir, width=15, height=15)
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
        make_option("--annotation_file", type="character"),
        make_option("--selected_events_file", type="character"),
        make_option("--diff_result_sample_file", type="character"),
        make_option("--gtex_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    annotation_file = args[["annotation_file"]]
    selected_events_file = args[["selected_events_file"]]
    diff_result_sample_file = args[["diff_result_sample_file"]]
    gtex_file = args[["gtex_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    diff_result_sample = read_tsv(diff_result_sample_file)
    gtex = read_tsv(gtex_file)
    annot = read_tsv(annotation_file)
    selected_events = readLines(selected_events_file)
    
    # correlate median psi of different TCGA tissues with GTEx
    gtex_median = gtex %>%
        column_to_rownames("EVENT") %>%
        as.matrix() %>%
        apply(.,1,median,na.rm=TRUE) %>%
        enframe("EVENT","gtex_psi_median")
    
    diff_result_sample = diff_result_sample %>%
        filter(cancer_type %in% CANCERS_OI) %>%
        left_join(gtex_median, by="EVENT") %>%
        mutate(
            gtex_median_diff = `psi__condition_a-median` - gtex_psi_median,
            psi__is_significant = psi__padj<THRESH_FDR & abs(psi__median_diff)>THRESH_MEDIAN_DIFF
        )
    
    correlations = diff_result_sample %>%
        drop_na(`psi__condition_b-median`, gtex_psi_median) %>%
        group_by(cancer_type) %>%
        summarize(
            correlation_pearson = cor(`psi__condition_b-median`, gtex_psi_median, method="pearson"),
            correlation_spearman = cor(`psi__condition_b-median`, gtex_psi_median, method="spearman"),
            n_events = n()
        ) %>%
        ungroup()
    
    # plot
    plts = make_plots(diff_result_sample, correlations)
    
    # make figdata
    figdata = make_figdata(diff_result_sample, correlations)

    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
