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

require(tidyverse)
require(ggpubr)
require(cowplot)
require(extrafont)

ROOT = here::here()
source(file.path(ROOT,"src","R","utils.R"))


PAL_SINGLE_ACCENT = "orange"
PAL_SINGLE_LIGHT = "#6AC2BF"
PAL_SINGLE_DARK = "#716454"
PAL_DUAL = c(PAL_SINGLE_DARK, PAL_SINGLE_LIGHT)
LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"


# Development
# -----------
PREP_DIR = file.path(ROOT,"data","prep")
RAW_DIR = file.path(ROOT,"data","raw")
RESULTS_DIR = file.path(ROOT,"results","model_splicing_dependency")
models_file = file.path(RESULTS_DIR,"files","achilles","models_gene_dependency-EX","model_summaries.tsv.gz")
cancer_events_file = file.path(ROOT,"support","cancer_events.tsv")
figs_dir = file.path(RESULTS_DIR,"figures","validation_achilles")


# ccle_stats_file = file.path(PREP_DIR,"stats","CCLE.tsv.gz")
# rnai_file = file.path(PREP_DIR,"demeter2","CCLE.tsv.gz")
# genexpr_file = file.path(PREP_DIR,"genexpr_tpm","CCLE.tsv.gz")
# splicing_file = file.path(PREP_DIR,"event_psi","CCLE-EX.tsv.gz")
# spldep_file = file.path(RESULTS_DIR,"files","splicing_dependency-EX","mean.tsv.gz")
# msigdb_dir = file.path(ROOT,"data","raw","MSigDB","msigdb_v7.4","msigdb_v7.4_files_to_download_locally","msigdb_v7.4_GMTs")
# protein_impact_file = file.path(ROOT,"data","raw","VastDB","PROT_IMPACT-hg38-v3.tab.gz")
# gene_mut_freq_file = file.path(ROOT,"data","prep","gene_mutation_freq","CCLE.tsv.gz")
# event_mut_freq_file = file.path(ROOT,"data","prep","event_mutation_freq","CCLE-EX.tsv.gz")
# indices_file = file.path(RESULTS_DIR,"files","correlation_spldep_indices-EX.tsv.gz")
# metadata_file = file.path(PREP_DIR,"metadata","CCLE.tsv.gz")
# event_info_file = file.path(RAW_DIR,"VastDB","EVENT_INFO-hg38_noseqs.tsv")

##### FUNCTIONS #####
plot_model_selection = function(models, cancer_events){
    plts = list()         
    
    ctl_pos_events = cancer_events %>% pull(EVENT) %>% unique()
    plts[["model_sel-lr_pvalue_ctl_pos"]] = models %>%
        mutate(is_ctl_pos = EVENT %in% ctl_pos_events) %>%
        ggviolin(x="is_ctl_pos", y="lr_pvalue", trim = TRUE, 
                 fill="is_ctl_pos", color=NA, palette = c("grey","darkgreen")) + 
        geom_boxplot(width=0.1) +
        stat_compare_means(method="wilcox.test", size=FONT_SIZE, family=FONT_FAMILY) +
        guides(fill="none") + 
        labs(x="Is Positive Control", y="LR Test p-value")
    
    return(plts)
}


make_plots = function(models, cancer_events){
    plts = list(
        plot_model_selection(models, cancer_events)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(models, cancer_events){
    figdata = list(
        "model_selection" = list(
            "model_summaries"= models,
            "cancer_events" = cancer_events,
        )
    )
    return(figdata)
}


save_plt = function(plts, plt_name, extension=".pdf", 
                    directory="", dpi=350, format=TRUE,
                    width = par("din")
                    3[1], height = par("din")[2]){
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
    save_plt(plts, "model_sel-lr_pvalue_ctl_pos", ".pdf", figs_dir, width=4, height=4)
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
    cancer_events_file = args$cancer_events_file
    figs_dir = args$figs_dir
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    models = read_tsv(models_file)
    cancer_events = read_tsv(cancer_events_file)
    
    # add info to models
    models = models %>% 
        mutate(event_gene = paste0(EVENT,"_",GENE),
               event_type = gsub("Hsa","",gsub("[^a-zA-Z]", "",EVENT))) %>% 
        dplyr::select(-c(event_mean,event_std,gene_mean,gene_std))
    
    # plot
    plts = make_plots(models, cancer_events)

    # make figdata
    figdata = make_figdata(models, cancer_events)
    
    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
