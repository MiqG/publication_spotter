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
# PREP_DIR = file.path(ROOT,"data","prep")
# RAW_DIR = file.path(ROOT,"data","raw")
# RESULTS_DIR = file.path(ROOT,"results","model_splicing_dependency")
# models_achilles_file = file.path(RESULTS_DIR,"files","achilles","models_gene_dependency-EX","model_summaries.tsv.gz")
# models_file = file.path(RESULTS_DIR,"files","models_gene_dependency-EX","model_summaries.tsv.gz")
# cancer_events_file = file.path(ROOT,"support","cancer_events.tsv")
# figs_dir = file.path(RESULTS_DIR,"figures","validation_achilles")


##### FUNCTIONS #####
plot_model_selection = function(models, models_achilles, cancer_events){
    plts = list()         
    
    ctl_pos_events = cancer_events %>% pull(EVENT) %>% unique()
   
    X = models %>%
        bind_rows(models_achilles) %>%
        mutate(is_ctl_pos = EVENT %in% ctl_pos_events,
               lab = sprintf("%s & %s", is_ctl_pos, model_type))
    
    comparisons = list(
        c("TRUE & KO (CRISPR/Cas9)","TRUE & KD (shRNA)"),
        c("FALSE & KO (CRISPR/Cas9)","TRUE & KO (CRISPR/Cas9)"),
        c("FALSE & KD (shRNA)","TRUE & KD (shRNA)")    
    )
    
    plts[["model_sel-lr_pvalue_ctl_pos_all"]] = X %>%
        ggplot(aes(x=lab, y=lr_pvalue)) +
        geom_violin(aes(fill=is_ctl_pos, color=model_type), trim=TRUE) +
        geom_boxplot(width=0.1, outlier.size=0.1) +
        fill_palette(PAL_DUAL) +
        color_palette(c("white","black")) +
        theme_pubr(x.text.angle=30) +
        stat_compare_means(method="wilcox.test", comparisons=comparisons, 
                           size=FONT_SIZE, family=FONT_FAMILY) + 
        labs(x="Is Positive Control & Dependency Type", y="LR Test p-value")
    
    
    plts[["model_sel-lr_pvalue_ctl_pos"]] = models_achilles %>%
        mutate(is_ctl_pos = EVENT %in% ctl_pos_events) %>%
        ggviolin(x="is_ctl_pos", y="lr_pvalue", trim = TRUE, 
                 fill="is_ctl_pos", color=NA, palette = PAL_DUAL) + 
        geom_boxplot(width=0.1) +
        stat_compare_means(method="wilcox.test", family=FONT_FAMILY, size=FONT_SIZE) +
        guides(fill="none") + 
        labs(x="Is Positive Control", y="LR Test p-value")
    
    models_achilles %>%
        mutate(is_ctl_pos = EVENT %in% ctl_pos_events) %>%
        group_by(is_ctl_pos) %>%
        summarize(median(lr_pvalue, na.rm=TRUE))
    
    models %>%
        mutate(is_ctl_pos = EVENT %in% ctl_pos_events) %>%
        group_by(is_ctl_pos) %>%
        summarize(median(lr_pvalue, na.rm=TRUE))
    
    return(plts)
}


make_plots = function(models, models_achilles, cancer_events){
    plts = list(
        plot_model_selection(models, models_achilles, cancer_events)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(models_achilles, cancer_events){
    figdata = list(
        "model_selection" = list(
            "model_summaries" = models_achilles,
            "cancer_events" = cancer_events,
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
    save_plt(plts, "model_sel-lr_pvalue_ctl_pos_all", ".pdf", figs_dir, width=5, height=7)
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
    models_achilles = read_tsv(models_achilles_file)
    cancer_events = read_tsv(cancer_events_file)
    
    # add info to models
    models = models %>% 
        mutate(event_gene = paste0(EVENT,"_",GENE),
               event_type = gsub("Hsa","",gsub("[^a-zA-Z]", "",EVENT))) %>% 
        dplyr::select(-c(event_mean,event_std,gene_mean,gene_std)) %>%
        mutate(model_type="KD (shRNA)")
    
    models_achilles = models_achilles %>% 
        mutate(event_gene = paste0(EVENT,"_",GENE),
               event_type = gsub("Hsa","",gsub("[^a-zA-Z]", "",EVENT))) %>% 
        dplyr::select(-c(event_mean,event_std,gene_mean,gene_std)) %>%
        mutate(model_type="KO (CRISPR/Cas9)")
    
    # plot
    plts = make_plots(models, models_achilles, cancer_events)

    # make figdata
    figdata = make_figdata(models_achilles, cancer_events)
    
    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
