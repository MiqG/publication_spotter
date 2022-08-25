#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# 
# Outline
# -------
# 1. Real dependencies of our selected exons in different conditions
# 2. How do our predicted splicing dependencies compare to the real ones?

require(tidyverse)
require(ggrepel)
require(ggpubr)
require(cowplot)
require(extrafont)

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

CELL_TYPES = data.frame(
    cell_line = c("HeLa", "PC9"),
    sampleID = c("SRR7946515_1", "SRR7946516_1")
)

COMPARISONS_OI = c("0d-vs-14d")

# formatting
PAL_SINGLE_ACCENT = "orange"
PAL_SINGLE_LIGHT = "#8AD5E7"#"#6AC2BF"
PAL_SINGLE_DARK = "#E6C9A8"#"#716454"
PAL_DUAL = c(PAL_SINGLE_DARK, PAL_SINGLE_LIGHT)
LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

# Development
# -----------
# RAW_DIR = file.path(ROOT,"data","raw")
# PREP_DIR = file.path(ROOT,"data","prep")
# RESULTS_DIR = file.path(ROOT,"results","model_splicing_dependency")
# crispr_file = file.path(PREP_DIR,'Thomas2020','crispr_screen.tsv.gz')
# psi_file = file.path(RAW_DIR,'articles','Thomas2020','vast_out','PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz')
# selected_events_file = file.path(RESULTS_DIR,"files","selected_models-EX.txt")
# spldep_file = file.path(RESULTS_DIR,'files','Thomas2020','splicing_dependency-EX','mean.tsv.gz')
# figs_dir = file.path(RESULTS_DIR,'figures','validation_crispr_screen')


##### FUNCTIONS #####
plot_summary = function(crispr, selected_events){
    X = crispr %>% 
        filter(EVENT %in% selected_events) %>% 
        mutate(log10_pvalue=-log10(pvalue),
               event_gene = paste0(EVENT,"_",geneName))
    n_events = X %>% distinct(EVENT) %>% nrow()
    
    plts = list()
    plts[["summary-scatters"]] = X %>% 
        ggscatter(x="event_gene", y="fitness_score", palette=PAL_DUAL,
                  size="log10_pvalue", color="cell_line") + 
        facet_wrap(~comparison) + 
        geom_hline(yintercept=0, linetype="dashed", size=LINE_SIZE) + 
        labs(x=sprintf("Event & Gene (n=%s)", n_events), y="Event Dependency", 
             color="Cell Line", size="-log10(p-value)") + 
        coord_flip() +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        scale_size(range=c(0.5,3))
    
    return(plts)
}


plot_predictions = function(crispr, selected_events, harm){
    X = crispr %>% 
        filter(EVENT %in% selected_events) %>% 
        left_join(harm, by=c("EVENT"="index","cell_line")) %>%
        drop_na(sign_harm, fitness_score) %>%
        group_by(comparison, cell_line) %>%
        mutate(cell_line = sprintf("%s (n=%s)", cell_line, n())) %>%
        ungroup() %>%
        arrange(cell_line, comparison)
    
    plts = list()
    plts[["predictions-scatters"]] = X %>% 
        ggscatter(x="sign_harm", y="fitness_score", size=1, color="cell_line", palette=PAL_DUAL) + 
        facet_wrap(~comparison+cell_line, ncol=2) + 
        stat_cor(method="pearson", size=FONT_SIZE, family=FONT_FAMILY) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY),
              aspect.ratio=1) +
        guides(color="none") + 
        labs(x="Harm Score", y="Event Dependency") +
        geom_smooth(method="lm", linetype="dashed", color="black", size=LINE_SIZE)
    
    return(plts)
}


make_plots = function(crispr, selected_events, spldep){
    plts = list(
        plot_summary(crispr, selected_events),
        plot_predictions(crispr, selected_events, spldep)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(crispr, selected_events, spldep){
    figdata = list(
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
    save_plt(plts, "summary-scatters", ".pdf", figs_dir, width=5.75, height=7)
    save_plt(plts, "predictions-scatters", ".pdf", figs_dir, width=6, height=6)
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
    figs_dir = args$figs_dir
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    crispr = read_tsv(crispr_file)
    selected_events = readLines(selected_events_file)
    spldep = read_tsv(spldep_file)
    psi = read_tsv(psi_file)
    
    # subset experiments of interest
    crispr = crispr %>% filter(comparison %in% COMPARISONS_OI)
    
    # compute harm score
    harm = spldep %>% 
        pivot_longer(-index, names_to="sampleID", values_to="spldep") %>% 
        left_join(psi %>% 
                      pivot_longer(-EVENT, names_to="sampleID", values_to="psi_ctl"), 
                  by=c("index"="EVENT","sampleID")) %>%  
        # the amount of change upon cutting out the EVENT
        mutate(harm_score=spldep*(-psi_ctl),
               sign_harm=(-1)*harm_score) %>%
        left_join(CELL_TYPES, by="sampleID")
    
    plts = make_plots(crispr, selected_events, harm)

    # make figdata
    # figdata = make_figdata(crispr, selected_events, spldep)
    
    # save
    save_plots(plts, figs_dir)
    # save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}