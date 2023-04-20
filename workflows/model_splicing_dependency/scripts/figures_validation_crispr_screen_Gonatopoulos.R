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

require(optparse)
require(tidyverse)
require(readxl)
require(ggrepel)
require(ggpubr)
require(cowplot)
require(extrafont)

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

CELL_TYPES = data.frame(
    cell_line = c("RPE1"),
    sampleID = c("SRR2925168")
)

COMPARISONS_OI = c("0d-vs-18d")

# formatting
PAL_SINGLE_ACCENT = "orange"
PAL_DUAL = "#0C0A3E"
LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"


# Development
# -----------
# RAW_DIR = file.path(ROOT,"data","raw")
# PREP_DIR = file.path(ROOT,"data","prep")
# RESULTS_DIR = file.path(ROOT,"results","model_splicing_dependency")
# crispr_file = file.path(RAW_DIR,"articles","Gonatopoulos-Pournatzis2020","exon_targeting_library_scores.xlsx")
# event_info_file = file.path(RAW_DIR,"VastDB","EVENT_INFO-hg38_noseqs.tsv")
# psi_file = file.path(RAW_DIR,'articles','Hart2015','vast_out','PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz')
# tpm_file = file.path(RAW_DIR,'articles','Hart2015','vast_out','TPM-hg38-12.tab.gz')
# selected_events_file = file.path(RESULTS_DIR,"files","selected_models-EX.txt")
# spldep_file = file.path(RESULTS_DIR,'files','Hart2015','splicing_dependency-EX','mean.tsv.gz')

# figs_dir = file.path(RESULTS_DIR,'figures','validation_crispr_screen','Gonatopoulos-Pournatzis2020')


##### FUNCTIONS #####
plot_summary = function(crispr, avail_events){
    X = crispr %>% 
        filter(EVENT %in% avail_events) %>%
        drop_na(fitness_score)
    n_events = X %>% distinct(EVENT) %>% nrow()
    
    plts = list()
    plts[["summary-scatters"]] = X %>% 
        ggscatter(x="event_gene", y="fitness_score", palette=PAL_DUAL,
                  size=1, color="cell_line") + 
        facet_wrap(~comparison) + 
        geom_hline(yintercept=0, linetype="dashed", size=LINE_SIZE) + 
        labs(x=sprintf("Event & Gene (n=%s)", n_events), y="Obs. Delta Cell Prolif.", color="Cell Line") + 
        coord_flip() +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    return(plts)
}


plot_predictions = function(crispr, avail_events, harm){
    X = crispr %>% 
        filter(EVENT %in% avail_events) %>% 
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
        labs(x="Harm Score", y="Obs. Delta Cell Prolif.") +
        geom_smooth(method="lm", linetype="dashed", color="black", size=LINE_SIZE)
    
    return(plts)
}


make_plots = function(crispr, avail_events, spldep){
    plts = list(
        plot_summary(crispr, avail_events),
        plot_predictions(crispr, avail_events, spldep)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(crispr, avail_events, harm){
    crispr_screen = crispr %>% 
        filter(EVENT %in% avail_events) %>%
        drop_na(fitness_score)
    
    figdata = list(
        "model_validation" = list(
            "crispr_screen" = crispr_screen,
            "spotter_results" = harm
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
    save_plt(plts, "summary-scatters", ".pdf", figs_dir, width=5.75, height=8)
    save_plt(plts, "predictions-scatters", ".pdf", figs_dir, width=4.45, height=4.45)
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
        make_option("--crispr_file", type="character"),
        make_option("--event_info_file", type="character"),
        make_option("--psi_file", type="character"),
        make_option("--tpm_file", type="character"),
        make_option("--selected_events_file", type="character"),
        make_option("--spldep_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    crispr_file = args[["crispr_file"]]
    event_info_file = args[["event_info_file"]]
    psi_file = args[["psi_file"]]
    tpm_file = args[["tpm_file"]]
    selected_events_file = args[["selected_events_file"]]
    spldep_file = args[["spldep_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    crispr = read_xlsx(crispr_file, skip = 1)
    selected_events = readLines(selected_events_file)
    spldep = read_tsv(spldep_file)
    psi = read_tsv(psi_file)
    tpm = read_tsv(tpm_file)
    event_info = read_tsv(event_info_file) %>% distinct(GENE,EVENT)
    
    # preprocess
    crispr = crispr %>%
        dplyr::rename(EVENT = Event, GENE = Gene) %>%
        mutate(
            fitness_score = as.numeric(ess.RPE1.log2FC),
            comparison = "0d-vs-18d",
            cell_line = "RPE1",
            event_gene = paste0(EVENT,"_",GENE)
        )
    
    # subset experiments of interest
    crispr = crispr %>% filter(comparison %in% COMPARISONS_OI)
    
    # compute harm score
    harm = spldep %>% 
        pivot_longer(-index, names_to="sampleID", values_to="spldep") %>% 
        left_join(psi %>% 
                      pivot_longer(-EVENT, names_to="sampleID", values_to="psi_ctl"), 
                  by=c("index"="EVENT","sampleID")) %>%  
        left_join(tpm %>% 
                      dplyr::select(-ID) %>%
                      pivot_longer(-NAME, names_to="sampleID", values_to="tpm_ctl") %>%
                      left_join(event_info, by=c("NAME"="GENE")), 
                  by=c("index"="EVENT","sampleID")) %>%  
        # the amount of change upon cutting out the EVENT
        mutate(harm_score=spldep*(-psi_ctl),
               sign_harm=(-1)*harm_score) %>%
        left_join(CELL_TYPES, by="sampleID")

    # available events
    avail_events = intersect(
        crispr %>% filter(EVENT %in% selected_events) %>% drop_na(fitness_score) %>% pull(EVENT),
        harm %>% filter(index %in% selected_events) %>% drop_na(sign_harm, cell_line) %>% pull(index)
    )
    
    # plot
    plts = make_plots(crispr, avail_events, harm)

    # make figdata
    figdata = make_figdata(crispr, avail_events, harm)
    
    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
