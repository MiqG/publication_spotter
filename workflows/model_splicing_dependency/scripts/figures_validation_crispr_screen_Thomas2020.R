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
require(ggrepel)
require(ggpubr)
require(cowplot)
require(extrafont)

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
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,"data","raw")
# PREP_DIR = file.path(ROOT,"data","prep")
# RESULTS_DIR = file.path(ROOT,"results","model_splicing_dependency")
# crispr_file = file.path(PREP_DIR,'Thomas2020','crispr_screen.tsv.gz')
# event_info_file = file.path(RAW_DIR,"VastDB","EVENT_INFO-hg38_noseqs.tsv")
# psi_file = file.path(RAW_DIR,'articles','Thomas2020','vast_out','PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz')
# tpm_file = file.path(RAW_DIR,'articles','Thomas2020','vast_out','TPM-hg38-2.tab.gz')
# selected_events_file = file.path(RESULTS_DIR,"files","selected_models-EX.txt")
# spldep_file = file.path(RESULTS_DIR,'files','Thomas2020','splicing_dependency-EX','mean.tsv.gz')
# figs_dir = file.path(RESULTS_DIR,'figures','validation_crispr_screen','Thomas2020')
# ccle_splicing_file = file.path(PREP_DIR,"event_psi","CCLE-EX.tsv.gz")
# ccle_spldep_file = file.path(RESULTS_DIR,"files","splicing_dependency-EX","mean.tsv.gz")
# figdata_gonato_spotter_file = file.path(RESULTS_DIR,'figures','validation_crispr_screen','Gonatopoulos-Pournatzis2020','figdata','model_validation','spotter_results.tsv.gz')
# figdata_gonato_crispr_file = file.path(RESULTS_DIR,'figures','validation_crispr_screen','Gonatopoulos-Pournatzis2020','figdata','model_validation','crispr_screen.tsv.gz')


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
        labs(x=sprintf("Event & Gene (n=%s)", n_events), y="Obs. Delta Cell Prolif.", 
             color="Cell Line", size="-log10(p-value)") + 
        coord_flip() +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        scale_size(range=c(0.5,3))
    
    return(plts)
}


plot_predictions = function(crispr, selected_events, harm, figdata_gonato){
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
        labs(x="Harm Score", y="Obs. Delta Cell Prolif.") +
        geom_smooth(method="lm", linetype="dashed", color="black", size=LINE_SIZE)
    
    # highlight the low dynamic range of HeLa
    plts[["predictions-scatters-low_dynrange"]] = X %>%
        #bind_rows(figdata_gonato) %>%
#         group_by(cell_line) %>%
#         mutate(
#             fitness_score = (fitness_score - mean(fitness_score)) / sd(fitness_score),
#             #sign_harm = (sign_harm - mean(sign_harm)) / sd(sign_harm)
#         ) %>%
#         ungroup() %>%
        ggplot(aes(x=sign_harm, y=fitness_score)) +
        geom_smooth(aes(color=cell_line), method="lm", linetype="dashed", size=LINE_SIZE, alpha=0.1) +
        geom_point(aes(color=cell_line), . %>% filter(!str_detect(cell_line, "HeLa")), size=1) +
        geom_point(aes(color=cell_line), . %>% filter(str_detect(cell_line, "HeLa")), size=1) +
        stat_cor(aes(color=cell_line), method="pearson", size=FONT_SIZE, family=FONT_FAMILY) +
        color_palette(PAL_DUAL) + 
        theme_pubr() +
        facet_wrap(~comparison) + 
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY), aspect.ratio=1) +
        labs(x="Harm Score", y="Obs. Delta Cell Prolif.", color="Cell Line")

    
    return(plts)
}


plot_predictions_ccle = function(crispr, selected_events, ccle_harm){
    X = crispr %>% 
        filter(EVENT %in% selected_events) %>% 
        left_join(ccle_harm, by=c("EVENT"="index","cell_line")) %>%
        drop_na(harm_score, fitness_score) %>%
        group_by(comparison, cell_line) %>%
        mutate(cell_line = sprintf("%s (n=%s)", cell_line, n())) %>%
        ungroup() %>%
        arrange(cell_line, comparison)
    
    plts = list()
    plts[["predictions_ccle-scatters"]] = X %>% 
        ggscatter(x="harm_score", y="fitness_score", size=1, color="cell_line", palette=PAL_DUAL) + 
        facet_wrap(~comparison+cell_line, ncol=2) + 
        stat_cor(method="pearson", size=FONT_SIZE, family=FONT_FAMILY) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY),
              aspect.ratio=1) +
        guides(color="none") + 
        labs(x="Harm Score", y="Obs. Delta Cell Prolif.") +
        geom_smooth(method="lm", linetype="dashed", color="black", size=LINE_SIZE)
    
    return(plts)
}


make_plots = function(crispr, selected_events, harm, ccle_harm){
    plts = list(
        plot_summary(crispr, selected_events),
        plot_predictions(crispr, selected_events, harm),
        plot_predictions_ccle(crispr, selected_events, ccle_harm)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(crispr, selected_events, harm, ccle_harm){
    # crispr screen
    crispr_screen = crispr %>% 
        filter(EVENT %in% selected_events) %>% 
        mutate(log10_pvalue=-log10(pvalue),
               event_gene = paste0(EVENT,"_",geneName))
    
    figdata = list(
        "model_validation" = list(
            "crispr_screen" = crispr_screen,
            "spotter_results" = harm,
            "spotter_results_ccle" = ccle_harm
        )
    )
    return(figdata)
}

make_source_data = function(plts){
    
    source_data = list(
        # FIGURE 2
        ## Fig. 2b left
        "fig02b_left" = plts[["predictions-scatters"]][["data"]] %>% 
            dplyr::select(event,EVENT,geneName,sign_harm,fitness_score,comparison,cell_line) %>%
            filter(cell_line=="PC9 (n=16)"),
        
        # SUPPLEMENTARY FIGURE 5
        ## Sup. Fig. 5a
        "supfig05a" = plts[["summary-scatters"]][["data"]]%>% 
            dplyr::select(event,EVENT,geneName,fitness_score,pvalue,comparison,cell_line),
        
        ## Sup. Fig. 5c
        "supfig05c" = plts[["predictions-scatters-low_dynrange"]][["data"]] %>%
            dplyr::select(event,EVENT,geneName,sign_harm,fitness_score,comparison,cell_line)
    )
    
    return(source_data)
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
    save_plt(plts, "summary-scatters", ".pdf", figs_dir, width=5.75, height=7.5)
    save_plt(plts, "predictions-scatters", ".pdf", figs_dir, width=6, height=6)
    save_plt(plts, "predictions-scatters-low_dynrange", ".pdf", figs_dir, width=3, height=6)
    save_plt(plts, "predictions_ccle-scatters", ".pdf", figs_dir, width=6, height=6)
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

save_source_data = function(source_data, dir){
    d = file.path(dir,"figdata",'source_data')
    dir.create(d, recursive=TRUE)
    lapply(names(source_data), function(nm){
        df = source_data[[nm]]
        filename = file.path(d, paste0(nm,'.tsv.gz'))
        write_tsv(df, filename)
        print(filename)
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
        make_option("--ccle_splicing_file", type="character"),
        make_option("--ccle_spldep_file", type="character"),
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
    ccle_splicing_file = args[["ccle_splicing_file"]]
    ccle_spldep_file = args[["ccle_spldep_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    crispr = read_tsv(crispr_file)
    selected_events = readLines(selected_events_file)
    spldep = read_tsv(spldep_file)
    psi = read_tsv(psi_file)
    tpm = read_tsv(tpm_file)
    event_info = read_tsv(event_info_file) %>% distinct(GENE,EVENT)
    
    ccle_spldep = read_tsv(ccle_spldep_file)
    ccle_splicing = read_tsv(ccle_splicing_file)
    
    #figdata_gonato_spotter = read_tsv(figdata_gonato_spotter_file)
    #figdata_gonato_crispr = read_tsv(figdata_gonato_crispr_file)
    
    # subset experiments of interest
    crispr = crispr %>% filter(comparison %in% COMPARISONS_OI)
    
    # subset ccle
    ccle_spldep = ccle_spldep %>% 
        filter(index %in% (crispr%>%pull(EVENT))) %>% 
        dplyr::select(c(index,`ACH-001086`,`ACH-000779`)) %>% 
        dplyr::rename(HeLa=`ACH-001086`, PC9=`ACH-000779`)
    
    ccle_splicing = ccle_splicing %>% 
        filter(EVENT %in% (crispr%>%pull(EVENT))) %>% 
        dplyr::select(c(EVENT,`ACH-001086`,`ACH-000779`)) %>% 
        dplyr::rename(HeLa=`ACH-001086`, PC9=`ACH-000779`)
    
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
    
    ccle_harm = ccle_spldep %>% 
        pivot_longer(-index, names_to="cell_line", values_to="spldep") %>% 
        left_join(
            ccle_splicing %>% 
                pivot_longer(-EVENT, names_to="cell_line", values_to="psi_ctl"), 
            by=c("index"="EVENT","cell_line")
        ) %>%  
        # the amount of change upon cutting out the EVENT
        mutate(harm_score=(-1)*(0 - psi_ctl)*spldep) 
    
    # available events
    avail_events = intersect(
        crispr %>% filter(EVENT %in% selected_events) %>% drop_na(fitness_score) %>% pull(EVENT),
        harm %>% filter(index %in% selected_events) %>% drop_na(sign_harm, cell_line) %>% pull(index)
    )
    
    # prepare figdata Gonatopoulos
    #     figdata_gonato = figdata_gonato_crispr %>%
    #         distinct(EVENT, GENE, event_gene, cell_line, fitness_score) %>%
    #         left_join(
    #             figdata_gonato_spotter %>% 
    #                 distinct(index, sign_harm, cell_line), by=c("EVENT"="index", "cell_line")
    #         )
    
    # plot
    plts = make_plots(crispr, selected_events, harm, ccle_harm)

    # make figdata
    figdata = make_figdata(crispr, selected_events, harm, ccle_harm)
    
    # make source data
    source_data = make_source_data(plts)
    
    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
    save_source_data(source_data, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}