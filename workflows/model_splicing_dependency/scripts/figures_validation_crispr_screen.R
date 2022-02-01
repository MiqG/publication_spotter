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

# Development
# -----------
# RAW_DIR = file.path(ROOT,"data","raw")
# PREP_DIR = file.path(ROOT,"data","prep")
# RESULTS_DIR = file.path(ROOT,"results","model_splicing_dependency")
# crispr_file = file.path(PREP_DIR,'Thomas2020','crispr_screen.tsv.gz')
# selected_events_file = file.path(RESULTS_DIR,"files","selected_models-EX.txt")
# spldep_file = file.path(RESULTS_DIR,'files','Thomas2020','splicing_dependency-EX','mean.tsv.gz')
# figs_dir = file.path(RESULTS_DIR,'figures','validation_crispr_screen')


##### FUNCTIONS #####
plot_summary = function(crispr, selected_events){
    X = crispr %>% 
        filter(EVENT %in% selected_events) %>% 
        mutate(log10_pvalue=-log10(pvalue),
               event_gene = paste0(EVENT,"_",geneName))
    
    plts = list()
    plts[["summary-scatters"]] = X %>% 
        ggscatter(x="event_gene", y="fitness_score", palette="nejm",
                  size="log10_pvalue", color="cell_line") + 
        facet_wrap(~comparison) + 
        geom_hline(yintercept=0, linetype="dashed") + 
        labs(x="Event & Gene", y="Fitness ScoreS", 
             color="Cell Line", size="-log10(p-value)") + 
        coord_flip() +
        theme(strip.text.x = element_text(size=6, family="Arial")) +
        scale_size(range=c(0.5,3))
    
    return(plts)
}


plot_predictions = function(crispr, selected_events, spldep){
    X = crispr %>% 
        filter(EVENT %in% selected_events) %>% 
        left_join(
            spldep %>% 
                pivot_longer(-index, names_to="sampleID", values_to="pred") %>% 
                drop_na() %>% 
                left_join(CELL_TYPES, by="sampleID"), 
                          by=c("EVENT"="index","cell_line")
        )
    
    plts = list()
    plts[["predictions-scatters"]] = X %>% 
        ggscatter(x="fitness_score", y="pred", size=1, alpha=0.5) + 
        facet_wrap(~cell_line+comparison, ncol=2) + 
        stat_cor(method="pearson", size=1, family="Arial") + 
        geom_text(aes(x=-1.35, y=-0.5, label=lab), 
                  X %>% 
                  drop_na(fitness_score,pred) %>% 
                  count(cell_line,comparison) %>% 
                  mutate(lab=paste0("n=",n)),
                  size=1, family="Arial") +
        theme(strip.text.x = element_text(size=6, family="Arial"))
    
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
                    font.x=8, font.y=8, font.legend=8,
                    font.tickslab=6, font.family="Arial")    
    }
    filename = file.path(directory,paste0(plt_name,extension))
    save_plot(filename, plt, base_width=width, base_height=height, dpi=dpi, units="cm")
}


save_plots = function(plts, figs_dir){
    save_plt(plts, "summary-scatters", ".pdf", figs_dir, width=8, height=10)
    save_plt(plts, "predictions-scatters", ".pdf", figs_dir, width=8, height=9)
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
    
    plts = make_plots(crispr, selected_events, spldep)

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