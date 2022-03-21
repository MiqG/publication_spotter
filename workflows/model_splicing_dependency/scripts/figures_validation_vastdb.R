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
# - can we predict the proliferation effects of KDs from the changes in splicing?
# - if we define exon sets from KDs, is any KD enriched in our selected events?

require(tidyverse)
require(ggpubr)
require(cowplot)
require(scattermore)
require(ggrepel)
require(extrafont)
require(tidytext)
require(clusterProfiler)
require(ComplexHeatmap)
require(ggplotify)
require(gridExtra)
require(circlize)

ROOT = here::here()
source(file.path(ROOT,"src","R","utils.R"))

# vaiables 
THRESH_DPSI = 0
RANDOM_SEED = 1234

# Development
# -----------
# RAW_DIR = file.path(ROOT,"data","raw")
# PREP_DIR = file.path(ROOT,"data","prep")
# RESULTS_DIR = file.path(ROOT,"results","model_splicing_dependency")
# psi_file = file.path(RAW_DIR,"VastDB","PSI_TABLE-hg38-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz")
# genexpr_file = file.path(RAW_DIR,"VastDB","EXPRESSION_TABLE-hg38.tab.gz")
# selected_events_file = file.path(RESULTS_DIR,"files","selected_models-EX.txt")
# event_info_file = file.path(RAW_DIR,"VastDB","EVENT_INFO-hg38_noseqs.tsv")
# figs_dir = file.path(RESULTS_DIR,'figures','validation_vastdb')


##### FUNCTIONS #####
plot_vastdb_validation = function(psi, index){
    
    plts = list()
    
    cols_order = index %>% pull(sample_clean)
    
    mat = psi %>% column_to_rownames("EVENT")
    mat = apply(mat, 1, function(x){
        x[is.na(x)] = median(x, na.rm=TRUE)
        return(x)
    })
    mat = scale(mat)
    mat = t(mat)
    mat = mat[rowSums(is.na(mat))!=ncol(mat),]
    
    col_fun = colorRamp2(c(-1, 0, 1), c("green", "white", "red"))
    
    Heatmap(mat[,cols_order], cluster_columns = FALSE, name="mat", row_km=3, show_row_names=FALSE, col=col_fun) %>%
        draw() %>%
        grid.grabExpr() %>%
        as.ggplot()
    
    return(plts)
}


make_plots = function(psi){
    plts = list(
        plot_vastdb_validation(psi)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(){
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
    # correlations of predictions
    save_plt(plts, "encore_val-thresh_vs_pearsons", ".pdf", figs_dir, width=6, height=12)
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
    genexpr = read_tsv(genexpr_file)
    psi = read_tsv(psi_file)
    selected_events = readLines(selected_events_file)
    
    psi = psi %>% filter(EVENT %in% selected_events)
    index = genexpr %>% 
        filter(NAME=="MKI67") %>% 
        dplyr::select(-ID) %>% 
        column_to_rownames("NAME") %>% 
        t() %>% as.data.frame() %>% 
        rownames_to_column("sampleID") %>% 
        filter(str_detect(sampleID,"cRPKM")) %>% 
        arrange(MKI67) %>% 
        mutate(sample_clean = gsub("-.*","",sampleID))
    
    plts = make_plots(psi, index)

    # make figdata
    # figdata = make_figdata()
    
    # save
    save_plots(plts, figs_dir)
    # save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
