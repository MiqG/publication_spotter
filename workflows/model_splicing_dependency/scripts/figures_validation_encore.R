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
require(tidytext)

ROOT = here::here()
source(file.path(ROOT,"src","R","utils.R"))

# vaiables 
THRESH_DPSI = 10

# Development
# -----------
# RAW_DIR = file.path(ROOT,"data","raw")
# PREP_DIR = file.path(ROOT,"data","prep")
# RESULTS_DIR = file.path(ROOT,"results","model_splicing_dependency")
# rnai_file = file.path(PREP_DIR,"demeter2","CCLE.tsv.gz")
# delta_psi_file = file.path(RESULTS_DIR,'files','ENCORE','delta_psi-EX.tsv.gz')
# norm_spldep_file = file.path(RESULTS_DIR,"files","ENCORE","norm_spldep-EX.tsv.gz")
# selected_events_file = file.path(RESULTS_DIR,"files","selected_models-EX.txt")
# event_info_file = file.path(RAW_DIR,"VastDB","EVENT_INFO-hg38_noseqs.tsv")
# metadata_file = file.path(PREP_DIR,'metadata','ENCORE.tsv.gz')

##### FUNCTIONS #####
plot_encore_validation = function(){
    plts = list()
    
    cells_oi = metadata %>% pull(DepMap_ID) %>% unique()
    genes_kd = metadata %>% pull(KD) %>% unique() %>% na.omit()
    event_annot = event_info %>% distinct(GENE, EVENT)
    
    genedep = rnai %>% 
        filter(index%in%genes_kd) %>% 
        dplyr::select(c("index",cells_oi)) %>% 
        pivot_longer(-index, names_to="DepMap_ID", values_to="demeter2") %>% 
        rename(KD=index) %>%
        drop_na()
    
    dpsi = delta_psi %>%
        filter(EVENT%in%selected_events) %>%
        pivot_longer(-EVENT, names_to="sampleID", values_to="deltaPSI") %>%
        filter(abs(deltaPSI)>THRESH_DPSI) %>%
        rename(index=EVENT) %>%
        drop_na()
    
    X = norm_spldep %>% 
        filter(index%in%selected_events) %>% 
        pivot_longer(-index, names_to="sampleID", values_to="nspldep") %>%
        filter(is.finite(nspldep)) %>%
        left_join(dpsi, by=c("index","sampleID")) %>%
        left_join(
            metadata %>% mutate(lab=paste0(KD,"_",replicate)),
            by="sampleID"
        ) %>%
        left_join(genedep, by=c("KD","DepMap_ID")) %>% 
        mutate(sign_nspldep = sign(deltaPSI) * nspldep) %>% 
        drop_na(deltaPSI)
    
    # the most harmful change
    X %>% group_by(cell_line, KD) %>% slice_max(nspldep, n=1) %>% ggscatter(x="sign_nspldep",y="demeter2",color="KD",palette=get_palette("simpsons",length(genes_kd))) + guides(color="none") + stat_cor(method="spearman") +  geom_smooth(method="lm", linetype="dashed", color="black") + facet_wrap(~cell_line, ncol=1)
    
    # a sum of changes
    X %>% group_by(cell_line, KD, demeter2) %>% filter(nspldep>0) %>% summarize(pred=sum(sign_nspldep))  %>% ggscatter(x="pred",y="demeter2",color="KD",palette=get_palette("simpsons",length(genes_kd))) + guides(color="none") + stat_cor(method="spearman") +  geom_smooth(method="lm", linetype="dashed", color="black") + facet_wrap(~cell_line, ncol=1, scales="free")
    
    # correlation between cell lines gene dependencies
    genedep %>% left_join(metadata %>% distinct(DepMap_ID,cell_line)) %>% pivot_wider(id_cols="KD", names_from="cell_line", values_from="demeter2") %>% ggscatter(x="K562", y="HepG2") + stat_cor(method="spearman") + geom_smooth(method="lm", linetype="dashed", color="black")
    
    return(plts)
}


make_plots = function(){
    plts = list(
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(models, rnai_stats, cancer_events, 
                        eval_pvalue, eval_corr, 
                        enrichment, indices, spldep_stats, 
                        gene_mut_freq, event_mut_freq,
                        rnai, spldep, splicing, genexpr){
    # prep enrichments
    df_enrichs = do.call(rbind,
        lapply(names(enrichment), function(e){
            res = as.data.frame(enrichment[[e]])
            res[["ontology"]] = e
            return(res)
        })
    )
    
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
    # model selection
    ## controls
    save_plt(plts, "model_sel-deps_sorted_vs_std_ctl_neg", ".pdf", figs_dir, width=5, height=5)
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
    metadata = read_tsv(metadata_file)
    delta_psi = read_tsv(delta_psi_file)
    rnai = read_tsv(rnai_file)
    norm_spldep = read_tsv(norm_spldep_file)
    selected_events = readLines(selected_events_file)
    event_info = read_tsv(event_info_file)
    
    plts = make_plots()

    # make figdata
    figdata = make_figdata()
    
    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}