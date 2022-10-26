require(tidyverse)
require(ggpubr)
require(cowplot)
require(ggrepel)
require(extrafont)

# formatting
PAL_SINGLE_ACCENT = "orange"
PAL_SINGLE_LIGHT = "#efb300ff"#"#6AC2BF"
PAL_SINGLE_DARK = "#007e67ff"
PAL_SINGLE_NEUTRAL = "#716454"
PAL_DUAL = c(PAL_SINGLE_DARK, PAL_SINGLE_LIGHT)
PAL_FDR_DARK = "#005AB5"
PAL_FDR_LIGHT = "#DC3220"

LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"


# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,'support')
# RESULTS_DIR = file.path(ROOT,'results','splicing_dependency_validation')

# annotation_file = file.path(RAW_DIR,'VastDB','event_annotation-Hs2.tsv.gz')
# event_info_file = file.path(RAW_DIR,'VastDB','EVENT_INFO-hg38_noseqs.tsv')

# validation_splicing_file = file.path(RAW_DIR,"experiments","validation_therapeutic_potential","20220928-psi-aso.tsv")
# validation_spldep_file = file.path(RESULTS_DIR,"files","splicing_dependency-EX","mean.tsv.gz")
# validation_od_file = file.path(RAW_DIR,"experiments","validation_therapeutic_potential","clonogenic_assay-od-merged.tsv")

# figs_dir = file.path(RESULTS_DIR,'figures','validation')


##### FUNCTIONS #####
plot_validation = function(validation_clonogenic, validation_harm_scores){
    
    PAL_REPLICATES = setNames(get_palette("jco",2), 1:2)
    
    plts = list()
    
    
    plts[["validation-od"]] = validation_clonogenic %>%
        mutate(event_gene = fct_reorder(event_gene, -od, mean)) %>%
        ggplot(aes(x=event_gene, y=od)) +
        geom_boxplot(width=0.5, outlier.shape=NA) +
        geom_jitter(aes(color=as.factor(replicate_biological)), size=1, width=0.1) +
        color_palette(PAL_REPLICATES) +
        labs(x="Condition", y="OD570", color="Replicate") +
        theme_pubr(x.text.angle = 70) +
        stat_compare_means(ref.group="CONTROL_NEG", method="t.test", label="p.format", 
                           size=FONT_SIZE, family=FONT_FAMILY)
    
    od_ctl = validation_clonogenic %>%
        filter(event_gene=="WATER") %>%
        group_by(CCLE_Name, replicate_biological) %>%
        summarize(od_ctl = median(od)) %>%
        ungroup()
    
    plts[["validation-od_vs_harm"]] = validation_clonogenic %>%
        left_join(od_ctl, by=c("CCLE_Name","replicate_biological")) %>%
        mutate(od_fc = log2(od / od_ctl)) %>%
        #mutate(od_fc = od - od_ctl) %>%
        group_by(CCLE_Name, event_gene, replicate_biological) %>%
        summarize(
            od_fc = median(od_fc),
            replicate_biological = as.factor(replicate_biological)
        ) %>%
        ungroup() %>%
        distinct() %>%
        left_join(validation_harm_scores, by=c("CCLE_Name","event_gene")) %>%
        drop_na(harm_score, od_fc) %>%
        ggplot(aes(x=harm_score, y=od_fc)) +
        geom_smooth(method="lm", linetype="dashed", color="black", size=LINE_SIZE) +
        geom_point(aes(color=replicate_biological)) +
        geom_text_repel(aes(label=event_gene), size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1) +
        stat_cor(method="pearson", size=FONT_SIZE, family=FONT_FAMILY) +
        facet_wrap(~CCLE_Name+replicate_biological) +
        color_palette(PAL_REPLICATES) +
        theme_pubr() +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY), aspect.ratio=1) +
        labs(x="Harm Score", y="Obs. Cell Prolif. Effect", color="Replicate")
    
    return(plts)
}


make_plots = function(validation_clonogenic, validation_harm_scores){
    plts = list(
        plot_validation(validation_clonogenic, validation_harm_scores)
    )
    plts = do.call(c,plts)
    return(plts)
}


save_plt = function(plts, plt_name, extension='.pdf', 
                    directory='', dpi=350, format=FALSE,
                    width = par("din")[1], height = par("din")[2]){
    plt = plts[[plt_name]]
    if (format){
        plt = ggpar(plt, font.title=10, font.subtitle=10, font.caption=10, 
                    font.x=8, font.y=8, font.legend=8,
                    font.tickslab=6, font.family='Arial')    
    }
    filename = file.path(directory,paste0(plt_name,extension))
    save_plot(filename, plt, base_width=width, base_height=height, dpi=dpi, units='cm')
}


save_plots = function(plts, figs_dir){
    save_plt(plts, "validation-od", '.pdf', figs_dir, width=18, height=12)
    save_plt(plts, "validation-od_vs_harm", '.pdf', figs_dir, width=12, height=9)
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
    figs_dir = args$figs_dir
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    event_info = read_tsv(event_info_file)
    annot = read_tsv(annotation_file) %>%
        mutate(event_gene = paste0(EVENT,"_",GENE))
    events_genes = annot %>%
        dplyr::select(EVENT,GENE,ENSEMBL,event_gene)
    
    validation_psi = read_tsv(validation_splicing_file) %>%
        left_join(events_genes %>% distinct(EVENT,event_gene), by="EVENT") %>%
        dplyr::select(-EVENT)
    validation_spldep = read_tsv(validation_spldep_file) %>%
        left_join(events_genes %>% distinct(EVENT,event_gene), by=c("index"="EVENT")) %>%
        dplyr::select(-index)
    validation_clonogenic = read_tsv(validation_od_file) %>%
        group_by(CCLE_Name, event_gene, replicate_technical, replicate_biological) %>%
        summarize(od = mean(od)) %>% # summarize OD replicates
        ungroup()
    
    # compute harm scores
    validation_harm_scores = validation_spldep %>%
        pivot_longer(-event_gene, names_to="DepMap_ID", values_to="spldep") %>%
        left_join(validation_psi, by=c("DepMap_ID","event_gene")) %>%
        mutate(
            delta_psi = psi_treated - psi_untreated,
            harm_score = (-1) * delta_psi * spldep
        )
    
    # plot
    plts = make_plots(validation_clonogenic, validation_harm_scores)
    
    # make figdata
    #figdata = make_figdata()

    # save
    save_plots(plts, figs_dir)
    #save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}