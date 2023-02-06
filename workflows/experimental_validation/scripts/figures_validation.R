require(optparse)
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

PAL_REPLICATES = setNames(get_palette("Accent",3), 1:3)

LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,'support')
# RESULTS_DIR = file.path(ROOT,'results','experimental_validation')

# annotation_file = file.path(RAW_DIR,'VastDB','event_annotation-Hs2.tsv.gz')
# event_info_file = file.path(RAW_DIR,'VastDB','EVENT_INFO-hg38_noseqs.tsv')

# splicing_file = file.path(PREP_DIR,"event_psi","inhouse-EX.tsv.gz")
# genexpr_file = file.path(PREP_DIR,"genexpr_tpm","inhouse.tsv.gz")
# spldep_file = file.path(RESULTS_DIR,"files","splicing_dependency-EX","mean.tsv.gz")

# validation_splicing_file = file.path(RAW_DIR,"experiments","validation_therapeutic_potential","20220928-psi-aso.tsv")
# validation_od_file = file.path(RAW_DIR,"experiments","validation_therapeutic_potential","clonogenic_assay-od-merged.tsv")

# figs_dir = file.path(RESULTS_DIR,'figures','validation')


##### FUNCTIONS #####
plot_validation = function(validation_clonogenic, validation_harm_scores){
    
    plts = list()
    
    plts[["validation-od_raw"]] = validation_clonogenic %>%
        mutate(event_gene = fct_reorder(event_gene, -od, mean)) %>%
        ggplot(aes(x=event_gene, y=od)) +
        geom_boxplot(width=0.5, outlier.shape=NA) +
        geom_jitter(aes(color=as.factor(replicate_biological)), size=0.5, width=0.1) +
        color_palette(PAL_REPLICATES) +
        labs(x="Condition", y="Cell Viability (OD570)", color="Replicate") +
        theme_pubr(x.text.angle = 70) +
        stat_compare_means(ref.group="CONTROL_NEG", method="t.test", label="p.signif", 
                           size=FONT_SIZE, family=FONT_FAMILY) + 
        facet_wrap(~CCLE_Name, ncol=1) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    plts[["validation-od_norm"]] = validation_clonogenic %>%
        mutate(event_gene = fct_reorder(event_gene, -od_norm, mean)) %>%
        ggplot(aes(x=event_gene, y=od_norm)) +
        geom_boxplot(width=0.5, outlier.shape=NA) +
        geom_jitter(aes(color=as.factor(replicate_biological)), size=0.5, width=0.1) +
        color_palette(PAL_REPLICATES) +
        labs(x="Condition", y="Cell Viability (Norm. OD570)", color="Replicate") +
        theme_pubr(x.text.angle = 70) +
        stat_compare_means(ref.group="CONTROL_NEG", method="t.test", label="p.signif", 
                           size=FONT_SIZE, family=FONT_FAMILY) + 
        facet_wrap(~CCLE_Name, ncol=1) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    plts[["validation-od_norm_averaged"]] = validation_clonogenic %>%
        distinct(CCLE_Name, replicate_biological, event_gene, od_norm_avg) %>%
        mutate(event_gene = fct_reorder(event_gene, -od_norm_avg, mean)) %>%
        ggplot(aes(x=event_gene, y=od_norm_avg)) +
        geom_boxplot(width=0.5, outlier.shape=NA) +
        geom_jitter(aes(color=as.factor(replicate_biological)), size=0.5, width=0.1) +
        color_palette(PAL_REPLICATES) +
        labs(x="Condition", y="Cell Viability (Norm. OD570)", color="Replicate") +
        theme_pubr(x.text.angle = 70) +
        stat_compare_means(ref.group="CONTROL_NEG", method="t.test", label="p.signif", 
                           size=FONT_SIZE, family=FONT_FAMILY) + 
        facet_wrap(~CCLE_Name, ncol=1) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    

    od_ctl = validation_clonogenic %>%
        filter(event_gene=="WATER") %>%
        group_by(CCLE_Name, replicate_biological) %>%
        summarize(od_ctl = mean(od)) %>%
        ungroup()
    
    x = validation_clonogenic %>%
        left_join(od_ctl, by=c("CCLE_Name","replicate_biological")) %>%
        mutate(od_fc = log2(od / od_ctl)) %>%
        group_by(CCLE_Name, event_gene, replicate_biological) %>%
        summarize(
            od_fc = median(od_fc),
            replicate_biological = as.factor(replicate_biological)
        ) %>%
        ungroup() %>%
        distinct() %>%
        left_join(validation_harm_scores, by=c("CCLE_Name","event_gene")) %>%
        drop_na(harm_score, od_fc)
    
    plts[["validation-od_vs_harm"]] = x %>%
        ggplot(aes(x=harm_score, y=od_fc)) +
        geom_smooth(method="lm", linetype="dashed", color="black", size=LINE_SIZE) +
        geom_point(aes(color=event_gene)) +
        #geom_text_repel(aes(label=event_gene), size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1) +
        stat_cor(method="pearson", size=FONT_SIZE, family=FONT_FAMILY) +
        facet_wrap(~CCLE_Name+replicate_biological, ncol=2, scales="free_y") +
        color_palette(get_palette("simpsons",8)) +
        theme_pubr() +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY), aspect.ratio=1) +
        labs(x="Harm Score", y="Obs. Cell Prolif. Effect", color="Replicate")
    
    plts[["validation-od_vs_harm_combined"]] = x %>%
        # average over biological replicates
        # group_by(CCLE_Name, event_gene) %>%
        # summarize(
        #     harm_score = mean(harm_score, na.rm=TRUE),
        #     od_fc = mean(od_fc, na.rm=TRUE)
        # ) %>%
        # ungroup() %>%
        ggplot(aes(x=harm_score, y=od_fc)) +
        geom_smooth(aes(color=replicate_biological, fill=replicate_biological), method="lm", linetype="dashed", size=LINE_SIZE) +
        geom_point(aes(color=replicate_biological)) +
        stat_cor(aes(color=replicate_biological), method="pearson", size=FONT_SIZE, family=FONT_FAMILY) +
        facet_wrap(~CCLE_Name, ncol=2, scales="free_y") +
        color_palette(PAL_REPLICATES) +
        fill_palette(PAL_REPLICATES) +
        theme_pubr() +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY), aspect.ratio=1) +
        labs(x="Harm Score", y="Obs. Cell Prolif. Effect", color="Event & Gene")
    
    return(plts)
}


make_plots = function(validation_clonogenic, validation_harm_scores){
    plts = list(
        plot_validation(validation_clonogenic, validation_harm_scores)
    )
    plts = do.call(c,plts)
    return(plts)
}

make_figdata = function(validation_clonogenic, validation_harm_scores){

    figdata = list(
        "experiments" = list(
            "validation_clonogenic" = validation_clonogenic,
            "validation_harm_scores" = validation_harm_scores
        )
    )
}

save_plt = function(plts, plt_name, extension='.pdf', 
                    directory='', dpi=350, format=TRUE,
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
    save_plt(plts, "validation-od_raw", '.pdf', figs_dir, width=5, height=10)
    save_plt(plts, "validation-od_norm", '.pdf', figs_dir, width=5, height=10)
    save_plt(plts, "validation-od_norm_averaged", '.pdf', figs_dir, width=5, height=10)
    save_plt(plts, "validation-od_vs_harm", '.pdf', figs_dir, width=12, height=15)
    save_plt(plts, "validation-od_vs_harm_combined", '.pdf', figs_dir, width=8, height=6)

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
        make_option("--event_info_file", type="character"),
        make_option("--validation_splicing_file", type="character"),
        make_option("--splicing_ccle_file", type="character"),
        make_option("--genexpr_ccle_file", type="character"),
        make_option("--spldep_ccle_file", type="character"),
        make_option("--validation_spldep_file", type="character"),
        make_option("--validation_od_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    annotation_file = args[["annotation_file"]]
    event_info_file = args[["event_info_file"]]
    validation_splicing_file = args[["validation_splicing_file"]]
    splicing_ccle_file = args[["splicing_ccle_file"]]
    genexpr_ccle_file = args[["genexpr_ccle_file"]]
    spldep_ccle_file = args[["spldep_ccle_file"]]
    validation_spldep_file = args[["validation_spldep_file"]]
    validation_od_file = args[["validation_od_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    event_info = read_tsv(event_info_file)
    annot = read_tsv(annotation_file) %>%
        mutate(event_gene = paste0(EVENT,"_",GENE))
    events_genes = annot %>%
        dplyr::select(EVENT,GENE,ENSEMBL,event_gene)
    
    splicing = read_tsv(splicing_file)
    genexpr = read_tsv(genexpr_file)
    spldep = read_tsv(spldep_file) %>% dplyr::rename(EVENT=index)
    
    validation_psi = read_tsv(validation_splicing_file) %>%
        left_join(events_genes %>% distinct(EVENT,event_gene), by="EVENT")
    validation_clonogenic = read_tsv(validation_od_file) %>%
        group_by(DepMap_ID, CCLE_Name, event_gene, replicate_technical, replicate_biological) %>%
        summarize(od = mean(od)) %>% # summarize OD replicates
        ungroup() %>%
        separate(event_gene, c("EVENT","GENE"), remove=FALSE)
    
    # compute harm scores
    ## using experimentally measured delta PSIs
    validation_harm_scores = spldep %>%
        pivot_longer(-EVENT, names_to="DepMap_ID", values_to="spldep") %>%
        left_join(validation_psi, by=c("DepMap_ID","EVENT")) %>% 
        # compute harm score
        mutate(
            delta_psi = psi_treated - psi_untreated,
            harm_score = (-1) * delta_psi * spldep
        )
    
    # add missing PSI and Spldep from inhouse cell lines
    missing = validation_clonogenic %>%
        distinct(DepMap_ID, event_gene, EVENT, GENE) %>%
        left_join(validation_harm_scores %>% distinct(CCLE_Name, DepMap_ID), by="DepMap_ID") %>%
        filter(!(DepMap_ID %in% validation_psi[["DepMap_ID"]]))
    
    missing = missing %>%
        # add initial PSI from inhouse RNA seq
        left_join(
            splicing %>% 
                dplyr::select(c(EVENT,unique(missing[["DepMap_ID"]]))) %>%
                pivot_longer(-EVENT, names_to="DepMap_ID", values_to="psi_untreated"),
            by=c("EVENT","DepMap_ID")
        ) %>%
        # add initial gene expression from inhouse RNA seq
        left_join(
            genexpr %>% 
                dplyr::select(c(ID,unique(missing[["DepMap_ID"]]))) %>%
                left_join(events_genes%>% distinct(ENSEMBL,GENE), by=c("ID"="ENSEMBL")) %>%
                pivot_longer(-c("ID","GENE"), names_to="DepMap_ID", values_to="tpm_untreated"),
            by=c("GENE","DepMap_ID")
        ) %>%
        # assume maximum exclusion
        mutate(psi_treated = 0) %>%
        # get splicing dependency from inhouse RNA seq
        left_join(
            spldep %>%
                dplyr::select(c(EVENT,unique(missing[["DepMap_ID"]]))) %>%
                pivot_longer(-EVENT, names_to="DepMap_ID", values_to="spldep"),
            by=c("EVENT","DepMap_ID")
        ) %>%
        # compute harm score
        mutate(
            delta_psi = psi_treated - psi_untreated,
            harm_score = (-1) * delta_psi * spldep
        ) %>%
        drop_na(harm_score) %>%
        # add CCLE Name
        dplyr::select(-CCLE_Name) %>%
        left_join(
            validation_clonogenic %>% distinct(CCLE_Name, DepMap_ID),
            by = "DepMap_ID"
        )
    
    validation_harm_scores = validation_harm_scores %>%
        bind_rows(missing)
    
    # normalize OD by biological replicate
    od_ctl_neg = validation_clonogenic %>%
        filter(event_gene=="CONTROL_NEG") %>%
        group_by(CCLE_Name, replicate_biological) %>%
        summarize(od_ctl_neg = mean(od)) %>%
        ungroup()
    
    validation_clonogenic = validation_clonogenic %>%
        left_join(od_ctl_neg, by=c("CCLE_Name","replicate_biological")) %>%
        mutate(
            od_norm = od / od_ctl_neg,
            replicate_biological = factor(replicate_biological, levels=1:3)
        ) %>%
        group_by(CCLE_Name, replicate_biological, event_gene) %>%
        mutate(od_norm_avg = mean(od_norm)) %>%
        ungroup()
    
    # drop KRAS because we do not know if the SSO worked
    validation_clonogenic = validation_clonogenic %>% filter(event_gene != "HsaEX0034998_KRAS")
    validation_harm_scores = validation_harm_scores %>% filter(event_gene != "HsaEX0034998_KRAS")
    
    # plot
    plts = make_plots(validation_clonogenic, validation_harm_scores)
    
    # make figdata
    figdata = make_figdata(validation_clonogenic, validation_harm_scores)

    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}