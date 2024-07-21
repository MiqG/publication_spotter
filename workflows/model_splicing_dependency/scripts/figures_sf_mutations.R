#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#

require(optparse)
require(tidyverse)
require(ggpubr)
require(cowplot)
require(extrafont)
require(tidytext)

LINE_SIZE = 0.25
FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

THRESH_PVAL = 0.05
SFS_OI = c(
    "RBM10",
    "SF3B1",
    "SRSF1",
    "SRSF2",
    "ZRSR2",
    "U2AF1"
)

PAL_DUAL = c("#C4D6B0","#477998")

# Development
# -----------
# ROOT = here::here()
# PREP_DIR = file.path(ROOT,"data","prep")
# RAW_DIR = file.path(ROOT,"data","raw")
# RESULTS_DIR = file.path(ROOT,"results","model_splicing_dependency")
# snv_file = file.path(RAW_DIR,'DepMap','achilles_ccle','CCLE_mutations.csv')
# max_harm_file = file.path(RESULTS_DIR,"files","splicing_dependency-EX","max_harm_score-mean.tsv.gz")
# selected_exons_file = file.path(RESULTS_DIR,'files','selected_models-EX.txt')
# event_info_file = file.path(RAW_DIR,"VastDB","EVENT_INFO-hg38_noseqs.tsv")
# cancer_events_file = file.path(ROOT,"support","cancer_events.tsv")
# ascanceratlas_file = file.path(RAW_DIR,"ASCancerAtlas","CASE_all-VastDB_mapped.tsv.gz")
# figs_dir = file.path(RESULTS_DIR,"figures","sf_mutations")

##### FUNCTIONS #####
compute_diff_max_harm = function(max_harm, mutations){
    X = max_harm %>%
        pivot_longer(-event_gene, names_to="DepMap_ID", values_to="max_harm_score") %>%
        drop_na(max_harm_score)
    
    mutations_oi = mutations %>%
        count(Hugo_Symbol, Variant_Classification) %>%
        arrange(n) %>%
        filter(n>=10)
    
    results = lapply(1:nrow(mutations_oi), function(row_i){
        # add mutation info
        sf_oi = mutations_oi[row_i,] %>% pull(Hugo_Symbol)
        variant_class = mutations_oi[row_i,] %>% pull(Variant_Classification)
        mutated_samples = mutations %>%
            filter(Variant_Classification==variant_class & Hugo_Symbol==sf_oi) %>%
            pull(DepMap_ID)
        x = X %>% mutate(is_mutated = DepMap_ID %in% mutated_samples)
        
        # do not consider those events without both categories
        events_oi = x %>%
            distinct(event_gene, is_mutated) %>%
            count(event_gene) %>%
            filter(n==2) %>%
            pull(event_gene)
        x = x %>% filter(event_gene %in% events_oi)
        
        # test
        result = compare_means(
            max_harm_score ~ is_mutated, x, method="wilcox.test", 
            group.by="event_gene", p.adjust.method="fdr"
        )
        
        # add more info
        result = result %>%
            left_join(
                x %>%
                group_by(event_gene, is_mutated) %>%
                summarize(
                    max_harm_median = median(max_harm_score),
                    max_harm_q25 = quantile(max_harm_score, 0.25),
                    max_harm_q75 = quantile(max_harm_score, 0.75),
                    n_obs = n()
                ) %>%
                ungroup(),
                by = "event_gene"
            ) %>%
            mutate(
                mutated_sf = sf_oi,
                mutation_effect = variant_class
            )
        return(result)
    }) %>%
    do.call(rbind, .)
    
    return(results)
}


plot_diff_max_harm = function(diff_max_harm){
    plts = list()
    
    X = diff_max_harm 
    mutated_max_harms = X %>%
        pivot_wider(
            id_cols=c("event_gene","mutated_sf","mutation_effect"), names_from="is_mutated", 
            values_from="max_harm_median", names_prefix = "max_harm_median_"
        ) %>%
        mutate(max_harm_median_diff = abs(max_harm_median_TRUE) - abs(max_harm_median_FALSE))
    X = X %>%
        left_join(mutated_max_harms, by=c("event_gene","mutated_sf","mutation_effect"))
    
    # how many significant differences per sf and mutation effect?
    plts[["diff_max_harm-freqs_vs_sf_vs_mut_effect-bar"]] = X %>%
        filter(is_significant & is_known_driver) %>%
        distinct(mutated_sf, mutation_effect, event_gene) %>%
        count(mutated_sf, mutation_effect) %>%
        mutate(mutated_sf = factor(mutated_sf, levels=SFS_OI)) %>%
        ggbarplot(
            x="mutated_sf", y="n", fill="mutation_effect", 
            color=NA, palette="jco", label=TRUE, 
            lab.size=FONT_SIZE, lab.family=FONT_FAMILY
        ) +
        theme_pubr(x.text.angle=70) +
        labs(x="Mutated SF", y="Count", fill="Mutation Effect")
    
    # which are the top differences per sf and mutation effect?
    plts[["diff_max_harm-top_diffs-bar"]] = X %>% 
        filter(is_significant & is_known_driver & mutation_effect=="Missense_Mutation") %>%
        mutate(
            event_gene = as.factor(event_gene),
            name = reorder_within(
                event_gene, max_harm_median_diff, list(mutated_sf, mutation_effect)
            )
        ) %>%
        ggplot(aes(x=event_gene, y=max_harm_median, group=is_mutated)) +
        geom_col(aes(fill=is_mutated), color=NA, position=position_dodge(0.9)) +
        geom_errorbar(
            aes(ymin=max_harm_q25, ymax=max_harm_q75), position=position_dodge(.9),
            width=0.2
        ) +
        fill_palette(PAL_DUAL) +
        theme_pubr(x.text.angle=70) +
        scale_x_reordered() +
        facet_grid(~mutated_sf, scales="free_x", space="free") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Event & Gene", y="median(Max. Harm Score)", fill="Is Mutated")
    
    return(plts)
}


make_plots = function(diff_max_harm){
    plts = list(
        plot_diff_max_harm(diff_max_harm)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(diff_max_harm){
    
    figdata = list(
        "sf_mutations" = list(
            "diff_max_harm" = diff_max_harm %>% 
                filter(is_significant & is_known_driver)
        )
    )
    return(figdata)
}

make_source_data = function(plts){
    
    source_data = list(
        # SUPPLEMENTARY FIGURE 4
        ## Sup. Fig. 4d
        "supfig04d" = plts[["diff_max_harm-top_diffs-bar"]][["data"]]
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
    save_plt(plts, "diff_max_harm-freqs_vs_sf_vs_mut_effect-bar", ".pdf", figs_dir, width=4, height=6)
    save_plt(plts, "diff_max_harm-top_diffs-bar", ".pdf", figs_dir, width=11, height=9)
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
        make_option("--snv_file", type="character"),
        make_option("--max_harm_file", type="character"),
        make_option("--selected_exons_file", type="character"),
        make_option("--cancer_events_file", type="character"),
        make_option("--ascanceratlas_file", type="character"),
        make_option("--event_info_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    snv_file = args[["snv_file"]]
    max_harm_file = args[["max_harm_file"]]
    selected_exons_file = args[["selected_exons_file"]]
    cancer_events_file = args[["cancer_events_file"]]
    ascanceratlas_file = args[["ascanceratlas_file"]]
    event_info_file = args[["event_info_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    snv = read_csv(snv_file)
    max_harm = read_tsv(max_harm_file)
    selected_exons = readLines(selected_exons_file)
    event_info = read_tsv(event_info_file) %>%
        mutate(event_gene=paste0(EVENT,"_",GENE))
    cancer_events = read_tsv(cancer_events_file)
    ascanceratlas = read_tsv(ascanceratlas_file)

    # prep
    ## prep cancer events
    ascanceratlas = ascanceratlas %>%
        dplyr::rename(EVENT = EVENT_perf) %>%
        left_join(
            event_info %>% distinct(EVENT,GENE),
            by=c("EVENT")
        ) %>%
        drop_na(EVENT) %>%
        mutate(source="ASCancerAtlas")
    
    cols_oi = c('EVENT','ENSEMBL','GENE','length','source_authors',
                'source_year','source_doi','description','event_id','source')
    cancer_events = cancer_events %>%
        mutate(source="HandCurated") %>%
        bind_rows(ascanceratlas) %>%
        dplyr::select(all_of(cols_oi)) %>%
        distinct() %>%
        mutate(event_gene=paste0(EVENT,"_",GENE))

    ## subset max harms
    max_harm = max_harm %>% 
        filter(index %in% selected_exons) %>%
        left_join(event_info[,c("EVENT","event_gene")], by=c("index"="EVENT")) %>%
        dplyr::select(-index) 
    
    ## mutations in splicing factors of interest
    mutations = snv %>% 
        # one mutation per SF and cancer cell line
        distinct(DepMap_ID, Hugo_Symbol, Variant_Classification) %>%
        filter(Variant_Classification != "Silent")
    mutations = mutations %>% filter(Hugo_Symbol %in% SFS_OI)
    
    diff_max_harm = compute_diff_max_harm(max_harm, mutations) 
    
    ## add info
    diff_max_harm = diff_max_harm %>%
        mutate(
            is_significant = p < THRESH_PVAL,
            is_known_driver = event_gene %in% cancer_events[["event_gene"]]
        )
    
    # plot
    plts = make_plots(diff_max_harm)

    # make figdata
    figdata = make_figdata(diff_max_harm)
    
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
