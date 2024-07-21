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
require(scattermore)

# variables
RANDOM_SEED = 1234

# formatting
PAL_SINGLE_ACCENT = "orange"
PAL_DUAL = "#0C0A3E"
LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,"data","raw")
# PREP_DIR = file.path(ROOT,"data","prep")
# RESULTS_DIR = file.path(ROOT,"results","model_splicing_dependency")

# selected_events_file = file.path(RESULTS_DIR,"files","selected_models-EX.txt")
# psi_file = file.path(PREP_DIR,"event_psi","CCLE-EX.tsv.gz")
# total_reads_file = file.path(PREP_DIR,"event_total_reads","CCLE-EX.tsv.gz")
# spldep_file = file.path(RESULTS_DIR,"files","splicing_dependency-EX","mean.tsv.gz")
# psi_simulated_file = file.path(RESULTS_DIR,"files","psi_uncertainty","psi_simulated","CCLE-EX.tsv.gz")
# spldep_simulated_file = file.path(RESULTS_DIR,"files","psi_uncertainty",'splicing_dependency-EX','CCLE','mean.tsv.gz')
# figs_dir = file.path(RESULTS_DIR,'figures','psi_uncertainty')

##### FUNCTIONS #####
simulate_psi_vs_read_count = function(n_reads_size, n_simulation=1000){
    set.seed(RANDOM_SEED)
    
    psis = seq(0,1,0.05)
    simulations = lapply(psis, function(psi){
        simulation = rbinom(n=n_simulation, size=n_reads_size, prob=psi) / n_reads_size
        simulation = simulation * 100 # rescale
        simulation_summary = data.frame(
            psi_prob = psi*100, # rescale
            n_reads = n_reads_size,
            simulation_mean = mean(simulation),
            simulation_median = median(simulation),
            simulation_std = sd(simulation),
            simulation_q25 = quantile(simulation, probs=0.25),
            simulation_q75 = quantile(simulation, probs=0.75)
        ) %>%
        mutate(
            simulation_iqr = simulation_q75 - simulation_q25,
            simulation_error = mean(psi*100 - simulation)
        )
    }) %>% bind_rows()
    
    return(simulations)
}


plot_simulations = function(simulations){
    plts = list()
    
    X = simulations
    
    plts[["simulations-psi_vs_std_vs_n_reads-line"]] = X %>%
        mutate(n_reads = as.factor(n_reads), simulation_error=abs(simulation_error)) %>%
        ggline(x="psi_prob", y="simulation_error", color="n_reads", numeric.x.axis=TRUE, point.size=0.25) +
        labs(x="PSI", y="|PSI Real - PSI with Simulated Error|", color="N Reads to Compute PSI")

    return(plts)
}


plot_uncertainty_quant = function(uncertainty_quant){
    plts = list()
    
    X = uncertainty_quant
    
    plts[["uncertainty_quant-distr_total_reads-violin"]] = X %>%
        mutate(total_reads = log10(total_reads+1)) %>%
        ggviolin(x="sampleID", y="total_reads", fill=PAL_SINGLE_ACCENT, color=NA, trim=TRUE) + 
        geom_boxplot(fill=NA, outlier.size=0.1, width=0.2) +
        labs(x="Sample", y="log10(Total Exon Reads+1)") +
        theme_pubr(x.text.angle=70)
    
    plts[["uncertainty_quant-psi_real_vs_simulated-scatter"]] = X %>%
        ggplot(aes(x=psi_real, y=psi_simulated)) +
        geom_scattermore(pixels = c(1000,1000), pointsize=15, alpha=0.5, color="black") +
        stat_cor(method="pearson", size=FONT_SIZE, family=FONT_FAMILY) +
        theme_pubr() +
        facet_wrap(~sampleID) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="PSI Real", y="PSI with Simulated Error")
    
    plts[["uncertainty_quant-total_reads_vs_diff_psi-scatter"]] = X %>%
        mutate(total_reads = log10(total_reads+1)) %>%
        ggplot(aes(x=total_reads, y=abs(diff_psi))) +
        geom_scattermore(pixels = c(1000,1000), pointsize=15, alpha=0.5, color="black") +
        stat_cor(method="pearson", size=FONT_SIZE, family=FONT_FAMILY) +
        geom_vline(xintercept=log10(10+1), size=LINE_SIZE, color="black", linetype="dashed") +
        theme_pubr() +
        facet_wrap(~sampleID) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="log10(Total Exon Reads+1)", y="|PSI Real - PSI with Simulated Error|")
    
    
    plts[["uncertainty_quant-spldep_real_vs_simulated-scatter"]] = X %>%
        ggplot(aes(x=spldep_real, y=spldep_simulated)) +
        geom_scattermore(pixels = c(1000,1000), pointsize=15, alpha=0.5, color="black") +
        stat_cor(method="pearson", size=FONT_SIZE, family=FONT_FAMILY) +
        theme_pubr() +
        facet_wrap(~sampleID) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Spl. Dep. Real", y="Spl. Dep. with Simulated Error")
    
    
    plts[["uncertainty_quant-total_reads_vs_diff_spldep-scatter"]] = X %>%
        mutate(total_reads = log10(total_reads+1)) %>%
        ggplot(aes(x=total_reads, y=abs(diff_spldep))) +
        geom_scattermore(pixels = c(1000,1000), pointsize=15, alpha=0.5, color="black") +
        stat_cor(method="pearson", size=FONT_SIZE, family=FONT_FAMILY) +
        geom_vline(xintercept=log10(10+1), size=LINE_SIZE, color="black", linetype="dashed") +
        theme_pubr() +
        facet_wrap(~sampleID) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="log10(Total Exon Reads+1)", y="|Spl. Dep. Real - Spl. Dep. with Simulated Error|")
    
    
    plts[["uncertainty_quant-diff_psi_vs_diff_spldep-scatter"]] = X %>%
        ggplot(aes(x=abs(diff_psi), y=abs(diff_spldep))) +
        geom_scattermore(pixels = c(1000,1000), pointsize=15, alpha=0.5, color="black") +
        stat_cor(method="pearson", size=FONT_SIZE, family=FONT_FAMILY) +
        theme_pubr() +
        facet_wrap(~sampleID) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="|PSI Real - PSI with Simulated Error|", y="|Spl. Dep. Real - Spl. Dep. with Simulated Error|")
    
    
    return(plts)
}


make_plots = function(simulations, uncertainty_quant){
    plts = list(
        plot_simulations(simulations),
        plot_uncertainty_quant(uncertainty_quant)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(simulations, uncertainty_quant){
    figdata = list(
        "psi_uncertainty" = list(
            "simulations" = simulations,
            "uncertainty_quant" = uncertainty_quant
        )
    )
    return(figdata)
}

make_source_data = function(plts){
    
    source_data = list(
        # SUPPLEMENTARY FIGURE 2
        ## Sup. Fig. 2a
        "supfig02a" = plts[["simulations-psi_vs_std_vs_n_reads-line"]][["data"]],
        ## Sup. Fig. 2b
        "supfig02b" = plts[["uncertainty_quant-distr_total_reads-violin"]][["data"]],
        ## Sup. Fig. 2c
        "supfig02c" = plts[["uncertainty_quant-total_reads_vs_diff_spldep-scatter"]][["data"]] %>% dplyr::select(-c(psi_real,psi_simulated,spldep_simulated,spldep_real)),
        ## Sup. Fig. 2d
        "supfig02d" = plts[["uncertainty_quant-psi_real_vs_simulated-scatter"]][["data"]] %>% dplyr::select(-c(total_reads,diff_psi,diff_spldep))
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
    save_plt(plts, "simulations-psi_vs_std_vs_n_reads-line", ".pdf", figs_dir, width=5, height=5)
    
    save_plt(plts, "uncertainty_quant-distr_total_reads-violin", ".pdf", figs_dir, width=3, height=5.5)
    save_plt(plts, "uncertainty_quant-psi_real_vs_simulated-scatter", ".pdf", figs_dir, width=5, height=5.5)
    save_plt(plts, "uncertainty_quant-total_reads_vs_diff_psi-scatter", ".pdf", figs_dir, width=5, height=5.5)
    save_plt(plts, "uncertainty_quant-spldep_real_vs_simulated-scatter", ".pdf", figs_dir, width=5, height=5.5)
    save_plt(plts, "uncertainty_quant-total_reads_vs_diff_spldep-scatter", ".pdf", figs_dir, width=5, height=5.5)
    save_plt(plts, "uncertainty_quant-diff_psi_vs_diff_spldep-scatter", ".pdf", figs_dir, width=5, height=5.5)
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
        make_option("--psi_file", type="character"),
        make_option("--total_reads_file", type="character"),
        make_option("--psi_simulated_file", type="character"),
        make_option("--spldep_file", type="character"),
        make_option("--spldep_simulated_file", type="character"),
        make_option("--selected_events_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    psi_file = args[["psi_file"]]
    total_reads_file = args[["total_reads_file"]]
    psi_simulated_file = args[["psi_simulated_file"]]
    spldep_file = args[["spldep_file"]]
    spldep_simulated_file = args[["spldep_simulated_file"]]
    selected_events_file = args[["selected_events_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    selected_events = readLines(selected_events_file)
    ## samples with lowest and highest read counts
    psi = read_tsv(psi_file, col_select=c("EVENT","ACH-000934","ACH-000143"))
    total_reads = read_tsv(total_reads_file, col_select=c("EVENT","ACH-000934","ACH-000143"))
    spldep = read_tsv(spldep_file, col_select=c("index","ACH-000934","ACH-000143"))
    psi_simulated = read_tsv(psi_simulated_file, col_select=c("EVENT","ACH-000934","ACH-000143"))
    spldep_simulated = read_tsv(spldep_simulated_file, col_select=c("index","ACH-000934","ACH-000143"))
    
    # preprocess
    uncertainty_quant = psi %>%
        filter(EVENT %in% selected_events) %>%
        pivot_longer(-EVENT, names_to="sampleID", values_to="psi_real") %>%
        left_join(
            total_reads %>%
            filter(EVENT %in% selected_events) %>%
            pivot_longer(-EVENT, names_to="sampleID", values_to="total_reads"),
            by=c("EVENT","sampleID")
        ) %>%
        left_join(
            spldep %>%
            dplyr::rename(EVENT=index) %>%
            filter(EVENT %in% selected_events) %>%
            pivot_longer(-EVENT, names_to="sampleID", values_to="spldep_real"),
            by=c("EVENT","sampleID")
        ) %>%
        left_join(
            psi_simulated %>%
            filter(EVENT %in% selected_events) %>%
            pivot_longer(-EVENT, names_to="sampleID", values_to="psi_simulated"),
            by=c("EVENT","sampleID")
        ) %>%
        left_join(
            spldep_simulated %>%
            dplyr::rename(EVENT=index) %>%
            filter(EVENT %in% selected_events) %>%
            pivot_longer(-EVENT, names_to="sampleID", values_to="spldep_simulated"),
            by=c("EVENT","sampleID")
        ) %>%
        mutate(
            diff_psi = psi_real - psi_simulated,
            diff_spldep = spldep_real - spldep_simulated,
            sampleID = factor(sampleID, levels=c("ACH-000934","ACH-000143"))
        )
    
    # simulate
    n_reads = c(2,10,100,1000)
    simulations = lapply(n_reads, simulate_psi_vs_read_count) %>% bind_rows()
    
    # plot
    plts = make_plots(simulations, uncertainty_quant)

    # make figdata
    figdata = make_figdata(simulations, uncertainty_quant)
    
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
