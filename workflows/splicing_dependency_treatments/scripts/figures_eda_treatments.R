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
# EDA associations
# - pvalues
# - n drug/events significantly associated
# - association agreement between GDSC1 and GDSC2
# - overlaps with cell lines used to model splicing dependency
# - overlaps with cell lines used in GDSC1 and GDSC2
# 
# Rankings
# - significant associations with drug target(s)
# - significant associations distance from drug target(s)
# - significant associations better than with drug target(s)
# - drugs that share target(s) ? 
# - significant associations with drug whose target is not significantly associated
#   at similar overall ranking levels
#
# Association profiles
# - splicing dependency coefficients separate drugs in clusters that separate 
#   their targets (splicing dependencies that differentiate the clusters?)

require(argparse)
require(tidyverse)
require(ggpubr)
require(scattermore)

ROOT = here::here()

# variables
THRESH_FDR = 0.1
THRESH_PVALUE = 0.05
THRESH_NOBS = 20
RANDOM_SEED = 1234

# formatting
PAL_SINGLE_ACCENT = "orange"
PAL_SINGLE_LIGHT = "#6AC2BF"
PAL_SINGLE_DARK = "#716454"
LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"


# Development
# -----------
# RAW_DIR = file.path(ROOT,"data","raw")
# PREP_DIR = file.path(ROOT,"data","prep")
# RESULTS_DIR = file.path(ROOT,'results','splicing_dependency_treatments')
# drug_targets_file = file.path(PREP_DIR,'drug_screens','drug_targets.tsv.gz')
# protein_impact_file = file.path(ROOT,"data","raw","VastDB","PROT_IMPACT-hg38-v3.tab.gz")
# event_info_file = file.path(RAW_DIR,"VastDB","EVENT_INFO-hg38_noseqs.tsv")
# gene_info_file = file.path(RAW_DIR,"ENSEMBL","gene_annotation-hg38.tsv.gz")
# metadata_file = file.path(PREP_DIR,"metadata","Zhang2022.tsv.gz")
# spldep_file = file.path(RESULTS_DIR,'files','Zhang2022','splicing_dependency-EX','mean.tsv.gz')
# splicing_file = file.path(PREP_DIR,"event_psi","Zhang2022-EX.tsv.gz")
# estimated_response_file = file.path(RESULTS_DIR,"files","Zhang2022","estimated_drug_response_by_drug-EX.tsv.gz")

# figs_dir = file.path(RESULTS_DIR,"figures","eda_treatments")

##### FUNCTIONS #####
plot_eda_metadata = function(metadata){
    X = metadata 
    X[["sample_type_num"]] = as.numeric(X[["sample_type"]]) # cannot do it with dplyr, weird
    X = X %>% 
        group_by(patientID, treatment) %>%
        mutate(sort_var = sum(sample_type_num)) %>%
        ungroup() %>%
        group_by(patientID) %>%
        mutate(sort_var = max(sort_var)) %>%
        ungroup() %>%
        count(patientID, treatment, sample_type, sample_location, sort_var) %>%
        mutate(patientID = fct_reorder(patientID, sort_var))
        
    plts = list()
    
    # 
    plts[["eda_metadata-sample_info-balloon"]] = X %>% 
        ggplot(aes(x=patientID, y=sample_type, group=patientID)) + 
        geom_point(aes(color=treatment, size=n)) + 
        geom_line(aes(color=treatment)) + 
        facet_grid(treatment~sample_location, scales="free_y") +
        labs(x="Patient", y="Sample Type", color="Treatment") + 
        coord_flip() + 
        theme_pubr(x.text.angle = 45)
    
    return(plts)
}

make_plots = function(){
    plts = list(
        plot_eda_associations(models, drug_screen),
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(){
    
    figdata = list(
        "drug_event_assoc" = list(
            "model_summaries" = models,
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
                    font.tickslab=6, font.family=FONT_FAMILY)    
    }
    filename = file.path(directory,paste0(plt_name,extension))
    save_plot(filename, plt, base_width=width, base_height=height, dpi=dpi, units="cm")
}


save_plots = function(plts, figs_dir){
    # drug-event associations
    save_plt(plts, "associations-lr_pvalues", ".pdf", figs_dir, width=5, height=5)
    
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
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    metadata = read_tsv(metadata_file)
    spldep = read_tsv(spldep_file)
    splicing = read_tsv(splicing_file)
    estimated_response = read_tsv(estimated_response_file)
    
    drug_targets = read_tsv(drug_targets_file)
    protein_impact = read_tsv(protein_impact_file)
    event_info = read_tsv(event_info_file)
    gene_info = read_tsv(gene_info_file)
    
    # prep inputs
    metadata = metadata %>%
        dplyr::rename(
            sampleID = RNASeq_ID,
            patientID = `EOC ID`,
            treatment = `Treatment strategy`,
            treatment_outcome = `Primary therapy outcome`,
            sample_type = `Stage of sampling`,
            sample_location = Sample,
            relapse_type = `Relapse type`
            
        ) %>%
        mutate(
            sample_type = factor(sample_type, levels = c("Pre-chemo","Post-chemo","Relapse"))
        ) %>%
        separate_rows(sampleID, sep="; ") %>%
        drop_na(sampleID) %>%
        distinct()
    
    # make plots
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