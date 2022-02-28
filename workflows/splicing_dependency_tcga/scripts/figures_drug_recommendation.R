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
# Drug Recommendations
# - combine partial predictions of drug IC50 into one
# - would we recommend drugs given to patients? By cancer type and subtype.
# - for those cases in which we would give another drug, which would it be?

require(tidyverse)
require(ggpubr)
require(cowplot)
require(extrafont)
require(writexl)

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

# Development
# ----------- 
# PREP_DIR = file.path(ROOT,'data','prep')
# RAW_DIR = file.path(ROOT,'data','raw')
# RESULTS_DIR = file.path(ROOT,'results','splicing_dependency_tcga')
# MODELS_DIR = file.path(ROOT,"results","splicing_dependency_drugs")
# drug_targets_file = file.path(PREP_DIR,'drug_screens','drug_targets.tsv.gz')
# metadata_file = file.path(PREP_DIR,'metadata','PANCAN.tsv.gz')
# drug_treatments_file = file.path(PREP_DIR,'drug_treatments','PANCAN.tsv.gz')
# figs_dir = file.path(RESULTS_DIR,'figures','drug_recommendation')
# drug_response_file = file.path(RESULTS_DIR,'files','PANCAN','estimated_drug_response_by_drug-EX.tsv.gz')
# metadata_subtypes_file = file.path(PREP_DIR,'metadata','PANCAN_subtypes.tsv.gz')
# selected_drugs_file = file.path(MODELS_DIR,'files','selected_models-EX.txt')

##### FUNCTIONS #####
plot_drug_recommendations = function(drug_treatments, drug_response, metadata){
    # Would we recommend some of the drugs given to patients if 
    # we would use our recommendations?
    
    ## treatments for each patient across cancer types
    treatments_clean = drug_treatments %>% 
        distinct(bcr_patient_barcode, cancer_type, DRUG_NAME)
    avail_drugs = treatments_clean %>% pull(DRUG_NAME) %>% unique()
    
    ## ranking predicted drug responses per sample
    drug_recs = drug_response %>%
        group_by(sample) %>%
        arrange(predicted_ic50) %>%
        mutate(ranking = row_number(),
               ranking_ratio = ranking / n(),
               bcr_patient_barcode = substr(sample, 1, 12)) %>%
        ungroup()
    
    ## for each patient, keep only treatment with the best ranking
    X = treatments_clean %>%
        left_join(drug_recs,
            by=c('bcr_patient_barcode','cancer_type','DRUG_NAME')) %>%
        drop_na() %>%
        distinct(drug_screen, cancer_type, sample, ranking, bcr_patient_barcode,
                 ranking_ratio, predicted_ic50, DRUG_NAME, DRUG_NAME_CLEAN) %>% 
        group_by(cancer_type, sample) %>% 
        # only consider the best case scenario
        slice_min(order_by=ranking_ratio, n=1) %>%
        ungroup() %>%
        left_join(metadata %>% distinct(sampleID, cancer_subtype), 
                  by=c("sample"="sampleID"))
    
    plts = list()
    plts[['drug_rec-sample_counts_by_cancer']] = X %>% 
        ungroup() %>% 
        count(cancer_type) %>% 
        ggbarplot(x="cancer_type", y="n", fill="cancer_type", color=NA, 
                  palette=get_palette("Paired", length(unique(X[["cancer_type"]]))), 
                  label=TRUE, lab.size=1) + 
        guides(fill="none") + 
        labs(x="Cancer Type", y="Count") + 
        coord_flip()
    
    plts[['drug_rec-sample_counts_by_cancer_and_treatment']] = X %>% 
        ungroup() %>% 
        group_by(cancer_type, DRUG_NAME_CLEAN) %>% 
        summarize(n = n()) %>%
        mutate(perc = n / sum(n)) %>%
        # filter(DRUG_NAME %in% X[["DRUG_NAME"]]) %>%
        ggballoonplot(x="cancer_type", y="DRUG_NAME_CLEAN", size="perc", fill="perc") +
        gradient_fill("red") +
        scale_size(range=c(0.5,3)) +
        labs(x="Cancer Type", y="Drug")
    
    plts[['drug_rec-drugs_ranking_ratios']] = X %>% 
        ggboxplot(x="cancer_type", y="ranking_ratio", 
                  fill="cancer_type", outlier.size=0.1, 
                  palette=get_palette("Paired", length(unique(X[["cancer_type"]])))) + 
        guides(fill="none") + 
        labs(x="Cancer Type", y="Ranking Ratio") + 
        coord_flip()
    
    cancers_w_subtypes = c("BRCA","READ","GBM","LGG","UCEC")
    plts[['drug_rec-sample_counts_by_cancer-by_subtype']] = X %>% 
        filter(!str_detect(cancer_subtype,"_STN") & !str_detect(cancer_subtype,"_NA")) %>%
        filter(cancer_type %in% cancers_w_subtypes) %>%
        ungroup() %>% 
        count(cancer_subtype) %>% 
        ggbarplot(x="cancer_subtype", y="n", fill="cancer_subtype", color=NA, 
                  palette=get_palette("Paired", 15), 
                  label=TRUE, lab.size=1) + 
        guides(fill="none") + 
        labs(x="Cancer Subtype", y="Count") + 
        coord_flip()
    
    plts[['drug_rec-drugs_ranking_ratios-by_subtype']] = X %>% 
        filter(!str_detect(cancer_subtype,"_STN") & !str_detect(cancer_subtype,"_NA")) %>%
        filter(cancer_type %in% cancers_w_subtypes) %>%
        ggboxplot(x="cancer_subtype", y="ranking_ratio", 
                  fill="cancer_subtype", outlier.size=0.1, 
                  palette=get_palette("Paired", 15)) + 
        guides(fill="none") + 
        labs(x="Cancer Subtype", y="Ranking Ratio") + 
        coord_flip()
    
    # combinations recommended drug cancer not tested?
    # patients with triple negative breast cancer were not given the drug 
    # we recommended, what would we recommend to them?
    cancer_oi = "KIRC"
    x = X %>% filter(cancer_type==cancer_oi)
    plts[["drug_rec-drugs_ranking_ratios-KIRC_given"]] = x %>% 
        ggboxplot(x="DRUG_NAME_CLEAN", y="ranking_ratio", outlier.size=0.1,
                  color="DRUG_NAME_CLEAN", 
                  palette=get_palette("Dark2",length(unique(x[["DRUG_NAME_CLEAN"]])))) + 
        labs(x="Drug", y="Ranking Ratio") + 
        guides(color="none") + 
        coord_flip()
    
    samples_oi = X %>% filter(cancer_type==cancer_oi) %>% pull(sample)
    x = lapply(1:9, function(nn){
          res = drug_recs %>%
            filter(sample%in%samples_oi & DRUG_NAME%in%avail_drugs) %>% 
            group_by(sample) %>% 
            arrange(ranking_ratio) %>%
            slice(nn) %>%
            ungroup()
        return(res)
    })
    x = do.call(rbind,x)
    x = x %>% 
        group_by(sample, DRUG_NAME_CLEAN) %>%
        slice_min(ranking, n=1) %>%
        ungroup() %>%
        group_by(sample) %>%
        arrange(ranking_ratio) %>%
        mutate(ranking=row_number())
    
    n_drugs = x %>% pull(DRUG_NAME_CLEAN) %>% unique() %>% length()
    drugs_given = treatments_clean %>% filter(cancer_type==cancer_oi) %>% pull(DRUG_NAME) %>% unique()
    drugs_recommended = x %>% pull(DRUG_NAME) %>% unique()
    new_recs = setdiff(drugs_recommended, drugs_given)
    plts[["drug_rec-drugs_ranking_ratios-KIRC_best"]] = x %>%
        mutate(lab = ifelse(DRUG_NAME %in% new_recs, 
                            paste0("*",DRUG_NAME_CLEAN), DRUG_NAME_CLEAN)) %>%
        ggboxplot(x="lab", y="ranking_ratio", outlier.size=0.1,
                  color="lab", palette=get_palette("Dark2", n_drugs)) + 
        labs(x="Drug", y="Ranking Ratio") + 
        guides(color="none") + 
        coord_flip() +
        facet_wrap(~ranking) +
        theme(strip.text.x = element_text(size=6, family='Arial'))
    
    return(plts)
}


make_plots = function(drug_treatments, drug_response, metadata){
    plts = list(
        plot_drug_recommendations(drug_treatments, drug_response, metadata)
    )
    plts = do.call(c,plts)
    return(plts)
}


save_plt = function(plts, plt_name, extension='.pdf', 
                    directory='', dpi=350, format=TRUE,
                    width = par("din")[1], height = par("din")[2]){
    print(plt_name)
    plt = plts[[plt_name]]
    if (format){
        plt = ggpar(plt, font.title=10, font.subtitle=10, font.caption=10, 
                    font.x=8, font.y=8, font.legend=8,
                    font.tickslab=6, font.family='Arial')    
    }
    filename = file.path(directory,paste0(plt_name,extension))
    save_plot(filename, plt, base_width=width, base_height=height, dpi=dpi, units='cm')
}


make_figdata = function(drug_treatments, drug_response, metadata){
    figdata = list(
        "drug_recommendation" = list(
            "drug_treatments" = drug_treatments,
            "drug_response" = drug_response,
            "metadata" = metadata
        )
    )
    return(figdata)
}


save_plots = function(plts, figs_dir){
    save_plt(plts, 'drug_rec-sample_counts_by_cancer', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'drug_rec-sample_counts_by_cancer_and_treatment', '.pdf', figs_dir, width=11, height=12)
    save_plt(plts, 'drug_rec-drugs_ranking_ratios', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'drug_rec-sample_counts_by_cancer-by_subtype', '.pdf', figs_dir, width=6, height=5)
    save_plt(plts, 'drug_rec-drugs_ranking_ratios-by_subtype', '.pdf', figs_dir, width=6, height=5)
    save_plt(plts, 'drug_rec-drugs_ranking_ratios-KIRC_given', '.pdf', figs_dir, width=7, height=5)
    save_plt(plts, 'drug_rec-drugs_ranking_ratios-KIRC_best', '.pdf', figs_dir, width=7, height=10)
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
    models_file = args$models_file
    figs_dir = args$figs_dir
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    drug_targets = read_tsv(drug_targets_file) %>% 
        mutate(DRUG_NAME=toupper(DRUG_NAME)) %>%
        dplyr::select(DRUG_ID,DRUG_NAME,TARGET,TARGET_PATHWAY) %>%
        distinct()
    drug_treatments = read_tsv(drug_treatments_file) %>% 
        mutate(DRUG_NAME=toupper(DRUG_NAME)) %>% 
        left_join(drug_targets, by="DRUG_NAME") %>% ## there's an overlap of 44
        mutate(pharmaceutical_tx_started_days_to = as.numeric(pharmaceutical_tx_started_days_to))
    drug_response = read_tsv(drug_response_file) %>%
        mutate(DRUG_ID=as.numeric(gsub("_.*","",ID))) %>%
        left_join(drug_targets %>% distinct(DRUG_ID,DRUG_NAME), by="DRUG_ID") %>%
        mutate(DRUG_NAME_CLEAN = ID) %>%
        separate(DRUG_NAME_CLEAN, c("drug_id","min_conc","max_conc"), sep="_") %>%
        mutate(DRUG_NAME_CLEAN = sprintf("%s (%s-%s)", DRUG_NAME, min_conc, max_conc))
        
    metadata_subtypes = read_tsv(metadata_subtypes_file) %>% 
        mutate(cancer_subtype=paste0(cancer,"_",cancer_subtype))
    metadata = read_tsv(metadata_file) %>%
        left_join(metadata_subtypes, by= c("sampleID", "cancer", "sample_type"))
    
    selected_drugs = readLines(selected_drugs_file)
    
    # subset for selected drugs
    drug_response = drug_response %>% filter(ID %in% selected_drugs)
    
    # make plots
    plts = make_plots(drug_treatments, drug_response, metadata)
    
    # make figdata
    figdata = make_figdata(drug_treatments, drug_response, metadata)
    
    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}