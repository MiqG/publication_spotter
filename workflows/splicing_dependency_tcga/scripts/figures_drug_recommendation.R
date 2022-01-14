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
require(survival)
require(gtools)
require(extrafont)

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

# variables
THRESH_LR_FDR = 0.1

# Development
# ----------- 
# PREP_DIR = file.path(ROOT,'data','prep')
# RAW_DIR = file.path(ROOT,'data','raw')
# RESULTS_DIR = file.path(ROOT,'results','splicing_dependency_tcga')
# models_file = file.path(ROOT,'results','splicing_dependency_drugs','files','model_summaries_drug_response-EX.tsv.gz')
# spldep_file = file.path(RESULTS_DIR,'files','LGG','splicing_dependency-EX','mean.tsv.gz')
# drug_targets_file = file.path(RAW_DIR,'GDSC','screened_compunds_rel_8.2.csv')
# metadata_file = file.path(PREP_DIR,'metadata','PANCAN.tsv.gz')
# drug_treatments_file = file.path(PREP_DIR,'drug_treatments','PANCAN.tsv.gz')
# figs_dir = file.path(RESULTS_DIR,'figures','drug_recommendation')
# drug_response_file = file.path(RESULTS_DIR,'files','PANCAN','estimated_drug_response_by_drug-EX.tsv.gz')
# metadata_subtypes_file = file.path(PREP_DIR,'metadata','PANCAN_subtypes.tsv.gz')

##### FUNCTIONS #####
plot_drug_recommendations = function(drug_treatments, drug_response, metadata){
    # Would we recommend some of the drugs given to patients if 
    # we would use our recommendations?
    
    ## treatments for each patient across cancer types
    treatments_clean = drug_treatments %>% 
        distinct(bcr_patient_barcode, cancer_type, DRUG_NAME)
    
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
                 ranking_ratio, predicted_ic50, DRUG_NAME) %>% 
        group_by(cancer_type, sample, DRUG_NAME) %>% 
        # only consider the best case scenario
        slice_min(order_by=ranking_ratio, n=1) %>%
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
        group_by(cancer_type, DRUG_NAME) %>% 
        summarize(n = n()) %>%
        mutate(perc = n / sum(n)) %>%
        # filter(DRUG_NAME %in% X[["DRUG_NAME"]]) %>%
        ggballoonplot(x="cancer_type", y="DRUG_NAME", size="perc", fill="perc") +
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
    
    ## Does the drug recommendation ranking associate with patient survival?
    x = X %>% 
        ungroup() %>%
        left_join(metadata %>% distinct(sampleID, PFI, PFI.time),
                  by=c("sample"="sampleID"))
    
    result = lapply(unique(X[["cancer_type"]]), function(cancer_oi){
        test = coxph(Surv(PFI.time, PFI) ~ ranking_ratio, 
                     data=x %>% filter(cancer_type==cancer_oi))
        test = summary(test)$coefficients
        df = data.frame(
            cancer_type = cancer_oi,
            cox_coef = test["ranking_ratio","coef"],
            cox_std = test["ranking_ratio","se(coef)"],
            cox_zscore = test["ranking_ratio","z"],
            cox_pvalue = test["ranking_ratio","Pr(>|z|)"]
        )
        return(df)
    })
    result = do.call(rbind, result)
    
    plts[['drug_rec-drug_ranking_coxph']] = result %>% 
        mutate(lab_pvals = stars.pval(cox_pvalue)) %>%
        ggplot(aes(x=cancer_type, y=cox_coef, color=cancer_type)) + 
        geom_point(size=0.5) + 
        geom_errorbar(aes(ymin=cox_coef - cox_std, ymax=cox_coef + cox_std), size=0.25) + 
        geom_text(aes(y=2.5, label=lab_pvals), color="black", size=2) +
        theme_pubr() + 
        geom_hline(yintercept=0, linetype='dashed') + 
        color_palette(get_palette("Paired", length(unique(X[["cancer_type"]])))) + 
        guides(color="none") + 
        labs(x="Cancer Type", y="Cox Coefficient") +
        coord_flip()

    # are ranking ratios associated with the order of the treatment?
    x = treatments_clean %>%
        left_join(drug_recs,
            by=c('bcr_patient_barcode','cancer_type','DRUG_NAME')) %>%
        drop_na() %>%
        distinct(drug_screen, cancer_type, sample, ranking, bcr_patient_barcode,
                 ranking_ratio, predicted_ic50, DRUG_NAME) %>% 
        left_join(metadata %>% distinct(sampleID, cancer_subtype), 
                  by=c("sample"="sampleID"))
    
    treatments_prep = drug_treatments %>% 
        distinct(bcr_patient_barcode, cancer_type, 
                 DRUG_NAME, pharmaceutical_tx_started_days_to) %>%
        # sort by how late the treatment started
        drop_na() %>%
        arrange(cancer_type, bcr_patient_barcode, 
                pharmaceutical_tx_started_days_to) %>%
        # if a treatment was given more than once, keep the earliest it was given
        group_by(bcr_patient_barcode, DRUG_NAME) %>%
        slice_min(pharmaceutical_tx_started_days_to, n=1) %>%
        ungroup() %>%
        # for each patient, in which order was the drug given?
        group_by(cancer_type, bcr_patient_barcode) %>%
        mutate(tx_order = row_number(),
               tx_order_ratio = tx_order / n(),
               tx_mult = n()>1,
               tx_n = n())
    x = X %>%
        ungroup() %>%
        left_join(treatments_prep, by=c("bcr_patient_barcode", "DRUG_NAME", "cancer_type"))
    x %>% filter(cancer_type=="BRCA") %>% ggscatter(x="tx_order", y="ranking_ratio", color="DRUG_NAME", size="tx_n") + stat_cor(method="spearman")# + facet_wrap(~cancer_subtype)
    
    # combinations recommended drug cancer not tested?
    
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
    save_plt(plts, 'drug_rec-sample_counts_by_cancer', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'drug_rec-sample_counts_by_cancer_and_treatment', '.pdf', figs_dir, width=11, height=12)
    save_plt(plts, 'drug_rec-drugs_ranking_ratios', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'drug_rec-sample_counts_by_cancer-by_subtype', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'drug_rec-drugs_ranking_ratios-by_subtype', '.pdf', figs_dir, width=5, height=5)
    save_plt(plts, 'drug_rec-drug_ranking_coxph', '.pdf', figs_dir, width=5, height=5)
}


make_figdata = function(){
    figdata = list(
    )
    return(figdata)
}


save_figdata = function(figdata, dir){
    lapply(names(figdata), function(x){
        filename = file.path(dir,'figdata',paste0(x,'.xlsx'))
        dir.create(dirname(filename), recursive=TRUE)
        write_xlsx(figdata[[x]], filename)
    })
}


main = function(){
    args = getParsedArgs()
    models_file = args$models_file
    figs_dir = args$figs_dir
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
#     models = read_tsv(models_file) %>% 
#         mutate(event_gene = paste0(EVENT,'_',GENE),
#                event_type = gsub('Hsa','',gsub("[^a-zA-Z]", "",EVENT)))
    drug_targets = read_csv(drug_targets_file) %>% 
        mutate(DRUG_NAME=toupper(DRUG_NAME)) %>%
        dplyr::select(DRUG_ID,DRUG_NAME,TARGET,TARGET_PATHWAY) %>%
        separate_rows(TARGET) %>%
        distinct()
    drug_treatments = read_tsv(drug_treatments_file) %>% 
        mutate(DRUG_NAME=toupper(DRUG_NAME)) %>% 
        left_join(drug_targets, by="DRUG_NAME") %>% ## there's an overlap of 44
        mutate(pharmaceutical_tx_started_days_to = as.numeric(pharmaceutical_tx_started_days_to))
    drug_response = read_tsv(drug_response_file) %>%
        left_join(drug_targets %>% distinct(DRUG_ID,DRUG_NAME), by="DRUG_ID")
    metadata_subtypes = read_tsv(metadata_subtypes_file) %>% 
        mutate(cancer_subtype=paste0(cancer,"_",cancer_subtype))
    metadata = read_tsv(metadata_file) %>%
        left_join(metadata_subtypes, by= c("sampleID", "cancer", "sample_type"))
    
#     spldep = read_tsv(spldep_file)
    
    # make plots
    plts = make_plots(drug_treatments, drug_response, metadata)
    
    # make figdata
    # figdata = make_figdata(results_enrich)
    
    # save
    save_plots(plts, figs_dir)
    # save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}