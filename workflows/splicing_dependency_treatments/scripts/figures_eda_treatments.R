#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# - clinical trial: https://clinicaltrials.gov/ct2/show/NCT01276574
# - NACT regimens (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6020442/; https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4206650/, https://www.cancer.org/cancer/ovarian-cancer/treating/chemotherapy.html): 
#   PACLITAXEL, DOCETAXEL, FLUOROURACIL, CYCLOPHOSPHAMIDE, CARBOPLATIN, CISPLATIN,
#   ALTRETAMINE, Capecitabine, Etoposide, Gemcitabine, Ifosfamide, Irinotecan,
#   DOXORUBICIN, MELPHALAN, PEMETREXED, TOPOTECAN, VINORELBINE
# - HGSC = high grade serous ovarian carcinoma; N/A = no information; a = death due other than cancer; NACT = neoadjuvant chemotherapy; PDS = primary debulking surgery; IDS = interval debulking surgery;
#   CP = platinum; DX = doxorubicin; TX = taxane; BV = bevacizumab; GEM = gemcitabine; 
#   PD = progressive disease; PR = partial response; SD = stable disease; b = ongoing follow-up; Relapse type = relapse type in ascites samples derived at progressive disease phase; IHC = immunohistochemistry; RNAseq = RNA sequencing; c = in vitro.
# 
# Outline
# -------
# 

require(optparse)
require(tidyverse)
require(ggpubr)
require(scattermore)
require(survival)
require(survminer)

ROOT = here::here()

# variables
NACT_DRUGS = c("PACLITAXEL", "DOCETAXEL", "5-FLUOROURACIL", 
               "CYCLOPHOSPHAMIDE", "CISPLATIN", "CARBOPLATIN", # carboplatin not found
               "ALTRETAMINE", "CAPECITABINE", "ETOPOSIDE", "GEMCITABINE", "IFOSFAMIDE", 
               "IRINOTECAN", "DOXORUBICIN", "MELPHALAN", "PEMETREXED", "TOPOTECAN", "VINORELBINE")
REGIMENS = c(
    "CISPLATIN" = "CP", # a.k.a. platinum (CISPLATIN or CARBOPLATIN)
    "CARBOPLATIN" = "CP",
    "DOXORUBICIN" = "DX",
    "PACLITAXEL" = "TX", # a.k.a. taxane (PACLITAXEL or DOCETAXEL)
    "DOCETAXEL" = "TX",
    "BEVACIZUMAB" = "BV", # a.k.a. bevacizumab
    "GEMCITABINE" = "GEM"
) %>% enframe("DRUG_NAME","chemotherapy_regimen")

TREATMENT_OUTCOMES = c(
    "Progressive Disease" = "PD",
    "Progressive Disease" = "PD after NACT",
    "Partial Response" = "PR",
    "Stable Disease" = "SD",
    "Chemotherapy Resistant" = "CR", # ?????
    "Progressive Disease" = "PDb" # on-going follow up
) %>% enframe("treatment_outcome_lab", "treatment_outcome")

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
# drug_screens_dir = file.path(PREP_DIR,'drug_screens')

# figs_dir = file.path(RESULTS_DIR,"figures","eda_treatments")

##### FUNCTIONS #####
load_drug_screens = function(drug_screens_dir){
    filenames = cbind(
        drug_screens_dir,
        expand_grid(c("train","test"),c("GDSC1.tsv.gz","GDSC2.tsv.gz")
                   )
    )
    dfs = apply(filenames, 1, function(x){
        x = as.vector(x)
        f = file.path(x[1],x[2],x[3])
        df = read_tsv(f)
        return(df)
    })
    dfs = do.call(rbind,dfs)
    return(dfs)
}


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
    
    # patient vs sample type vs treatment vs sample location vs regimen
    plts[["eda_metadata-sample_info-balloon"]] = X %>% 
        ggplot(aes(x=patientID, y=sample_type, group=patientID)) + 
        geom_point(aes(color=treatment, size=n)) + 
        geom_line(aes(color=treatment)) + 
        facet_grid(treatment~sample_location, scales="free_y") +
        labs(x="Patient", y="Sample Type", color="Treatment") + 
        coord_flip() + 
        theme_pubr(x.text.angle = 45)
    
    # responder vs non-responder
    X %>% count(sample_type)
    
    X %>% 
        group_by(patientID) %>%
        mutate(has_relapsed = any(sample_type == "Relapse")) %>%
        ungroup() %>%
        filter(sample_type=="Pre-chemo") %>%
        count(treatment, sample_location, has_relapsed)
    
    return(plts)
}


# differences between samples that responded or not
plot_drug_rec = function(metadata, estimated_response){
    X = metadata %>% 
        # add estimated_responses
        left_join(
            estimated_response,
            by=c("sampleID","DRUG_NAME")
        ) %>%
        drop_na(predicted_ic50) %>%
        group_by(
            patientID, treatment, sample_type,
            chemotherapy_regimen_lab, chemo_sensitivity, 
            drug_screen, PFI, PFI_status
        ) %>%
        summarize(predicted_ic50 = median(predicted_ic50)) %>%
        ungroup()
    
    plts = list()
    # Can we predict which patients will relapse before?
    plts[["drug_rec-"]] = X %>%
        filter(sample_type=="Pre-chemo") %>%
        ggplot(aes(x=chemo_sensitivity, y=predicted_ic50)) +
        #geom_violin(aes(fill=has_relapsed), color=NA) +
        geom_point(aes(color=chemotherapy_regimen_lab), position=position_jitter(0.2)) +
        geom_boxplot(width=0.1, outlier.shape=NA, outlier.size=0.1, fill=NA) +
        facet_wrap(~drug_screen, scales="free_y") +
        #guides(fill="none") +
        stat_compare_means(method="wilcox.test") +
        labs(x="Treatment Outcome", y="Predicted log(IC50)", color="Regimen") +
        theme_pubr(x.text.angle = 45)
    
    x = X %>%
        group_by(drug_screen, treatment, sample_type) %>%
        mutate(pred_sensitivity = case_when(
            predicted_ic50 < quantile(predicted_ic50, 0.10) ~ "Sensitive",
            predicted_ic50 > quantile(predicted_ic50, 0.90) ~ "Resistant",
            TRUE ~ "Ambiguous"
        )) %>%
        filter(sample_type=="Pre-chemo" & treatment=="PDS" & 
               drug_screen=="GDSC1" & pred_sensitivity!="Ambiguous")
    fit = survfit(Surv(PFI, PFI_status) ~ pred_sensitivity, data=x)
    ggsurvplot(fit, data=x, conf.int=TRUE, pval=TRUE, risk.table=TRUE)
    
    # deconvolve which exons determined the differential sensitivity: diff splicing of each cancer-driver exon vs Pearson correlation of their model.
    
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

    drug_screen = load_drug_screens(drug_screens_dir)
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
            relapse_type = `Relapse type`,
            chemotherapy_regimen = `Primary chemo regimen`,
            PFS = `Progression free survival / PFS (months)`,
            PFI = `Platinum free interval / PFI (days)`,
            OS = `Overall survival (months)`
        ) %>%
        mutate(
            sample_type = factor(sample_type, levels = c("Pre-chemo","Post-chemo","Relapse")),
            chemotherapy_regimen_lab = chemotherapy_regimen,
            PFI_status = ifelse(str_detect(PFI, "b"), 0, 1),
            PFI = as.numeric(gsub("b","", PFI)),
            chemo_sensitivity = ifelse(PFI<=(30*6), "Resistant", "Sensitive"),
        ) %>%
        separate_rows(sampleID, sep="; ") %>%
        separate_rows(chemotherapy_regimen, sep="-") %>%
        drop_na(sampleID) %>%
        left_join(REGIMENS, by="chemotherapy_regimen") %>%
        left_join(TREATMENT_OUTCOMES, by="treatment_outcome") %>%
        distinct()
    
    estimated_response = estimated_response %>%
        left_join(
            drug_screen %>% distinct(DRUG_NAME, DRUG_ID, ID),
            by="ID"
        ) %>%
        rowwise() %>%
        mutate(sampleID = paste(strsplit(sample, "_")[[1]][2:5], collapse="_")) %>%
        ungroup()    
    
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