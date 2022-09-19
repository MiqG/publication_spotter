#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# 
# 
# Notes
# -------
# - each patient was treated with a regimen.
# - clinical trial: https://clinicaltrials.gov/ct2/show/NCT01276574
# - NACT regimens (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6020442/; https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4206650/, https://www.cancer.org/cancer/ovarian-cancer/treating/chemotherapy.html): 
#   PACLITAXEL, DOCETAXEL, FLUOROURACIL, CYCLOPHOSPHAMIDE, CARBOPLATIN, CISPLATIN,
#   ALTRETAMINE, Capecitabine, Etoposide, Gemcitabine, Ifosfamide, Irinotecan,
#   DOXORUBICIN, MELPHALAN, PEMETREXED, TOPOTECAN, VINORELBINE
# - HGSC = high grade serous ovarian carcinoma; N/A = no information; a = death due other than cancer; NACT = neoadjuvant chemotherapy; PDS = primary debulking surgery; IDS = interval debulking surgery;
#   CP = platinum; DX = doxorubicin; TX = taxane; BV = bevacizumab; GEM = gemcitabine; 
#   PD = progressive disease; PR = partial response; SD = stable disease; b = ongoing follow-up; Relapse type = relapse type in ascites samples derived at progressive disease phase; IHC = immunohistochemistry; RNAseq = RNA sequencing; c = in vitro.

require(optparse)
require(tidyverse)
require(ggpubr)
require(ggrepel)
require(tidytext)
require(cowplot)
require(extrafont)

# variables
# NACT_DRUGS = c("PACLITAXEL", "DOCETAXEL", "5-FLUOROURACIL", 
#                "CYCLOPHOSPHAMIDE", "CISPLATIN", "CARBOPLATIN", # carboplatin not found
#                "ALTRETAMINE", "CAPECITABINE", "ETOPOSIDE", "GEMCITABINE", "IFOSFAMIDE", 
#                "IRINOTECAN", "DOXORUBICIN", "MELPHALAN", "PEMETREXED", "TOPOTECAN", "VINORELBINE")
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

THRESH_PFI = 30*6
THRESH_FDR = 0.05
THRESH_DPSI = 5

# formatting
PAL_SINGLE_ACCENT = "#CEB5B7" # #F5E0B7 #DBD3AD #9AADBF
PAL_SINGLE_LIGHT = "#bc6c25"
PAL_SINGLE_DARK = "#606c38"
PAL_RANGE = c(PAL_SINGLE_DARK, PAL_SINGLE_ACCENT, PAL_SINGLE_LIGHT)
PAL_REGIMEN = setNames(
    get_palette(PAL_RANGE, 4),
    c("CP","CP-GEM-BV","CP-TX","CP-TX-BV")
)
PAL_IS_TARGET = get_palette(c("#716454","#6AC2BF"), 4)[3:4]

LINE_SIZE = 0.25
FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"


# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,"data","raw")
# PREP_DIR = file.path(ROOT,"data","prep")
# RESULTS_DIR = file.path(ROOT,'results','splicing_dependency_treatments')
# MODELS_DIR = file.path(ROOT,"results","splicing_dependency_drugs")
# drug_targets_file = file.path(PREP_DIR,'drug_screens','drug_targets.tsv.gz')
# protein_impact_file = file.path(ROOT,"data","raw","VastDB","PROT_IMPACT-hg38-v3.tab.gz")
# event_info_file = file.path(RAW_DIR,"VastDB","EVENT_INFO-hg38_noseqs.tsv")
# gene_info_file = file.path(RAW_DIR,"ENSEMBL","gene_annotation-hg38.tsv.gz")
# metadata_file = file.path(PREP_DIR,"metadata","Zhang2022.tsv.gz")
# spldep_file = file.path(RESULTS_DIR,'files','Zhang2022','splicing_dependency-EX','mean.tsv.gz')
# splicing_file = file.path(PREP_DIR,"event_psi","Zhang2022-EX.tsv.gz")
# estimated_response_file = file.path(RESULTS_DIR,"files","Zhang2022","estimated_drug_response_by_drug-EX.tsv.gz")
# drug_screens_dir = file.path(PREP_DIR,'drug_screens')

# drug_models_file = file.path(MODELS_DIR,"files","model_summaries_drug_response-EX.tsv.gz")
# drug_models_selected_file = file.path(MODELS_DIR,"files","selected_models-EX.txt")

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
        
    plts = list()
    
    # patients PFI
    plts[["eda_metadata-patients-pfi"]] = X %>%
        distinct(patientID, PFI) %>% 
        gghistogram(x="PFI", fill=PAL_SINGLE_DARK, bins=20, color=NA) +
        geom_vline(xintercept=THRESH_PFI, linetype="dashed", size=LINE_SIZE) +
        labs(x="Progression-Free Interval (PFI)", y="Count") +
        theme(aspect.ratio=NULL)
    
    # patients vs sensitivity vs regimen
    plts[["eda_metadata-patients-sensitivity"]] = X %>%
        distinct(patientID, chemo_sensitivity, chemotherapy_regimen_lab, sample_type) %>%
        group_by(sample_type) %>%
        mutate(sample_type_lab=sprintf("%s (n=%s)", sample_type, n())) %>%
        ungroup() %>%
        count(chemo_sensitivity, chemotherapy_regimen_lab, sample_type_lab) %>%
        ggbarplot(
            x="chemo_sensitivity", y="n", color=NA,
            fill="chemotherapy_regimen_lab", position=position_dodge2(0.7, preserve="single"), 
            label=TRUE, lab.font=FONT_FAMILY, lab.size=FONT_SIZE, lab.hjust=0.5,
            palette=PAL_REGIMEN
        ) +
        facet_wrap(~sample_type_lab) +
        labs(x="Sensitivity to Chemotherapy", y="Count", fill="Regimen") +
        theme(aspect.ratio=NULL, strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    return(plts)
}


plot_drug_rec = function(metadata, estimated_response){
    
    X = metadata %>% 
        # add estimated_responses
        left_join(
            estimated_response,
            by=c("sampleID","DRUG_NAME")
        ) %>%
        drop_na(predicted_ic50) %>%
        # summarize patients with multiple treatments
        group_by(
            patientID, sample_type,
            chemotherapy_regimen_lab, chemo_sensitivity, 
            drug_screen, PFI
        ) %>%
        summarize(predicted_ic50 = median(predicted_ic50)) %>%
        ungroup()
    
    plts = list()
    # Can we predict which patients will relapse before?
    plts[["drug_rec-pred_ic50-boxplot"]] = X %>%
        filter(sample_type=="Pre-chemo") %>%
        ggplot(aes(x=chemo_sensitivity, y=predicted_ic50)) +
        geom_point(aes(color=chemotherapy_regimen_lab, shape=chemotherapy_regimen_lab), 
                   position=position_jitter(0.2), size=1) +
        geom_boxplot(width=0.2, outlier.shape=NA, outlier.size=0.1, fill=NA) +
        facet_wrap(~drug_screen, scales="free_y") +
        stat_compare_means(method="wilcox.test", family=FONT_FAMILY, size=FONT_SIZE) +
        color_palette(palette=PAL_REGIMEN) +
        labs(x="Treatment Outcome", y="Predicted log(IC50)", 
             color="Regimen", shape="Regimen") +
        theme_pubr(x.text.angle = 0) +
        theme(
            aspect.ratio = 1,
            strip.text.x = element_text(size=6, family=FONT_FAMILY)
        )
    
    return(plts)
}


plot_diff_sens = function(diff_splicing, event_info, protein_impact, drug_models, metadata){
    # deconvolve which exons determined the differential sensitivity: 
    # diff splicing of each cancer-driver exon vs Pearson correlation of their model.
    X = diff_splicing %>%
        mutate(log10_padj = -log10(p.adj)) %>%
        left_join(event_info %>% distinct(EVENT, GENE, event_gene), by="EVENT") %>%
        left_join(protein_impact, by="EVENT")
    
    plts = list()
    
    plts[["diff_sens-res_vs_sens-scatter"]] = X %>%
        ggplot(aes(x=delta_psi, y=log10_padj)) +
        geom_point(size=1, color=PAL_SINGLE_DARK, alpha=0.5) +
        geom_hline(yintercept=-log10(THRESH_FDR), size=LINE_SIZE, linetype="dashed") +
        geom_vline(xintercept=0, size=LINE_SIZE, linetype="dashed") +
        geom_text_repel(
            aes(label=event_gene), . %>% filter(p.adj<THRESH_FDR & abs(delta_psi)>THRESH_DPSI),
            size=FONT_SIZE, family=FONT_FAMILY, segment.size=LINE_SIZE, min.segment.length=0,
            force=10
        ) + 
        labs(x="Delta PSI", y="-log10(FDR)", subtitle="Resistant vs. Sensitive") +
        theme_pubr()
        
    
    x = X %>% 
        filter(p.adj<THRESH_FDR & abs(delta_psi)>THRESH_DPSI) %>% 
        left_join(
            drug_models %>% 
                filter(DRUG_NAME %in% metadata[["DRUG_NAME"]]), 
            by=c("EVENT","GENE")
        ) %>%
        left_join(
            drug_targets %>% distinct(DRUG_ID, TARGET) %>% mutate(is_target=TRUE),
            by = c("DRUG_ID","GENE"="TARGET")
        ) %>%
        mutate(is_target=replace_na(is_target, FALSE)) %>%
        left_join(metadata %>% distinct(DRUG_NAME, chemotherapy_regimen_lab), by="DRUG_NAME")
    
    plts[["diff_sens-event_contribution-bar"]] = x %>%
        group_by(event_gene, chemotherapy_regimen_lab, drug_screen, is_target) %>%
        summarize(pearson_correlation_median = median(pearson_correlation)) %>%
        ungroup() %>%
        ggbarplot(x="event_gene", y="pearson_correlation_median", color=NA,
                  fill="chemotherapy_regimen_lab", position=position_dodge(0.8)) +
        fill_palette(palette=PAL_REGIMEN) +
        facet_wrap(~drug_screen) +
        coord_flip() +
        labs(x="Event & Gene", y="median(Exon Contribution)", fill="Regimen") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    return(plts)
}


make_plots = function(metadata, estimated_response, diff_splicing, 
                      event_info, protein_impact, drug_models){
    plts = list(
        plot_eda_metadata(metadata),
        plot_drug_rec(metadata, estimated_response),
        plot_diff_sens(diff_splicing, event_info, protein_impact, drug_models, metadata)
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
    save_plt(plts, "eda_metadata-patients-pfi", ".pdf", figs_dir, width=3, height=4)
    save_plt(plts, "eda_metadata-patients-sensitivity", ".pdf", figs_dir, width=8, height=6)
    save_plt(plts, "drug_rec-pred_ic50-boxplot", ".pdf", figs_dir, width=5, height=5)
    save_plt(plts, "diff_sens-res_vs_sens-scatter", ".pdf", figs_dir, width=5, height=6)
    save_plt(plts, "diff_sens-event_contribution-bar", ".pdf", figs_dir, width=6.2, height=5.5)
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

    splicing = read_tsv(splicing_file)
    estimated_response = read_tsv(estimated_response_file)
    drug_models = read_tsv(drug_models_file)
    drug_screen = load_drug_screens(drug_screens_dir)
    drug_targets = read_tsv(drug_targets_file)
    
    protein_impact = read_tsv(protein_impact_file)
    event_info = read_tsv(event_info_file)
    gene_info = read_tsv(gene_info_file)
    
    # prep inputs
    event_info = event_info %>%
        mutate(event_gene = paste0(EVENT,"_",GENE))
    
    protein_impact = protein_impact %>%
            dplyr::rename(EVENT=EventID, term=ONTO) %>%
            dplyr::select(term,EVENT) %>%
            mutate(term_clean=gsub(" \\(.*","",term),
                   term_clean=gsub("ORF disruption upon sequence exclusion",
                                   "ORF disruption (exclusion)",term_clean),
                   term_clean=gsub("ORF disruption upon sequence inclusion",
                                   "ORF disruption (inclusion)",term_clean),
                   term_clean=gsub("In the CDS, with uncertain impact",
                                   "In the CDS (uncertain)",term_clean))
    
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
            chemo_sensitivity = ifelse(PFI <= THRESH_PFI, "Resistant", "Sensitive"),
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
    
    drug_models = drug_models %>%
        left_join(drug_screen %>% distinct(DRUG_NAME, DRUG_ID, ID), by=c("ID","DRUG_ID"))
    
    # differential splicing cancer-driver events: resistant vs sensitive patients
    events_oi = drug_models %>% pull(EVENT) %>% unique()
    X = splicing %>% 
        filter(EVENT %in% events_oi) %>%
        pivot_longer(-EVENT, names_to="sample", values_to="psi") %>%
        rowwise() %>%
        mutate(sampleID = paste(strsplit(sample, "_")[[1]][2:5], collapse="_")) %>%
        ungroup() %>%
        left_join(
            metadata %>% distinct(sampleID, patientID, chemo_sensitivity, sample_type),
            by = "sampleID"
        ) %>%
        drop_na(chemo_sensitivity, psi) %>% 
        filter(sample_type=="Pre-chemo")
    
    diff_splicing = compare_means(
            psi ~ chemo_sensitivity,
            X,
            group.by = "EVENT",
            method = "wilcox.test",
            p.adjust.method = "fdr"
        ) %>% left_join(
            X %>% 
            group_by(EVENT, chemo_sensitivity) %>%
            summarize(psi_med = median(psi)) %>%
            ungroup() %>%
            pivot_wider(id_cols=EVENT, names_from=chemo_sensitivity, 
                        names_prefix = "psi_med_", values_from=psi_med) %>%
            mutate(delta_psi = psi_med_Resistant - psi_med_Sensitive),
            by = "EVENT"
        )
    
    # make plots
    plts = make_plots(metadata, estimated_response, diff_splicing, 
                      event_info, protein_impact, drug_models)
    
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