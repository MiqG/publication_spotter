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
require(extrafont)
require(scattermore)
require(ComplexHeatmap)
require(stringi)

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

# Development
# -----------
# study_oi = "Kumar2022"
# study_oi = "Szenajch2020"
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# RESULTS_DIR = file.path(ROOT,'results','splicing_dependency_personalization')
# MODELS_DIR = file.path(ROOT,'results','splicing_dependency_drugs')
# metadata_file = file.path(RAW_DIR,'articles',study_oi,'metadata.tsv')

# spldep_file = file.path(RESULTS_DIR,'files',study_oi,'splicing_dependency-EX','mean.tsv.gz')
# psi_file = file.path(RAW_DIR,'articles',study_oi,'vast_out','PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz')
# drug_response_file = file.path(RESULTS_DIR,'files',study_oi,'estimated_drug_response_by_drug-EX.tsv.gz')

# selected_events_file = file.path(ROOT,'results','model_splicing_dependency','files','selected_models-EX.txt')
# event_info_file = file.path(RAW_DIR,'VastDB','EVENT_INFO-hg38_noseqs.tsv')
# selected_drugs_file = file.path(MODELS_DIR,'files','selected_models-EX.txt')
# drug_targets_file = file.path(PREP_DIR,'drug_screens','drug_targets.tsv.gz')

# figs_dir = file.path(RESULTS_DIR,'figures','personalization',study_oi)


##### FUNCTIONS #####
plot_Szenajch2020 = function(harm, drug_response, event_annot){
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7730278/
    # ovarian cancer cell line A2780 treated with paclitaxel (PTX) for 
    # check how vulnerabilities change across time
    # check how drug recommendations change across time
    # in their paper, paclitaxel-resistant cell lines become sensitive to cisplatin
    # check agreement between samples
       
    
    samples_oi = c("Control","4nM","8nM","16nM","32nM","64nM","128nM")
    
    # from table 1
    ic50 = data.frame(
        treatment_type = c(samples_oi,samples_oi),
        logIC50 = log(c(c(4.62, 11.47, 18.71, 77.54, 148.47, 173.86, 272.47),
                        c(4.58, 4.27,  2.12,  0.65,  0.44,   0.4,    0.83))),
        DRUG_NAME = c(rep("PACLITAXEL",7),rep("CISPLATIN",7))
    )
    
    plts = list()
    
    
    # drug recommendation
    ## real
    X = ic50 %>% 
        pivot_wider(id_cols="treatment_type", 
                    names_from="DRUG_NAME", 
                    values_from="logIC50")
    X %>% ggscatter(x="PACLITAXEL", y="CISPLATIN", label="treatment_type") + labs(x="Real Paclitaxel log(IC50)", y="Real Cisplatin log(IC50)") + stat_cor()
    
    
    X = drug_response %>%
        mutate(treatment_type = gsub("PTX","nM",treatment_type), 
               lab=paste0(treatment_type,"_",replicate)) %>%
        left_join(ic50, by=c("treatment_type","DRUG_NAME")) %>%
        group_by(lab) %>%
        # normalize predictions per sample
        mutate(zscore = scale(predicted_ic50)) %>%
        ungroup()
    
    drugs_oi = c("PACLITAXEL","CISPLATIN")
    
    # normalize
    X %>% ggboxplot(x="lab", y="predicted_ic50")
    X %>% ggboxplot(x="lab", y="zscore")
    
    # correspondence predictions
    X %>% 
        pivot_wider(id_cols="ID", names_from="lab", values_from="zscore") %>% 
        column_to_rownames("ID") %>% 
        cor(use="pairwise.complete.obs", method="pearson") %>% 
        Heatmap()
    
    # pred vs real
    X %>% filter(DRUG_NAME %in% drugs_oi) %>% ggerrorplot(x="logIC50", y="zscore", numeric.x.axis=TRUE, label="treatment_type", repel=TRUE) + facet_wrap(~DRUG_NAME, scales="free") + labs(x="Real log(IC50)", y="Pred. IC50") + stat_cor()
    
    X %>% filter(DRUG_NAME %in% drugs_oi) %>% group_by(DRUG_NAME, treatment_type) %>% summarize(zscore=mean(zscore)) %>% pivot_wider(id_cols = c("treatment_type"), names_from="DRUG_NAME", values_from="zscore", values_fn=mean) %>% ggscatter(x="PACLITAXEL", y="CISPLATIN", alpha=0.5, label="treatment_type", repel=TRUE) + stat_cor() + labs(x="Pred. Paclitaxel log(IC50)", y="Pred. Cisplatin log(IC50)")
    
    
    
    plts[["drug_rec-preds_ic50-examples"]] = X %>%
        filter(DRUG_NAME %in% drugs_oi) %>%
        ggboxplot(x="treatment_type", y="zscore",fill="DRUG_NAME", 
                  order=samples_oi, palette="rickandmorty", outlier.shape=NA,
                  add="point", add.params=list(size=0.5, group="DRUG_NAME")) +
        theme_pubr(x.text.angle = 70) +
        labs(x="Paclitaxel-Adapted Cell Line", y="Pred. log(IC50) Scaled", fill="Drug") +
        facet_wrap(~DRUG_NAME, scales="free_y")
    
    
    drugs_oi = X %>%
        group_by(sample) %>%
        slice_min(zscore, n=5) %>%
        ungroup() %>%
        pull(DRUG_NAME) %>% 
        unique()
    
    plts[["drug_rec-preds_ic50-top5"]] = X %>%
        filter(DRUG_NAME %in% drugs_oi) %>%
        ggboxplot(x="treatment_type", y="zscore",fill="DRUG_NAME", 
                  order=samples_oi, palette="rickandmorty", 
                  add="point", add.params=list(size=0.5, group="DRUG_NAME")) +
        theme_pubr(x.text.angle = 70) +
        labs(x="Paclitaxel-Adapted Cell Line", y="Pred. log(IC50) Scaled", fill="Drug")
    
    # vulnerabilities
    ## spldep
    events_oi = harm %>%
        group_by(sampleID) %>%
        slice_min(spldep, n=10) %>%
        ungroup() %>%
        pull(index) %>% 
        unique()
    
    harm %>% filter(index%in%events_oi) %>% mutate(treatment_type = gsub("PTX","nM",treatment_type), treatment_type=factor(treatment_type, levels=samples_oi), lab=paste0(treatment_type,"_",replicate)) %>% pivot_wider(id_cols="index", names_from="lab", values_from="spldep") %>% column_to_rownames("index") %>% filter_all(any_vars(!is.na(.))) %>% Heatmap(name="Spl. Dep.")
    
    harm %>% mutate(treatment_type = gsub("PTX","nM",treatment_type), treatment_type=factor(treatment_type, levels=samples_oi), lab=paste0(treatment_type,"_",replicate)) %>% pivot_wider(id_cols="index", names_from="lab", values_from="psi") %>% column_to_rownames("index") %>% cor(use="pairwise.complete.obs") %>% Heatmap()
    
    ## harm
    events_oi = harm %>%
        group_by(sampleID) %>%
        slice_min(harm_score, n=10) %>%
        ungroup() %>%
        pull(index) %>% 
        unique()
    
    harm %>% filter(index%in%events_oi) %>% mutate(treatment_type = gsub("PTX","nM",treatment_type), treatment_type=factor(treatment_type, levels=samples_oi), lab=paste0(treatment_type,"_",replicate)) %>% pivot_wider(id_cols="index", names_from="lab", values_from="harm_score") %>% column_to_rownames("index") %>% filter_all(any_vars(!is.na(.))) %>% Heatmap(name="Max. Harm")
    
    return(plts)
}


plot_Kumar2021 = function(harm, drug_response, event_annot){
    # https://www.biorxiv.org/content/10.1101/2021.11.26.470089v2.full
    # they evolved clones resistant to irinotecan
    # check agreement between samples
    
    samples_oi = c("SCC1","MSC1","MSC2","MSC3")
    plts = list()
    
    
    # vulnerabilities
    ## spldep
    events_oi = harm %>%
        group_by(sampleID) %>%
        slice_min(spldep, n=10) %>%
        ungroup() %>%
        pull(index) %>% 
        unique()
    
    harm %>% filter(index%in%events_oi) %>% pivot_wider(id_cols="index", names_from="description", values_from="spldep") %>% column_to_rownames("index") %>% filter_all(any_vars(!is.na(.))) %>% Heatmap(name="Spl. Dep.")
    
    ## harm
    events_oi = harm %>%
        group_by(sampleID) %>%
        slice_min(harm_score, n=10) %>%
        ungroup() %>%
        pull(index) %>% 
        unique()
    
    harm %>% group_by(sampleID) %>% arrange(harm_score) %>% mutate(harm_rank = row_number()) %>% ungroup() %>% filter(index%in%events_oi) %>% pivot_wider(id_cols="index", names_from="description", values_from="harm_rank") %>% column_to_rownames("index") %>% filter_all(any_vars(!is.na(.))) %>% Heatmap(name="Max. Harm")
    
    # drug recommendation
    X = drug_response    
    drugs_oi = c("IRINOTECAN")
    
    plts[["drug_rec-preds_ic50-examples"]] = X %>%
        filter(DRUG_NAME %in% drugs_oi) %>%
        ggboxplot(x="strain", y="predicted_ic50",fill="DRUG_NAME", 
                  order=samples_oi, palette="rickandmorty", 
                  add="point", add.params=list(size=0.5, group="DRUG_NAME")) +
        theme_pubr(x.text.angle = 70) +
        labs(x="Irinotecan-Adapted Cell Line", y="Pred. log(IC50)", fill="Drug")
    
    drugs_oi = drug_response %>%
        group_by(sample) %>%
        slice_min(predicted_ic50, n=10) %>%
        ungroup() %>%
        pull(DRUG_NAME_CLEAN) %>% 
        unique()
    
    plts[["drug_rec-preds_ic50-top5"]] = X %>%
        filter(DRUG_NAME_CLEAN %in% drugs_oi) %>%
        ggboxplot(x="strain", y="predicted_ic50",fill="DRUG_NAME", outlier.shape = NA,
                  order=samples_oi, palette="rickandmorty",
                  add="point", add.params=list(size=0.5, group="DRUG_NAME")) +
        theme_pubr(x.text.angle = 70) +
        labs(x="Paclitaxel-Adapted Cell Line", y="Pred. log(IC50)", fill="Drug")
    
    return(plts)
}


make_plots = function(embedding, metadata, diff_result, harm, 
                      drug_response, event_annot){
    plts = list(
        plot_personalization_example(diff_result, harm, 
                                     drug_response, event_annot)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(embedding){
    figdata = list(
        "personalization" = list(
            "sample_clustering" = embedding
        )
    )
    return(figdata)
}


save_plt = function(plts, plt_name, extension='.pdf', 
                    directory='', dpi=350, format=TRUE,
                    width = par("din")[1], height = par("din")[2]){
    print(plt_name)
    plt = plts[[plt_name]]
    if (format){
        plt = ggpar(plt, font.title=10, font.subtitle=10, font.caption=10, 
                    font.x=8, font.y=8, font.legend=8,
                    font.tickslab=6, font.family='Arial', device=cairo_pdf)    
    }
    filename = file.path(directory,paste0(plt_name,extension))
    save_plot(filename, plt, base_width=width, base_height=height, dpi=dpi, units='cm')
}


save_plots = function(plts, figs_dir){
    # patient clustering
    save_plt(plts, 'pat_clust-umap-cancers', '.pdf', figs_dir, width=9, height=9)
    save_plt(plts, 'pat_clust-umap-clusters', '.pdf', figs_dir, width=9, height=9)
    save_plt(plts, 'pat_clust-clusters_vs_cancers-balloon', '.pdf', figs_dir, width=9, height=8)
    save_plt(plts, 'pat_clust-clusters_vs_cancers-bar', '.pdf', figs_dir, width=5, height=12)
    
    # example of personalization
    save_plt(plts, 'pers-diffpsi-scatter', '.pdf', figs_dir, width=5.5, height=7)
    save_plt(plts, 'pers-harm-bars', '.pdf', figs_dir, width=10, height=7)
    save_plt(plts, 'pers-drug_rec-bars', '.pdf', figs_dir, width=9, height=7)
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
    selected_events = readLines(selected_events_file)
    selected_drugs = readLines(selected_drugs_file)
    spldep = read_tsv(spldep_file) %>% filter(index%in%selected_events)
    psi = read_tsv(psi_file) %>% filter(EVENT%in%selected_events)
    event_info = read_tsv(event_info_file)
    drug_response = read_tsv(drug_response_file) %>% filter(ID%in%selected_drugs)
    drug_targets = read_tsv(drug_targets_file)
    
    # prep
    drug_response = drug_response %>%
        mutate(DRUG_ID=as.numeric(gsub("_.*","",ID))) %>%
        left_join(drug_targets %>% distinct(DRUG_NAME,DRUG_ID), by="DRUG_ID") %>%
        mutate(DRUG_NAME_CLEAN = ID) %>%
        separate(DRUG_NAME_CLEAN, c("drug_id","min_conc","max_conc"), sep="_") %>%
        mutate(DRUG_NAME_CLEAN = sprintf("%s (%s-%s)", DRUG_NAME, min_conc, max_conc)) %>%
        left_join(metadata, by=c("sample"="sampleID"))
    
    harm = spldep %>%
        pivot_longer(cols=-index, names_to="sampleID", values_to="spldep") %>%
        left_join(
            psi %>% pivot_longer(cols=-EVENT, names_to="sampleID", values_to="psi"), 
            by=c("index"="EVENT","sampleID")) %>%
        mutate(harm_score = ifelse(
            spldep<0,
                (-1)*spldep*(0-psi), # harm if we remove exon completely
                (-1)*spldep*(100-psi))) %>% # harm if we insert exon completely
        left_join(metadata, by="sampleID") %>%
        drop_na()
    
    event_annot = event_info %>% 
        distinct(EVENT,GENE) %>% 
        mutate(event_gene=paste0(EVENT,"_",GENE))
    
    # plot
    plts = make_plots(harm, drug_response, event_annot)

    # make figdata
    figdata = make_figdata(embedding)

    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
