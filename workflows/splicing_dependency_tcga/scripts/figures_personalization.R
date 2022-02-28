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
# Patient clustering
# - clustering splicing dependency profiles across all cancers
# - one-vs-all differential analysis for each of the new clusters, 
#   specific vulnerabilities
# - find a cluster with many different cancer types and find their common
#   vulnerability; example with a cancer type with few patients.

require(tidyverse)
require(ggpubr)
require(cowplot)
require(extrafont)
require(scattermore)

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

# Development
# -----------
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# RESULTS_DIR = file.path(ROOT,'results','splicing_dependency_tcga')
# embedding_file = file.path(RESULTS_DIR,'files','PANCAN','embedded_splicing_dependency_mean-EX.tsv.gz')
# metadata_file = file.path(PREP_DIR,'metadata','PANCAN.tsv.gz')
# figs_dir = file.path(RESULTS_DIR,'figures','personalization')
# pers_spldep_file = file.path(RESULTS_DIR,'files','Seryakov2021','splicing_dependency-EX','mean.tsv.gz')
# pers_diff_deltas_file = file.path(RESULTS_DIR,'files','Seryakov2021','onediff-EX','deltas.tsv.gz')
# pers_diff_pvalues_file = file.path(RESULTS_DIR,'files','Seryakov2021','onediff-EX','pvalues.tsv.gz')
# pers_psi_file = file.path(RAW_DIR,'articles','Seryakov2021','vast_out','PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz')
# selected_events_file = file.path(ROOT,'results','model_splicing_dependency','files','selected_models-EX.txt')
# event_info_file = file.path(RAW_DIR,'VastDB','EVENT_INFO-hg38_noseqs.tsv')
# drug_response_file = file.path(RESULTS_DIR,'files','Seryakov2021','estimated_drug_response_by_drug-EX.tsv.gz')
# drug_targets_file = file.path(PREP_DIR,'drug_screens','drug_targets.tsv.gz')

##### FUNCTIONS #####
plot_patient_clusters = function(embedding, metadata){
    X = embedding %>% 
        left_join(metadata %>% distinct(sampleID, cancer, sample_type, OS, OS.time), 
                  by=c("index"="sampleID")) %>%
        filter(sample_type=="Primary Tumor")
    n_cancers = length(unique(X[["cancer"]]))
    n_clusters = length(unique(X[["leiden_labels"]]))
    
    plts = list()
    plts[["pat_clust-umap-cancers"]] = X %>% 
        ggplot(aes(x=UMAP0, y=UMAP1)) + 
        geom_scattermore(aes(color=cancer), pointsize=5, 
                         pixels=c(1000,1000), alpha=0.5) + 
        color_palette(get_palette("Paired", n_cancers)) + 
        theme_pubr() + 
        guides(color = guide_legend(override.aes = list(alpha = 1))) + 
        labs(color="Cancer Type") + 
        theme(aspect.ratio=1)
    
    plts[["pat_clust-umap-clusters"]] = X %>% 
        ggplot(aes(x=UMAP0, y=UMAP1)) + 
        geom_scattermore(aes(color=leiden_labels), pointsize=5, 
                         pixels=c(1000,1000), alpha=0.5) + 
        color_palette(get_palette("default", n_clusters)) + 
        theme_pubr() + 
        guides(color = guide_legend(override.aes = list(alpha = 1))) + 
        labs(color="Cluster") + 
        theme(aspect.ratio=1)
    
    plts[["pat_clust-clusters_vs_cancers-balloon"]] = X %>% 
        count(leiden_labels, cancer) %>%
        group_by(leiden_labels) %>%
        mutate(perc = n / sum(n)) %>%
        ggballoonplot(x="cancer", y="leiden_labels", size="perc", 
                      fill="tomato3", color="white") +
        scale_size(range=c(1,3)) + 
        labs(x="Cancer Type", y="Cluster")
    
    plts[["pat_clust-clusters_vs_cancers-bar"]] = X %>% 
        count(leiden_labels, cancer) %>%
        ggbarplot(x="leiden_labels", y="n", fill="cancer", position=position_dodge(0.9),
                  color=NA, palette=get_palette("Paired", n_cancers)) +
        labs(x="Cluster", y="Count", fill="Cancer Type") +
        coord_flip()
    
        
    return(plts)
}


plot_personalization_example = function(diff_result, harm, 
                                        drug_response, event_annot){
    
    plts = list()
    
    # differential analysis
    plts[["pers-diffpsi-scatter"]] = diff_result %>% 
        left_join(event_annot, by="EVENT") %>%
        mutate(log10_pvalue=-log10(pvalue)) %>%
        ggscatter(x="event_gene", y="deltapsi", 
                  size="log10_pvalue", color="#9BC1BC") +
        geom_hline(yintercept=0, linetype="dashed", size=1) +
        labs(x="Event & Gene", y="Delta PSI") + 
        scale_size(range=c(0.5,3)) +
        coord_flip()
    
    # vulnerabilities
    plts[["pers-harm-bars"]] = harm %>% 
        slice_max(abs(harm_score), n=15) %>%
        left_join(event_annot, by=c("index"="EVENT")) %>% 
        arrange(-spldep) %>% 
        pivot_longer(c(spldep, psi, harm_score)) %>% 
        mutate(name=factor(name, levels=c("spldep","psi","harm_score"))) %>% 
        ggbarplot(x="event_gene", y="value", color=NA, fill="#9BC1BC")  + 
        facet_wrap(~name, scales="free_x") + 
        labs(x="Event & Gene") + 
        coord_flip() +
        theme(strip.text.x = element_text(size=6, family='Arial'))
    
    # drug recommendation
    X = drug_response %>%
        group_by(DRUG_NAME_CLEAN) %>%
        slice_min(predicted_ic50, n=1) %>%
        ungroup()
    
    # ifosfamide (not in GDSC), DOXORUBICIN, PAZOPANIB did not work for the patient
    # oncobox recommended 
    drugs_oi = c("PAZOPANIB","DOXORUBICIN")
    x = X %>%
        arrange(predicted_ic50) %>%
        head(10) %>%
        mutate(dataset="top10") %>%
        bind_rows(
            X %>% filter(DRUG_NAME%in%drugs_oi) %>% mutate(dataset="treated")
        ) %>%
        mutate(lab=ifelse(DRUG_NAME%in%drugs_oi, 
                          paste0("*",DRUG_NAME_CLEAN), DRUG_NAME_CLEAN))
    
    plts[["pers-drug_rec-bars"]] = x %>%
        ggbarplot(x="lab", y="predicted_ic50", 
                  fill="#9BC1BC", color="drug_screen") +
        facet_wrap(~dataset) + 
        theme(strip.text.x = element_text(size=6, family='Arial')) +
        color_palette(c("black","white")) +
        labs(x="Drug", y="Predicted IC50") +
        coord_flip()
    
    return(plts)
}


make_plots = function(embedding, metadata, diff_result, harm, 
                      drug_response, event_annot){
    plts = list(
        plot_patient_clusters(embedding, metadata),
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
                    font.tickslab=6, font.family='Arial')    
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
    embedding = read_tsv(embedding_file) %>%
        mutate(leiden_labels=as.factor(leiden_labels))
    selected_events = readLines(selected_events_file)
    pers_spldep = read_tsv(pers_spldep_file) %>% filter(index%in%selected_events)
    pers_deltas = read_tsv(pers_diff_deltas_file) %>% filter(EVENT%in%selected_events)
    pers_pvalues = read_tsv(pers_diff_pvalues_file) %>% filter(EVENT%in%selected_events)
    pers_psi = read_tsv(pers_psi_file) %>% filter(EVENT%in%selected_events)
    event_info = read_tsv(event_info_file)
    drug_response = read_tsv(drug_response_file)
    drug_targets = read_tsv(drug_targets_file)
    
    # prep
    drug_response = drug_response %>%
        mutate(DRUG_ID=as.numeric(gsub("_.*","",ID))) %>%
        left_join(drug_targets %>% distinct(DRUG_NAME,DRUG_ID), by="DRUG_ID") %>%
        mutate(DRUG_NAME_CLEAN = ID) %>%
        separate(DRUG_NAME_CLEAN, c("drug_id","min_conc","max_conc"), sep="_") %>%
        mutate(DRUG_NAME_CLEAN = sprintf("%s (%s-%s)", DRUG_NAME, min_conc, max_conc))
    
    harm = pers_spldep %>%
        rename(spldep = SRR13664572_1) %>%
        left_join(pers_psi %>% rename(psi = SRR13664572_1), by=c("index"="EVENT")) %>%
        mutate(harm_score = ifelse(spldep<0,
                                   (-1)*spldep*(0-psi), # harm if we remove exon completely
                                   (-1)*spldep*(100-psi))) # harm if we insert exon completely
    
    diff_result = pers_deltas %>%
        rename(deltapsi = SRR13664572_1) %>%
        left_join(pers_pvalues %>% rename(pvalue = SRR13664572_1), by="EVENT") %>%
        mutate(padj = p.adjust(pvalue, method="fdr"))
    
    event_annot = event_info %>% 
        distinct(EVENT,GENE) %>% 
        mutate(event_gene=paste0(EVENT,"_",GENE))
    
    # plot
    plts = make_plots(embedding, metadata, diff_result, harm, 
                      drug_response, event_annot)

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
