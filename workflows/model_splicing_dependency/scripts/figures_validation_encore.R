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
# - can we predict the proliferation effects of KDs from the changes in splicing?
# - if we define exon sets from KDs, is any KD enriched in our selected events?

require(optparse)
require(tidyverse)
require(ggpubr)
require(cowplot)
require(scattermore)
require(ggrepel)
require(extrafont)
require(tidytext)
require(clusterProfiler)
require(ComplexHeatmap)
require(ggplotify)
require(gridExtra)

ROOT = here::here()
source(file.path(ROOT,"src","R","utils.R"))

# vaiables 
# THRESH_DPSI = 5
RANDOM_SEED = 1234

# formatting
PAL_SINGLE_ACCENT = "orange"
PAL_SINGLE_LIGHT = "#47A59C"#"#6AC2BF"
PAL_SINGLE_DARK = "#C86A73"#"#716454"
PAL_DUAL = c(PAL_SINGLE_DARK, PAL_SINGLE_LIGHT)
PAL_SINGLE_NEUTRAL = "#716454"
PAL_FDR_DARK = "#005AB5"
PAL_FDR_LIGHT = "#DC3220"

LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"


# Development
# -----------
# RAW_DIR = file.path(ROOT,"data","raw")
# PREP_DIR = file.path(ROOT,"data","prep")
# RESULTS_DIR = file.path(ROOT,"results","model_splicing_dependency")
# rnai_file = file.path(PREP_DIR,"demeter2","CCLE.tsv.gz")
# diff_tpm_file = file.path(RESULTS_DIR,'files','ENCORE','diff_tpm.tsv.gz')
# delta_psi_file = file.path(RESULTS_DIR,'files','ENCORE','delta_psi-EX.tsv.gz')
# spldep_file = file.path(RESULTS_DIR,'files','ENCORE','splicing_dependency-EX','mean.tsv.gz')
# harm_score_file = file.path(RESULTS_DIR,"files","ENCORE","harm_score-EX.tsv.gz")
# selected_events_file = file.path(RESULTS_DIR,"files","selected_models-EX.txt")
# event_info_file = file.path(RAW_DIR,"VastDB","EVENT_INFO-hg38_noseqs.tsv")
# gene_annot_file = file.path(RAW_DIR,"HGNC","gene_annotations.tsv.gz")
# metadata_file = file.path(PREP_DIR,'metadata','ENCORE.tsv.gz')
# figs_dir = file.path(RESULTS_DIR,'figures','validation_encore')
# crispr_file = file.path(PREP_DIR,'Thomas2020','crispr_screen.tsv.gz')

##### FUNCTIONS #####
prepare_data = function(
    metadata, event_info, rnai, diff_tpm, delta_psi, 
    harm_score, selected_events, events_crispr, spldep
){
    
    cells_oi = metadata %>% pull(DepMap_ID) %>% unique()
    event_annot = event_info %>% distinct(GENE, EVENT)   
    metadata = metadata %>% 
        mutate(
            KD = ifelse(KD=="GNB2L1","RACK1",KD), # GNB2L1 is RACK1 in the latest annotation
            kd_lab = paste0(KD,"_",cell_line)
        )
    genes_kd = metadata %>% pull(KD) %>% unique() %>% na.omit()
    
    genedep = rnai %>% 
        filter(index %in% genes_kd) %>% 
        dplyr::select(c("index",all_of(cells_oi))) %>% 
        pivot_longer(-index, names_to="DepMap_ID", values_to="demeter2") %>% 
        rename(KD=index) %>%
        distinct() %>%
        drop_na()
    
    dtpm = metadata %>%
        distinct(sampleID, KD) %>%
        drop_na(KD) %>% # remove negative controls
        left_join(
            diff_tpm %>% 
                filter(NAME %in% genes_kd) %>% 
                pivot_longer(-c(ID,NAME), names_to="sampleID", values_to="fcTPM"),
            by = c("KD"="NAME", "sampleID")
        )    
    
    dpsi = delta_psi %>%
        filter(EVENT%in%selected_events) %>%
        pivot_longer(-EVENT, names_to="sampleID", values_to="deltaPSI") %>%
        #filter(abs(deltaPSI)>THRESH_DPSI) %>%
        rename(index=EVENT) %>%
        drop_na()
    
    spd = spldep %>%
        filter(index%in%selected_events) %>%
        pivot_longer(-index, names_to="sampleID", values_to="spldep") %>%
        drop_na()
    
    X = harm_score %>% 
        filter(index%in%selected_events) %>% 
        pivot_longer(-index, names_to="sampleID", values_to="harm") %>%
        filter(is.finite(harm) & harm>0) %>% # only positive harm score
        left_join(dpsi, by=c("index","sampleID")) %>%
        left_join(spd, by=c("index","sampleID")) %>%
        left_join(
            metadata %>% 
                mutate(lab=paste0(KD,"_",replicate)) %>% 
                distinct(sampleID,DepMap_ID,KD,cell_line,lab),
            by="sampleID"
        ) %>%
        left_join(genedep, by=c("KD","DepMap_ID")) %>% 
        left_join(dtpm, by=c("KD","sampleID")) %>%
        drop_na(deltaPSI, demeter2) %>%
        # summarize replicates
        group_by(cell_line, KD, demeter2, index) %>%
        summarize(deltaPSI = mean(deltaPSI),
                  spldep = mean(spldep),
                  harm = (-1) * mean(harm),
                  fcTPM = mean(fcTPM)) %>%
                  # harm = (-1) * sign(deltaPSI) * spldep) %>%
        ungroup() %>%
        group_by(cell_line, KD) %>%
        arrange(cell_line, KD, harm) %>%
        mutate(harm_rank = row_number()) %>%
        ungroup()
    
    return(X)
}

compute_correls_top_max = function(df){
    X = df
    
    set.seed(RANDOM_SEED)
    
    threshs_dtpm = c(0.5, 0, -0.5, -1, -1.5, -2)
    correls = lapply(threshs_dtpm, function(thresh_dtpm){
        correl = lapply(seq(1,100,2), function(x){
            corr_real = X %>% 
                filter(fcTPM < thresh_dtpm) %>% # ideally, consider only KDs comparable to Demeter2...
                filter(harm_rank <= x) %>% 
                group_by(cell_line, KD, demeter2) %>% 
                summarize(pred = sum(harm),
                          nobs = n()) %>% 
                ungroup() %>% 
                group_by(cell_line) %>% 
                summarize(correlation = cor(demeter2, pred, method="pearson"), 
                          pvalue = cor.test(demeter2, pred, method="pearson")[["p.value"]],
                          log10_pvalue = -log10(pvalue),
                          thresh = x,
                          dataset = "real")

            corr_rev = X %>% 
                filter(fcTPM < thresh_dtpm) %>%
                group_by(cell_line, KD, demeter2) %>% 
                slice_min(-harm, n=x) %>% # less harmful changes
                summarize(pred = sum(harm),
                          nobs = n()) %>% 
                ungroup() %>% 
                group_by(cell_line) %>% 
                summarize(correlation = cor(demeter2, pred, method="pearson"), 
                          pvalue = cor.test(demeter2, pred, method="pearson")[["p.value"]],
                          log10_pvalue = -log10(pvalue),
                          thresh = x,
                          dataset = "reversed")

            corr_null = X %>% 
                filter(fcTPM < thresh_dtpm) %>%
                group_by(cell_line, KD, demeter2) %>% 
                slice_sample(n=x, replace=TRUE) %>%
                distinct() %>%
                slice_min(harm, n=x) %>% 
                summarize(pred = sum(harm),
                          nobs = n()) %>% 
                ungroup() %>% 
                group_by(cell_line) %>% 
                summarize(correlation = cor(demeter2, pred, method="pearson"), 
                          pvalue = cor.test(demeter2, pred, method="pearson")[["p.value"]],
                          log10_pvalue = -log10(pvalue),
                          thresh = x,
                          dataset = "random")

            corr = rbind(corr_real, corr_rev, corr_null)

            return(corr)
        })
        correl = do.call(rbind, correl)
        correl = correl %>%
            mutate(
                thresh_dtpm = thresh_dtpm,
                thresh_dtpm_lab = sprintf("fcTPM<%s",thresh_dtpm)
            )
        return(correl)
    })
    correls = do.call(rbind, correls)
    
    return(correls)
}

compute_correls_dyn_ranges = function(df){
    X = df
    correls = lapply(seq(0,150,10), function(range_min){
        corr_higher = X %>% 
            filter(harm_rank > range_min) %>% 
            group_by(cell_line, KD, demeter2) %>% 
            summarize(pred = sum(harm),
                      nobs = n()) %>% 
            ungroup() %>% 
            group_by(cell_line) %>% 
            summarize(correlation = cor(demeter2, pred, method="pearson"), 
                      pvalue=cor.test(demeter2, pred, method="pearson")[["p.value"]],
                      log10_pvalue = -log10(pvalue),
                      dyn_range = paste0(">",as.character(range_min)),
                      dataset = "higher")
        
        w_size = 10
        range_name = sprintf("(%s,%s]", 
                             as.character(range_min), as.character(range_min+w_size))
        corr_windows = X %>% 
            filter(harm_rank > range_min & harm_rank <= (range_min + w_size)) %>% 
            group_by(cell_line, KD, demeter2) %>% 
            summarize(pred = sum(harm),
                      nobs = n()) %>% 
            ungroup() %>% 
            group_by(cell_line) %>% 
            summarize(correlation = cor(demeter2, pred, method="pearson"), 
                      pvalue=cor.test(demeter2, pred, method="pearson")[["p.value"]],
                      log10_pvalue = -log10(pvalue),
                      dyn_range = range_name,
                      dataset = "windows")
        
        corr = rbind(corr_higher, corr_windows)
        return(corr)
    })
    correls = do.call(rbind, correls)
}

plot_encore_validation = function(
    df, correls_top_max, correls_dyn_ranges, 
    event_info, metadata, events_crispr
){
    plts = list()
    
    X = df
    
    plts[["encore_val-diff_genexpr_kd-violin"]] = X %>%
        distinct(cell_line, fcTPM) %>%
        ggviolin(x="cell_line", y="fcTPM", trim=TRUE,
                 fill="cell_line", color=NA, palette=PAL_DUAL) + 
        geom_boxplot(width=0.1, outlier.size=0.1) + 
        guides(fill="none") + 
        labs(x="Cell Line", y="log2FC Gene Expr. KDs")
    
    # how many selected events change in each KD?
    plts[["encore_val-n_selected_events-violin"]] = X %>% 
        count(cell_line, KD) %>% 
        ggviolin(x="cell_line", y="n", trim=TRUE,
                 fill="cell_line", color=NA, palette=PAL_DUAL) + 
        geom_boxplot(width=0.1, outlier.size=0.1) + 
        guides(fill="none") + 
        labs(x="Cell Line", y="N. Cancer-Driver Exons in KD")
    
    # distribution of Demeter2 scores in RBPs
    plts[["encore_val-demeter2-violin"]] = X %>% 
        distinct(cell_line, KD, demeter2) %>% 
        ggviolin(x="cell_line", y="demeter2", fill="cell_line", 
                 color=NA, palette=PAL_DUAL, trim=TRUE) + 
        geom_boxplot(width=0.1, outlier.size=0.1) + 
        guides(fill="none") + 
        labs(x="Cell Line", y="Gene Dependency")  +
        geom_text(aes(y=0.5, label=lab),
                  X %>% 
                      distinct(cell_line, KD) %>% 
                      count(cell_line) %>% 
                      mutate(lab=paste0("n=",n)),
                  size=FONT_SIZE, family=FONT_FAMILY)
    
    # how does correlation change summing different top maximum?
    plts[["encore_val-thresh_vs_pearsons"]] = correls_top_max %>% 
        filter(thresh_dtpm == max(thresh_dtpm)) %>%
        ggscatter(x="thresh", y="correlation", palette=PAL_DUAL,
                  size="log10_pvalue", color="cell_line", alpha=0.5) + 
        labs(x="N Harmful Exons", y="Pearson Correlation", 
             color="Cell Line", size="-log10(p-value)") +
        scale_size(range=c(0.1,1.5)) +
        theme(aspect.ratio=1) + 
        facet_wrap(~dataset, ncol=1) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    plts[["encore_val-thresh_vs_pearsons_vs_fcTPM"]] = correls_top_max %>% 
        mutate(thresh_dtpm_lab = factor(
            thresh_dtpm_lab, levels=paste0("fcTPM<", sort(unique(thresh_dtpm))))
        ) %>%
        filter(thresh_dtpm %in% c(0.5, -2)) %>%
        ggscatter(x="thresh", y="correlation", palette=PAL_DUAL,
                  size="log10_pvalue", color="cell_line", alpha=0.5) + 
        labs(x="N Harmful Exons", y="Pearson Correlation", 
             color="Cell Line", size="-log10(p-value)") +
        scale_size(range=c(0.1,1.5)) +
        theme(aspect.ratio=1) + 
        facet_grid(dataset ~ thresh_dtpm_lab) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    # evaluate predictions with different dynamic ranges
    plts[["encore_val-harm-hist"]] = X %>% 
        gghistogram(x="harm_rank", fill="cell_line", palette=PAL_DUAL, color=NA) + 
        facet_wrap(~cell_line) + 
        labs(x="Harm Score Rank", y="Count") + 
        guides(fill="none") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    plts[["encore_val-ranges_vs_pearsons"]] = correls_dyn_ranges %>% 
        ggscatter(x="dyn_range", y="correlation", palette=PAL_DUAL,
                  size="log10_pvalue", color="cell_line", alpha=0.8) + 
        labs(x="Harm Score Rank Range", y="Pearson Correlation", 
             color="Cell Line", size="-log10(p-value)") +
        scale_size(range=c(0.1,1.5)) +
        theme(aspect.ratio=1) + 
        theme_pubr(x.text.angle = 70) +
        facet_wrap(~dataset, scales="free_x", ncol=1) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        ylim(0,0.4)
    
    # what influences harm score rankings? DeltaPSI or Spl. Dep.?
    x = X %>% 
        mutate(bin=cut(harm_rank, breaks=seq(0,100,10))) %>% 
        drop_na(bin)
    
    plts[["encore_val-harm_rank_vs_dpsi"]] = x %>%
        ggboxplot(x="bin", y="deltaPSI", outlier.size=0.1) +
        labs(x="Harm Score Rank", y="Delta PSI") +
        theme_pubr(x.text.angle=70) +
        facet_wrap(~cell_line, ncol=1) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    plts[["encore_val-harm_rank_vs_spldep"]] = x %>%
        ggboxplot(x="bin", y="spldep", outlier.size=0.1) +
        labs(x="Harm Score Rank", y="Splicing Dependency") +
        theme_pubr(x.text.angle=70) +
        facet_wrap(~cell_line, ncol=1) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    plts[["encore_val-dpsi_vs_spldep"]] = x %>%
        ggscatter(x="spldep", y="deltaPSI", size=0.1, alpha=0.8) +
        labs(x="Splicing Dependency", y="Delta PSI") +
        facet_wrap(~cell_line, ncol=1) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    # the most harmful changes predict overall knockdown
    x = X %>% 
        group_by(cell_line, KD, demeter2) %>% 
        slice_min(harm, n=1) %>% 
        summarize(pred=sum(harm)) %>%
        group_by(cell_line) %>%
        mutate(cell_line_lab=sprintf("%s (n=%s)", cell_line, n()))
    plts[["encore_val-top1-scatters"]] = x %>%
        ggscatter(x="pred", y="demeter2", color="cell_line", size=1, alpha=0.5,
                  palette=PAL_DUAL) + 
        guides(color="none") + 
        stat_cor(method="pearson", size=FONT_SIZE, family=FONT_FAMILY) + 
        geom_smooth(method="lm", linetype="dashed", color="black", size=LINE_SIZE) + 
        facet_wrap(~cell_line_lab, ncol=1, scales="free") + 
        geom_text_repel(aes(label=KD),
                        x %>% slice_max(demeter2*pred, n=5), segment.size=0.1,
                        size=FONT_SIZE, family=FONT_FAMILY) +
        labs(x="Sum of Top Harm Scores", y="Gene Dependency") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    x = X %>% 
        group_by(cell_line, KD, demeter2) %>% 
        slice_min(harm, n=10) %>% 
        summarize(pred=sum(harm)) %>%
        group_by(cell_line) %>%
        mutate(cell_line_lab=sprintf("%s (n=%s)", cell_line, n()))
    plts[["encore_val-top10-scatters"]] = x %>%
        ggscatter(x="pred", y="demeter2", color="cell_line", size=1, alpha=0.5,
                  palette=PAL_DUAL) + 
        guides(color="none") + 
        stat_cor(method="pearson", size=FONT_SIZE, family=FONT_FAMILY) + 
        geom_smooth(method="lm", linetype="dashed", color="black", size=LINE_SIZE) + 
        facet_wrap(~cell_line_lab, ncol=1, scales="free") + 
        geom_text_repel(aes(label=KD),
                        x %>% slice_max(demeter2*pred, n=5), segment.size=0.1,
                        size=FONT_SIZE, family=FONT_FAMILY) +
        labs(x="Sum of Top Harm Scores", y="Gene Dependency") +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    # who are these top 10 harmful events?
    x = X %>% 
        group_by(cell_line, KD, demeter2) %>% 
        slice_min(harm, n=10) %>%
        ungroup() %>% 
        count(cell_line, index) %>%
        group_by(cell_line) %>%
        slice_max(n, n=10) %>%
        left_join(event_info, by=c("index"="EVENT")) %>%
        mutate(event_gene=paste0(index,"_",GENE))
    plts[["encore_val-top10-bars"]] = x %>% 
        ggbarplot(x="event_gene", y="n", fill="cell_line", 
                  color=NA, palette=PAL_DUAL, position=position_dodge(0.9)) + 
        labs(x="Event & Gene", y="N. Exon in Top 10 Harm", fill="Cell Line") + 
        coord_flip()
    
    # harm scores of Thomas 2020 exons?
    plts[["encore_val-harm-thomas2020"]] = X %>% 
        filter(index %in% events_crispr) %>%
        left_join(event_info, by=c("index"="EVENT")) %>%
        mutate(event_gene=paste0(index,"_",GENE)) %>%
        ggboxplot(x="event_gene", y="harm_rank", fill="cell_line", 
                  palette=PAL_DUAL, outlier.size=0.1) + 
        labs(x="Event & Gene", y="Harm Score Rank", fill="Cell Line") +    
        coord_flip()
    
    # correlation between cell lines gene dependencies
    plts[["encore_val-demeter2-scatter"]] = X %>%
        distinct(cell_line, KD, demeter2) %>% 
        pivot_wider(id_cols="KD", 
                    names_from="cell_line", 
                    values_from="demeter2") %>% 
        ggscatter(x="K562", y="HepG2", size=1, alpha=0.5, color=PAL_SINGLE_NEUTRAL) + 
        stat_cor(method="pearson", size=FONT_SIZE, family=FONT_FAMILY) + 
        geom_smooth(method="lm", linetype="dashed", color="black", size=LINE_SIZE)
    
    return(plts)
}


make_plots = function(
    df, correls_top_max, correls_dyn_ranges, 
    event_info, metadata, events_crispr
){
    plts = list(
        plot_encore_validation(
            df, correls_top_max, correls_dyn_ranges, 
            event_info, metadata, events_crispr
        )
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(
    df, correls_top_max, correls_dyn_ranges
){
    
    evaluation = correls_top_max %>% 
        filter(thresh_dtpm == max(thresh_dtpm)) %>%
        dplyr::select(-c(thresh_dtpm, thresh_dtpm_lab)) %>%
        bind_rows(correls_dyn_ranges)
    
    figdata = list(
        "model_validation" = list(
            "spotter_results" = df,
            "evaluation" = evaluation
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
                    font.tickslab=6, font.family=FONT_FAMILY, device=cairo_pdf)    
    }
    filename = file.path(directory,paste0(plt_name,extension))
    save_plot(filename, plt, base_width=width, base_height=height, dpi=dpi, units="cm")
}


save_plots = function(plts, figs_dir){
    # correlations of predictions
    save_plt(plts, "encore_val-diff_genexpr_kd-violin", ".pdf", figs_dir, width=5, height=5)
    save_plt(plts, "encore_val-thresh_vs_pearsons", ".pdf", figs_dir, width=6, height=12)
    save_plt(plts, "encore_val-thresh_vs_pearsons_vs_fcTPM", ".pdf", figs_dir, width=12, height=12)
    save_plt(plts, "encore_val-n_selected_events-violin", ".pdf", figs_dir, width=5, height=5)
    save_plt(plts, "encore_val-demeter2-violin", ".pdf", figs_dir, width=5, height=5)
    save_plt(plts, "encore_val-top1-scatters", ".pdf", figs_dir, width=5, height=10)
    save_plt(plts, "encore_val-top10-scatters", ".pdf", figs_dir, width=5, height=10)
    save_plt(plts, "encore_val-top10-bars", ".pdf", figs_dir, width=5, height=6.5)
    
    save_plt(plts, "encore_val-harm-hist", ".pdf", figs_dir, width=6, height=4)
    save_plt(plts, "encore_val-ranges_vs_pearsons", ".pdf", figs_dir, width=4.5, height=10.25)
    
    save_plt(plts, "encore_val-harm_rank_vs_dpsi", ".pdf", figs_dir, width=5, height=8)
    save_plt(plts, "encore_val-harm_rank_vs_spldep", ".pdf", figs_dir, width=5, height=8)
    save_plt(plts, "encore_val-dpsi_vs_spldep", ".pdf", figs_dir, width=5, height=8)
    
    # harm scores in Thomas 2020?
    save_plt(plts, "encore_val-harm-thomas2020", ".pdf", figs_dir, width=6, height=8)
    
    # controls
    save_plt(plts, "encore_val-demeter2-scatter", ".pdf", figs_dir, width=5, height=5)
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
        make_option("--rnai_file", type="character"),
        make_option("--diff_tpm_file", type="character"),
        make_option("--delta_psi_file", type="character"),
        make_option("--spldep_file", type="character"),
        make_option("--harm_score_file", type="character"),
        make_option("--selected_events_file", type="character"),
        make_option("--event_info_file", type="character"),
        make_option("--metadata_file", type="character"),
        make_option("--crispr_file", type="character"),
        make_option("--psi_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    rnai_file = args[["rnai_file"]]
    diff_tpm_file = args[["diff_tpm_file"]]
    delta_psi_file = args[["delta_psi_file"]]
    spldep_file = args[["spldep_file"]]
    harm_score_file = args[["harm_score_file"]]
    selected_events_file = args[["selected_events_file"]]
    event_info_file = args[["event_info_file"]]
    gene_annot_file = args[["gene_annot_file"]]
    metadata_file = args[["metadata_file"]]
    crispr_file = args[["crispr_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    metadata = read_tsv(metadata_file)
    spldep = read_tsv(spldep_file)
    diff_tpm = read_tsv(diff_tpm_file)
    delta_psi = read_tsv(delta_psi_file)
    rnai = read_tsv(rnai_file)
    harm_score = read_tsv(harm_score_file)
    selected_events = readLines(selected_events_file)
    event_info = read_tsv(event_info_file)
    crispr = read_tsv(crispr_file)
    
    events_crispr = crispr %>% pull(EVENT) %>% unique()
    
    df = prepare_data(
        metadata, event_info, rnai, diff_tpm, delta_psi, 
        harm_score, selected_events, 
        events_crispr, spldep
    )
    correls_top_max = compute_correls_top_max(df)
    correls_dyn_ranges = compute_correls_dyn_ranges(df)
    
    # make plots
    plts = make_plots(
            df, correls_top_max, correls_dyn_ranges, 
            event_info, metadata, events_crispr
    )

    # make figdata
    figdata = make_figdata(
        df, correls_top_max, correls_dyn_ranges
    )
    
    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}