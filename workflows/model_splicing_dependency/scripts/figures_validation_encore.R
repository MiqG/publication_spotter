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
require(scattermore)
require(ggrepel)
require(extrafont)
require(tidytext)
require(clusterProfiler)

ROOT = here::here()
source(file.path(ROOT,"src","R","utils.R"))

# vaiables 
THRESH_DPSI = 10

# Development
# -----------
# RAW_DIR = file.path(ROOT,"data","raw")
# PREP_DIR = file.path(ROOT,"data","prep")
# RESULTS_DIR = file.path(ROOT,"results","model_splicing_dependency")
# rnai_file = file.path(PREP_DIR,"demeter2","CCLE.tsv.gz")
# delta_psi_file = file.path(RESULTS_DIR,'files','ENCORE','delta_psi-EX.tsv.gz')
# harm_score_file = file.path(RESULTS_DIR,"files","ENCORE","harm_score-EX.tsv.gz")
# selected_events_file = file.path(RESULTS_DIR,"files","selected_models-EX.txt")
# event_info_file = file.path(RAW_DIR,"VastDB","EVENT_INFO-hg38_noseqs.tsv")
# metadata_file = file.path(PREP_DIR,'metadata','ENCORE.tsv.gz')
# ontology_file = file.path(RESULTS_DIR,'files','ENCORE','kd_gene_sets-EX.tsv.gz')


##### FUNCTIONS #####
plot_encore_validation = function(){
    plts = list()
    
    cells_oi = metadata %>% pull(DepMap_ID) %>% unique()
    genes_kd = metadata %>% pull(KD) %>% unique() %>% na.omit()
    event_annot = event_info %>% distinct(GENE, EVENT)
    metadata = metadata %>% mutate(kd_lab=paste0(KD,"_",cell_line))
    
    genedep = rnai %>% 
        filter(index%in%genes_kd) %>% 
        dplyr::select(c("index",all_of(cells_oi))) %>% 
        pivot_longer(-index, names_to="DepMap_ID", values_to="demeter2") %>% 
        rename(KD=index) %>%
        distinct() %>%
        drop_na()
    
    dpsi = delta_psi %>%
        filter(EVENT%in%selected_events) %>%
        pivot_longer(-EVENT, names_to="sampleID", values_to="deltaPSI") %>%
        filter(abs(deltaPSI)>THRESH_DPSI) %>%
        rename(index=EVENT) %>%
        drop_na()
    
    X = harm_score %>% 
        filter(index%in%selected_events) %>% 
        pivot_longer(-index, names_to="sampleID", values_to="harm") %>%
        filter(is.finite(harm) & harm>0) %>% # only positive harm score
        left_join(dpsi, by=c("index","sampleID")) %>%
        left_join(
            metadata %>% mutate(lab=paste0(KD,"_",replicate)) %>% distinct(sampleID,DepMap_ID,KD,cell_line),
            by="sampleID"
        ) %>%
        left_join(genedep, by=c("KD","DepMap_ID")) %>% 
        drop_na(deltaPSI, demeter2) %>%
        group_by(cell_line, KD, demeter2, index) %>%
        summarize(deltaPSI = mean(deltaPSI),
                  harm = mean(harm)) %>%
        ungroup() %>%
        mutate(sign_harm = sign(deltaPSI) * harm)
    
    # how does correlation change summing different top maximum?
    correls = lapply(1:25, function(x){
        corr = X %>% 
            group_by(cell_line, KD, demeter2) %>% 
            slice_max(harm, n=x) %>% 
            summarize(pred = sum(sign_harm),
                      nobs = n()) %>% 
            ungroup() %>% 
            group_by(cell_line) %>% 
            summarize(correlation = cor(demeter2, pred, method="pearson"), 
                      pvalue=cor.test(demeter2, pred, method="pearson")[["p.value"]],
                      log10_pvalue = -log10(pvalue),
                      thresh = x)
        return(corr)
    })
    correls = do.call(rbind, correls)
    
    plts[["encore_val-thresh_vs_pearsons"]] = correls %>% 
        ggscatter(x="thresh", y="correlation", palette="Set2",
                  size="log10_pvalue", color="cell_line") + 
        labs(x="N Most Harmful Exons", y="Pearson Correlation", 
             color="Cell Line", size="-log10(p-value)")
    
    # the most harmful changes predict overall knockdown
    x = X %>% 
        group_by(cell_line, KD, demeter2) %>% 
        slice_max(harm, n=1) %>% 
        summarize(pred=sum(sign_harm)) %>%
        group_by(cell_line) %>%
        mutate(cell_line_lab=sprintf("%s (n=%s)", cell_line, n()))
    plts[["encore_val-top1-scatters"]] = x %>%
        ggscatter(x="pred", y="demeter2", color="KD",
                  palette=get_palette("simpsons",length(genes_kd))) + 
        guides(color="none") + stat_cor(method="pearson") + 
        geom_smooth(method="lm", linetype="dashed", color="black") + 
        facet_wrap(~cell_line_lab, ncol=1, scales="free") + 
        geom_text_repel(aes(label=KD),
                        x %>% slice_max(demeter2*pred, n=5)) +
        labs(x="Sum of Top Harm Scores", y="Gene Dependency")
    
    x = X %>% 
        group_by(cell_line, KD, demeter2) %>% 
        slice_max(harm, n=10) %>% 
        summarize(pred=sum(sign_harm)) %>%
        group_by(cell_line) %>%
        mutate(cell_line_lab=sprintf("%s (n=%s)", cell_line, n()))
    plts[["encore_val-top10-scatters"]] = x %>%
        ggscatter(x="pred", y="demeter2", color="KD",
                  palette=get_palette("simpsons",length(genes_kd))) + 
        guides(color="none") + stat_cor(method="pearson") + 
        geom_smooth(method="lm", linetype="dashed", color="black") + 
        facet_wrap(~cell_line_lab, ncol=1, scales="free") + 
        geom_text_repel(aes(label=KD),
                        x %>% slice_max(demeter2*pred, n=5)) +
        labs(x="Sum of Top Harm Scores", y="Gene Dependency")
    
    # correlation between cell lines gene dependencies
    plts[["event_val-demeter2-scatter"]] = genedep %>% 
        left_join(metadata %>% distinct(DepMap_ID,cell_line), by="DepMap_ID") %>% 
        pivot_wider(id_cols="KD", 
                    names_from="cell_line", 
                    values_from="demeter2") %>% 
        ggscatter(x="K562", y="HepG2") + 
        stat_cor(method="pearson") + 
        geom_smooth(method="lm", linetype="dashed", color="black")
    
    # are some events changing in KDs enriched in our selected events? No.
    enrichment = enricher(
        selected_events, TERM2GENE=ontology, 
        universe=event_annot %>%  filter(str_detect(EVENT,"EX")) %>% 
          pull(EVENT) %>% unique(), 
        maxGSSize=Inf)
    
    plts[["event_val-enrichment-KD"]] = dotplot(enrichment) + theme_pubr()
    
    return(plts)
}


make_plots = function(){
    plts = list(
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(models, rnai_stats, cancer_events, 
                        eval_pvalue, eval_corr, 
                        enrichment, indices, spldep_stats, 
                        gene_mut_freq, event_mut_freq,
                        rnai, spldep, splicing, genexpr){
    # prep enrichments
    df_enrichs = do.call(rbind,
        lapply(names(enrichment), function(e){
            res = as.data.frame(enrichment[[e]])
            res[["ontology"]] = e
            return(res)
        })
    )
    
    figdata = list(
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
                    font.x=8, font.y=8, font.legend=8,
                    font.tickslab=6, font.family="Arial")    
    }
    filename = file.path(directory,paste0(plt_name,extension))
    save_plot(filename, plt, base_width=width, base_height=height, dpi=dpi, units="cm")
}


save_plots = function(plts, figs_dir){
    # model selection
    ## controls
    save_plt(plts, "model_sel-deps_sorted_vs_std_ctl_neg", ".pdf", figs_dir, width=5, height=5)
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
    delta_psi = read_tsv(delta_psi_file)
    rnai = read_tsv(rnai_file)
    harm_score = read_tsv(harm_score_file)
    selected_events = readLines(selected_events_file)
    event_info = read_tsv(event_info_file)
    ontology = read_tsv(ontology_file)
    
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