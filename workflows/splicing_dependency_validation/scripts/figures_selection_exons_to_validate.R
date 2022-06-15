# Script purpose
# --------------
# select 10 exons for validation in two cell lines:
# - 4 cell line specific responses
# - 4 response in both
# - (1 extra negative control with a random RNA sequence)
#
# Outline
# -------
# 1. They must be targetable exons from our TCGA analysis
# 2. How do they vary across cancer cell lines (median vs. standard deviation)
# 3. Discard ORF disruption, only plated cells, 10 exons in total
# 4. Add information on drugs
# 5. Check splicing correlation between exons

require(tidyverse)
require(ggpubr)
require(cowplot)
require(scattermore)
require(ggrepel)
require(extrafont)
require(ComplexHeatmap)
require(ggplotify)

ROOT = here::here()
source(file.path(ROOT,"src","R","utils.R"))

# variables
#MIN_N_OBS = 20
THRESH_LR_PVALUE = 0.025
THRESH_CORR = 0.2
SIZE_CTL = 100
THRESH_INDICES = 0.3 # correlation with sample indices
THRESH_FDR = 0.05
THRESH_MEDIAN_DIFF = 5
MIN_SAMPLES = 10

THRESH_GENEXPR = 2
THRESH_PSI = 25
TOP_N = 5

THRESH_DRUGS_FDR = 0.1
THRESH_DRUGS_NOBS = 20

# formatting
PAL_SINGLE_ACCENT = "orange"
PAL_SINGLE_LIGHT = "#efb300ff"#"#6AC2BF"
PAL_SINGLE_DARK = "#007e67ff"
PAL_SINGLE_NEUTRAL = "#716454"
PAL_DUAL = c(PAL_SINGLE_DARK, PAL_SINGLE_LIGHT)
PAL_FDR_DARK = "#005AB5"
PAL_FDR_LIGHT = "#DC3220"

LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

# Development
# -----------
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,'support')

# annotation_file = file.path(RAW_DIR,'VastDB','event_annotation-Hs2.tsv.gz')
# event_info_file = file.path(RAW_DIR,'VastDB','EVENT_INFO-hg38_noseqs.tsv')

# genexpr_file = file.path(PREP_DIR,"genexpr_tpm","CCLE.tsv.gz")
# splicing_file = file.path(PREP_DIR,"event_psi","CCLE-EX.tsv.gz")
# ccle_stats_file = file.path(PREP_DIR,"stats","CCLE.tsv.gz")
# ccle_metadata_file = file.path(PREP_DIR,"metadata","CCLE.tsv.gz")

# protein_impact_file = file.path(ROOT,"data","raw","VastDB","PROT_IMPACT-hg38-v3.tab.gz")
# cosmic_genes_file = file.path(ROOT,"data","raw","COSMIC","cancer_gene_census.tsv")
# encore_ontology_file = file.path(ROOT,"results","model_splicing_dependency",'files','ENCORE','kd_gene_sets-EX.tsv.gz')
# cancer_events_file = file.path(ROOT,"support","cancer_events.tsv")

# CCLE_DIR = file.path(ROOT,"results","model_splicing_dependency")
# spldep_file = file.path(CCLE_DIR,"files","splicing_dependency-EX","mean.tsv.gz")
# spldep_ccle_file = file.path(CCLE_DIR,'files','splicing_dependency-EX','mean.tsv.gz')
# selected_events_file = file.path(CCLE_DIR,'files','selected_models-EX.txt')

# DRUGS_DIR = file.path(ROOT,"results","splicing_dependency_drugs")
# drug_models_file = file.path(DRUGS_DIR,"files","model_summaries_drug_response-EX.tsv.gz")
# drug_targets_file = file.path(PREP_DIR,'drug_screens','drug_targets.tsv.gz') 

# TCGA_DIR = file.path(ROOT,'results','splicing_dependency_tcga')

# diff_result_sample_file = file.path(TCGA_DIR,'files','PANCAN','mannwhitneyu-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz')
# spldep_stats_file = file.path(TCGA_DIR,'files','PANCAN','summary_splicing_dependency-EX.tsv.gz')

# diff_result_subtypes_file = file.path(TCGA_DIR,'files','PANCAN_subtypes','mannwhitneyu-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz')
# spldep_stats_subtypes_file = file.path(TCGA_DIR,'files','PANCAN_subtypes','summary_splicing_dependency-EX.tsv.gz')

# available_cells_file = file.path(SUPPORT_DIR,"available_cancer_cells.tsv")

# RESULTS_DIR = file.path(ROOT,'results','splicing_dependency_validation')
# figs_dir = file.path(RESULTS_DIR,'figures','selection_exons_to_validate')


##### FUNCTIONS #####
from_matrix_to_edgelist = function(X){
    
    idx = upper.tri(X)
    
    df = data.frame(
        row = rownames(X)[row(X)[idx]],
        col = colnames(X)[col(X)[idx]],
        correlation = X[idx]
    )
    
    return(df)
}


prep_diff_result = function(diff_result, spldep_stats){
    diff_result = diff_result %>%
        rename_all(recode, index = "EVENT") %>%
        mutate(event_type = gsub('Hsa','',gsub("[^a-zA-Z]", "",EVENT))) %>%
        left_join(spldep_stats, by=c('cancer_type','EVENT')) %>%
        drop_na(event_gene) %>%
        group_by(cancer_type) %>%
        mutate(psi__padj = p.adjust(psi__pvalue, 'fdr'),
               psi__log10_padj = -log10(psi__padj)) %>%
        mutate(psi__is_significant = psi__padj<THRESH_FDR & 
                                     abs(psi__median_diff)>THRESH_MEDIAN_DIFF)
    return(diff_result)
}


compute_harm_score = function(spldep, splicing){
    # compute harm score
    ## H.S. = (-1) * SplDep * DeltaPSI
    ## where DeltaPSI is in the direction of inclusion if SplDep>0, or exclusion if SplDep<0
    common_samples = intersect(colnames(spldep), colnames(splicing))
    common_samples = setdiff(common_samples, "event_gene")
    common_events = intersect(spldep[["event_gene"]], splicing[["event_gene"]])
       
    spldep_mat = spldep %>% column_to_rownames("event_gene")
    spldep_mat = spldep_mat[common_events, common_samples]
    splicing_mat = splicing %>% column_to_rownames("event_gene")
    splicing_mat = splicing_mat[common_events, common_samples]
    
    psi_final = spldep_mat
    psi_final[psi_final<0] = 0 # remove oncoexons
    psi_final[psi_final>0] = 100 # include tumor-suppressor exons
    harm = (-1) * spldep_mat * (psi_final - splicing_mat)
    
    harm = harm %>% rownames_to_column("event_gene")
    
    return(harm)
}


get_stats = function(X, index_col){
    X = X %>% column_to_rownames(index_col)
    X_stats = data.frame(
            avg = apply(X,1,mean, na.rm=TRUE),
            med = apply(X,1,median, na.rm=TRUE),
            std = apply(X,1,sd, na.rm=TRUE),
            min = apply(X,1,min, na.rm=TRUE),
            max = apply(X,1,max, na.rm=TRUE),
            q75 = apply(X,1,quantile, na.rm=TRUE, probs=0.75),
            q25 = apply(X,1,quantile, na.rm=TRUE, probs=0.25),
            range = apply(X,1,quantile, na.rm=TRUE, probs=0.75) - 
                    apply(X,1,quantile, na.rm=TRUE, probs=0.25),
            n_missing = rowSums(is.na(X))
        ) %>% rownames_to_column(index_col)
    return(X_stats)
}


plot_eda_transcriptome = function(ccle_stats, genes_oi, events_oi){
    # medians vs. standard deviation across all cancer cell lines

    plts = list()
    
    # expression
    X = ccle_stats %>%
        distinct(GENE, gene_median, gene_std)
    
    plts[["eda_transcriptome-genexpr-scatter"]] = X %>% 
        ggplot(aes(x=gene_median, y=gene_std)) + 
        geom_scattermore(pixels=c(1000,1000), pointsize=2, alpha=0.5, color=PAL_SINGLE_ACCENT) + 
        theme_pubr() + 
        geom_text_repel(aes(label=GENE), 
                        X %>% filter(GENE %in% genes_oi), max.overlaps=50,
                        segment.size=0.1, size=FONT_SIZE, family=FONT_FAMILY) +
        labs(x="median(log2(TPM+1))", y="Standard Dev.") +
        theme(aspect.ratio=1)
    
    # splicing
    X = ccle_stats %>%
        mutate(event_gene = paste0(EVENT,"_",GENE)) %>%
        distinct(EVENT,event_gene, event_median, event_std)
    
    plts[["eda_transcriptome-splicing-scatter"]] = X %>% 
        ggplot(aes(x=event_median, y=event_std)) + 
        geom_scattermore(pixels=c(1000,1000), pointsize=2, alpha=0.5, color=PAL_SINGLE_ACCENT) + 
        theme_pubr() + 
        geom_text_repel(aes(label=event_gene), 
                        X %>% filter(event_gene %in% events_oi), max.overlaps=50,
                        segment.size=0.1, size=FONT_SIZE, family=FONT_FAMILY) +
        labs(x="median(PSI)", y="Standard Dev.") +
        theme(aspect.ratio=1)
    
    return(plts)
}


plot_selection_exons = function(ccle_harm_stats, ccle_harm, 
                                ccle_splicing, ccle_genexpr,
                                events_genes, events_oi, 
                                protein_impact, available_cells){
    plts = list()
    
    # summarize the selection of exons for our experiment
    # requirements
    # ------------
    # 1. [X] Exclusion is harmful
    # 2. [X] Gene expression: log2(TPM+1)>2.5
    # 3. [X] Splicing: PSI>25
    # 4. [X] Not ORF disruption
    # 5. [X] 10 exons in total
    
    # distributions of harm scores of events
    X = ccle_harm_stats %>%
        filter(event_gene %in% events_oi) %>%
        mutate(EVENT = gsub("_.*","", event_gene),
               GENE = gsub(".*_","", event_gene)) %>%
        left_join(protein_impact, by="EVENT") %>%
        filter(str_detect(impact_clean, "Alternative protein")) %>%
        filter(event_gene != "HsaEX1036699_SFPQ") # too long of an exon
    
    plts[["selection_exons-median_vs_std"]] = X %>% 
        ggplot(aes(x=med, y=std)) + 
        geom_scattermore(pixels=c(1000,1000), pointsize=3, alpha=0.5) + 
        theme_pubr() + 
        geom_text_repel(aes(label=event_gene, color=impact_clean), 
                        X, max.overlaps=50,
                        segment.size=0.1, size=FONT_SIZE, family=FONT_FAMILY) +
        labs(x="log10(median(Harm Score))", y="Standard Deviation", 
             color="Splicing Impact",
             title=sprintf("n=%s", nrow(X))) +
        color_palette("jco") +
        theme(aspect.ratio=1)
    
    # get selected events
    top_selection = X %>% distinct()
    top_genes = top_selection %>% pull(GENE) %>% unique()
    top_events = top_selection %>% pull(event_gene) %>% unique()
    
    # correlations between splicing events across all cancer cell lines (in general)
    mat = ccle_splicing %>% 
        filter(event_gene %in% top_events) %>% 
        column_to_rownames("event_gene") %>% 
        as.matrix() %>% 
        t() %>% 
        cor(method="spearman", use="pairwise.complete.obs") %>% 
        round(2)
    
    plts[["selection_exons-splicing-corr"]] = mat %>% # standardize for visualization
        Heatmap(
            name="Spearman Corr.",
            row_names_gp = gpar(fontsize=6, fontfamily=FONT_FAMILY),
            column_names_gp = gpar(fontsize=6, fontfamily=FONT_FAMILY),
            heatmap_legend_param = list(legend_gp = gpar(fontsize=6, fontfamily=FONT_FAMILY)),
            cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.text(sprintf("%.2f", mat[i, j]), x, y, 
                              gp = gpar(fontsize=6, fontfamily=FONT_FAMILY))
        }) %>% 
        draw() %>%
        grid.grabExpr() %>%
        as.ggplot()
    
    # filter available cancer cell lines
    ccle_genexpr = ccle_genexpr %>%
        dplyr::select(any_of(c(available_cells,"GENE")))
    ccle_splicing = ccle_splicing %>%
        dplyr::select(any_of(c(available_cells,"event_gene")))
    ccle_spldep = ccle_spldep %>%
        dplyr::select(any_of(c(available_cells,"event_gene")))
    
    # select cell pairs
    compliant_samples_genexpr = ccle_genexpr %>% 
        filter(GENE %in% top_genes) %>% 
        pivot_longer(-GENE, names_to="sampleID", values_to="expression") %>%
        group_by(sampleID) %>% 
        filter(all(median(expression, na.rm=TRUE) > 0)) %>% 
        ungroup() %>%
        pull(sampleID) %>%
        unique()
    
    compliant_samples_psi = ccle_splicing %>% 
        filter(event_gene %in% top_events) %>% 
        pivot_longer(-event_gene, names_to="sampleID", values_to="psi") %>%
        group_by(sampleID) %>% 
        filter(all(median(psi, na.rm=TRUE) > 0)) %>% 
        ungroup() %>%
        pull(sampleID) %>%
        unique()
    
    compliant_samples_harm = ccle_harm %>%
        filter(event_gene %in% top_events) %>%
        pivot_longer(-event_gene, names_to="sampleID", values_to="harm") %>%
        group_by(sampleID) %>%
        filter(sum(!is.na(harm)) > 2) %>% # at least in two samples
        ungroup() %>%
        pull(sampleID) %>%
        unique()
    
    compliant_samples = intersect(compliant_samples_genexpr, compliant_samples_psi)
    compliant_samples = intersect(compliant_samples, compliant_samples_harm)
    
    ## split the list in two (high-variant and low-variant)
    high_var = top_selection %>% slice_max(std, n=5) %>% pull(event_gene) # find most different
    low_var = top_selection %>% slice_min(std, n=8) %>% pull(event_gene) # find most similar
    
    ## find pair of cells with opposite high-variant and 
    ## similar low-variant harm scores
    ### do not compute correlations if we are missing more than 1 harm score
    mat = ccle_harm %>% column_to_rownames("event_gene")

    corr_high = mat[high_var,compliant_samples] %>% 
        cor(method="pearson", use="pairwise.complete.obs") %>% 
        from_matrix_to_edgelist()
    corr_low = mat[low_var,compliant_samples] %>% 
        cor(method="pearson", use="pairwise.complete.obs") %>% 
        from_matrix_to_edgelist()
    corrs = corr_high %>%
        left_join(corr_low, by=c("row","col"), suffix=c("_high","_low"))
    
    # correlations of high variant vs low variant
    plts[["selection_exons-corrs-scatter"]] = corrs %>%
        ggplot(aes(x=correlation_high, y=correlation_low)) +
        geom_scattermore(pixels=c(1000,1000), pointsize=3, alpha=0.5) + 
        theme_pubr() +
        labs(x="Corr. Highly Variant", y="Corr. Lowly Variant") +
        theme(aspect.ratio=1)
    
    # harm scores in selected cells
    for (i in 1:TOP_N) {
        cells_oi = corrs %>% 
            mutate(score = (-1)*correlation_high + correlation_low) %>% 
            arrange(-score) %>%
            filter(row_number() == i) %>% # take top N cell-cell pair
            dplyr::select(row,col) %>% 
            unlist() %>%
            as.vector()
    
        X = ccle_harm %>%
            filter(event_gene %in% top_events) %>%
            dplyr::select(all_of(c("event_gene", cells_oi))) %>%
            group_by(event_gene) %>%
            mutate(diff = abs(get(cells_oi[1]) - get(cells_oi[2])),
                   avg = mean(c(get(cells_oi[1]), get(cells_oi[2])), na.rm=TRUE)) %>%
            ungroup() %>%
            arrange(-avg) %>%
            mutate(event_harm = ifelse(row_number()<=4, "harmless", "harmful")) %>%
            arrange(-diff) %>%
            mutate(event_var = ifelse(row_number()<=4, "var", "ctt"),
                   label = paste0(event_harm,"_",event_var)) %>%
            arrange(label,-diff) %>% 
            dplyr::select(c(event_gene,cells_oi,label)) %>%
            left_join(annot %>% distinct(EVENT, event_gene), by="event_gene") %>%
            left_join(event_info %>% distinct(EVENT, LE_o), by="EVENT") %>%
            left_join(protein_impact %>% distinct(EVENT, impact_clean), by="EVENT") %>%
            mutate(event_length = cut(LE_o, breaks=c(1,50,100,250,500,1000,2000,3000,4000)))

        #colors_events = c("harmful_var"="red", "harmful_ctt"="yellow", "harmless_ctt"="grey")
        colors_lengths = setNames(get_palette("Greens",3), X %>% pull(event_length) %>% unique() %>% sort())
        colors_impact = setNames(get_palette("jco",1), X %>% pull(impact_clean) %>% unique() %>% sort())
        colors_annot = list(
            #label = colors_events,
            event_length = colors_lengths,
            impact_clean = colors_impact
        )
        annotation_row = HeatmapAnnotation(df = X %>% 
                                               column_to_rownames("event_gene") %>% 
                                               dplyr::select(impact_clean, event_length), 
                                           name="Event Info", which="row", col=colors_annot)


        colors_cells = setNames(get_palette("Dark2",2), 
                                ccle_metadata %>% filter(DepMap_ID %in% cells_oi) %>% pull(CCLE_Name))
        colors_samples = list(CCLE_Name = colors_cells)
        annotation_col = HeatmapAnnotation(df = ccle_metadata %>%  
                                               filter(DepMap_ID %in% cells_oi) %>% 
                                               mutate(DepMap_ID = factor(DepMap_ID, levels = cells_oi)) %>%
                                               arrange(DepMap_ID) %>%
                                               column_to_rownames("DepMap_ID") %>%
                                               dplyr::select(CCLE_Name),
                                           name = "Sample Info", which="col", col=colors_samples)


        mat = X %>% 
            column_to_rownames("event_gene") %>%
            dplyr::select(all_of(cells_oi))
        
        print(mat)
        
        plts[[paste0("selection_exons-harm-heatmap-",i)]] = scale(mat) %>% # standardize for visualization
            Heatmap(
                cluster_rows = FALSE,
                name=paste("Harm Score",i), 
                right_annotation = annotation_row,
                top_annotation = annotation_col,
                row_names_gp = gpar(fontsize=6, fontfamily=FONT_FAMILY),
                column_names_gp = gpar(fontsize=6, fontfamily=FONT_FAMILY),
                heatmap_legend_param = list(legend_gp = gpar(fontsize=6, fontfamily=FONT_FAMILY)),
                cell_fun = function(j, i, x, y, width, height, fill) {
                        grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize=6, fontfamily=FONT_FAMILY))
            }) %>% 
            draw() %>%
            grid.grabExpr() %>%
            as.ggplot()

        # where are our exons with respect to all in those cell lines?
        events_sel = top_selection %>% pull(event_gene)
        genes_sel = top_selection %>% pull(GENE)

        # gene expression of selected exons and cells
        X = ccle_genexpr %>%
            dplyr::select(c(GENE,cells_oi)) %>%
            pivot_longer(-GENE, names_to="sampleID", values_to="expression")

        plts[[paste0("selection_exons-genexpr-violin-",i)]] = X %>%
            ggviolin(x="sampleID", y="expression", trim=TRUE,
                     color=NA, fill=PAL_SINGLE_ACCENT) +
            geom_boxplot(width=0.1, outlier.size=0.1) +
            geom_text_repel(aes(label=GENE),
                            X %>% filter(GENE %in% genes_sel), max.overlaps=50,
                            segment.size=0.1, size=FONT_SIZE, family=FONT_FAMILY) +
            labs(x="Cell Line", y="log2(TPM + 1)")

        # splicing of selected exons and cells
        X = ccle_splicing %>%
            dplyr::select(c(event_gene,cells_oi)) %>%
            pivot_longer(-event_gene, names_to="sampleID", values_to="psi")

        plts[[paste0("selection_exons-splicing-violin-",i)]] = X %>%
            ggviolin(x="sampleID", y="psi", trim=TRUE,
                     color=NA, fill=PAL_SINGLE_ACCENT) +
            geom_boxplot(width=0.1, outlier.size=0.1) +
            geom_text_repel(aes(label=event_gene),
                            X %>% filter(event_gene %in% events_sel), max.overlaps=50,
                            segment.size=0.1, size=FONT_SIZE, family=FONT_FAMILY) +
            labs(x="Cell Line", y="PSI")
    }
    
    return(plts)
}


make_plots = function(ccle_stats, genes_oi, events_oi,
                      ccle_harm_stats, ccle_harm, ccle_splicing, ccle_genexpr,
                      events_genes, protein_impact, available_cells){
    plts = list(
        plot_eda_transcriptome(ccle_stats, genes_oi, events_oi),
        plot_selection_exons(ccle_harm_stats, ccle_harm, ccle_splicing, ccle_genexpr,
                             events_genes, events_oi, protein_impact, available_cells)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(events_oi, annot, event_info, protein_impact, ccle_harm_stats,
                        cancer_events, cosmic_genes, encore_ontology){
    
    # a table of exons to be considered for validation
    # as they have therapeutic potential
    exons_info = annot %>% 
        filter(event_gene %in% events_oi) %>% 
        left_join(event_info %>% dplyr::select(-GENE), by="EVENT") %>% 
        left_join(protein_impact, by="EVENT") %>%
        left_join(
            ccle_harm_stats %>% 
            rename_at(vars(-event_gene),function(x) paste0("harm_",x)),
            by = "event_gene") %>%
        mutate(in_true_positive_set = EVENT %in% cancer_events[["EVENT"]],
               in_cosmic_genes = GENE %in% cosmic_genes[["gene"]])
    
    ## add ontology        
    exons_info = exons_info %>%
        left_join(
            encore_ontology %>% 
            filter(EVENT %in% exons_info[["EVENT"]]) %>%
            group_by(EVENT) %>%
            arrange(term) %>%
            summarize(putative_upstream_rbps = str_c(term, collapse = "; ")),
            by = "EVENT"
        )
            
    ## add drug associations
    exons_info = exons_info %>%
        left_join(
            drug_models %>% 
            left_join(
                drug_targets %>% 
                distinct(DRUG_ID, DRUG_NAME), 
                by="DRUG_ID"
            ) %>% 
            distinct(EVENT, DRUG_NAME) %>%
            group_by(EVENT) %>%
            summarize(drug_associations = str_c(DRUG_NAME, collapse = "; ")),
            by="EVENT"
        )
            
    figdata = list(
        "selection_exons" = list(
            "exons_info" = exons_info
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
        plt = ggpar(plt, font.title=8, font.subtitle=8, font.caption=8, 
                    font.x=8, font.y=8, font.legend=6,
                    font.tickslab=6, font.family=FONT_FAMILY)   
    }
    filename = file.path(directory,paste0(plt_name,extension))
    save_plot(filename, plt, base_width=width, base_height=height, dpi=dpi, units='cm')
}


save_plots = function(plts, figs_dir){
    # expression and splicing of selected exons overview
    save_plt(plts, 'eda_transcriptome-genexpr-scatter', '.pdf', figs_dir, width=8, height=8)
    save_plt(plts, 'eda_transcriptome-splicing-scatter', '.pdf', figs_dir, width=8, height=8)
    
    # selection of exons and cells for validation
    save_plt(plts, 'selection_exons-median_vs_std', '.pdf', figs_dir, width=8, height=9)
    save_plt(plts, 'selection_exons-corrs-scatter', '.pdf', figs_dir, width=8, height=9)
    
    for (i in 1:TOP_N){
#         cm = 1/2.54
#         pdf(file.path(figs_dir,sprintf('selection_exons-harm-heatmap-%s.pdf',i)), width=16*cm, height=13*cm)
#         plts[[sprintf('selection_exons-harm-heatmap-%s',i)]] %>% draw(legend_labels_gp = gpar(fontsize=6)) # this does nothing
#         dev.off()
        save_plt(plts, paste0('selection_exons-harm-heatmap-',i), '.pdf', figs_dir, width=18, height=15, format=FALSE)
        save_plt(plts, paste0('selection_exons-genexpr-violin-',i), '.pdf', figs_dir, width=7, height=7)
        save_plt(plts, paste0('selection_exons-splicing-violin-',i), '.pdf', figs_dir, width=10, height=10)
    }
    
    save_plt(plts, "selection_exons-splicing-corr", '.pdf', figs_dir, width=15, height=13)
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
    diff_result_sample_file = args$diff_result_sample_file
    diff_result_response_file = args$diff_result_response_file
    selected_events_file = args$selected_events_file
    spldep_file = args$spldep_file
    psi_ccle_file = args$psi_ccle_file
    spldep_ccle_file = args$spldep_ccle_file
    figs_dir = args$figs_dir
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    protein_impact = read_tsv(protein_impact_file) %>%
        dplyr::rename(EVENT=EventID, impact=ONTO) %>%
        mutate(
            impact = gsub("ORF disruption upon sequence inclusion \\(Alt\\. Stop\\)",
                          "Alternative protein isoforms \\(Ref, Alt\\. Stop\\)", impact),
            impact_clean=gsub(" \\(.*","",impact),
            impact_clean=gsub("ORF disruption upon sequence exclusion",
                            "ORF disruption (exclusion)",impact_clean),
            impact_clean=gsub("ORF disruption upon sequence inclusion",
                            "ORF disruption (inclusion)",impact_clean),
            impact_clean=gsub("In the CDS, with uncertain impact",
                             "In the CDS (uncertain)",impact_clean)
        )
    
    encore_ontology = read_tsv(encore_ontology_file)
    cancer_events = read_tsv(cancer_events_file)
    cosmic_genes = read_tsv(cosmic_genes_file) %>%
            dplyr::select("Gene Symbol") %>%
            rename(gene = `Gene Symbol`) %>%
            mutate(term = "COSMIC_CENSUS") %>%
            dplyr::select(term,gene)
    
    event_info = read_tsv(event_info_file)
    annot = read_tsv(annotation_file) %>%
        mutate(event_gene = paste0(EVENT,"_",GENE))
    events_genes = annot %>%
        dplyr::select(EVENT,GENE,ENSEMBL,event_gene)
    
    selected_events = readLines(selected_events_file)
    selected_genes = annot %>% 
        filter(EVENT %in% selected_events) %>% 
        pull(GENE) %>% 
        unique()
    selected_event_genes = annot %>% 
        filter(EVENT %in% selected_events) %>% 
        pull(event_gene) %>% 
        unique()

    tcga_diff_result = read_tsv(diff_result_sample_file)
    tcga_spldep_stats = read_tsv(spldep_stats_file) %>% 
        filter(EVENT %in% selected_events) %>%
        left_join(events_genes, by='EVENT')
    
    tcga_diff_result_subtypes = read_tsv(diff_result_subtypes_file)
    tcga_spldep_stats_subtypes = read_tsv(spldep_stats_subtypes_file) %>% 
        filter(EVENT %in% selected_events) %>%
        left_join(events_genes, by='EVENT')
    
    ccle_metadata = read_tsv(ccle_metadata_file)
    ccle_stats = read_tsv(ccle_stats_file)
    ccle_splicing = read_tsv(splicing_file) %>%
        left_join(events_genes %>% distinct(EVENT,event_gene), by="EVENT") %>%
        dplyr::select(-EVENT)
    ccle_genexpr = read_tsv(genexpr_file) %>%
        mutate_if(is.numeric, function(x){ log2(x+1) }) %>%
        left_join(events_genes %>% distinct(GENE,ENSEMBL), by=c("ID"="ENSEMBL")) %>%
        dplyr::select(-ID)
    ccle_spldep = read_tsv(spldep_file) %>% 
        filter(index %in% selected_events) %>%
        left_join(events_genes %>% distinct(EVENT,event_gene), by=c("index"="EVENT")) %>%
        dplyr::select(-index)
    
    drug_models = read_tsv(drug_models_file) %>% 
        mutate(event_gene = paste0(EVENT,"_",GENE),
               event_type = gsub("Hsa","",gsub("[^a-zA-Z]", "",EVENT))) %>%
        filter(lr_padj < THRESH_DRUGS_FDR & n_obs > THRESH_DRUGS_NOBS)
    drug_targets = read_tsv(drug_targets_file) %>% 
        mutate(DRUG_NAME=toupper(DRUG_NAME)) %>%
        distinct(DRUG_ID,DRUG_NAME,TARGET,TARGET_PATHWAY)
    
    available_cells = read_tsv(available_cells_file) %>%
        drop_na() %>%
        pull(DepMap_ID) %>%
        unique()
    
    # subset cell types
#     not_oi = c("Unknown","Engineered","Non-Cancerous","Fibroblast")
#     cultures_oi = c("Adherent","not found",NA)
#     to_drop = ccle_metadata %>%
#         filter((primary_disease %in% not_oi) | !(culture_type %in% cultures_oi)) %>%
#         pull(DepMap_ID)
    
    # prep results differential analyses
    tcga_diff_result = prep_diff_result(tcga_diff_result, tcga_spldep_stats)
    tcga_diff_result_subtypes = prep_diff_result(tcga_diff_result_subtypes, tcga_spldep_stats_subtypes)
    
    # compute harm scores and their summary stats for available cell lines
    ccle_harm = compute_harm_score(ccle_spldep %>% 
                                       filter(event_gene %in% selected_event_genes) %>%
                                       dplyr::select(any_of(c(available_cells,"event_gene"))), 
                                   ccle_splicing %>% 
                                       filter(event_gene %in% selected_event_genes) %>%
                                       dplyr::select(any_of(c(available_cells,"event_gene"))))
    ccle_harm_stats = get_stats(ccle_harm, "event_gene")
    
    # events and genes differentially spliced and targetable
    ## by cancer types
    events_oi_types = tcga_diff_result %>%
        # differentially spliced
        filter(psi__is_significant) %>%
        # targetable (only exclusion of exon)
        filter((sign(psi__median_diff)>0 & sign(mean)<0)) %>%
        pull(event_gene) %>%
        unique()
    
    genes_oi_types = tcga_diff_result %>%
        # differentially spliced
        filter(psi__is_significant) %>%
        # targetable (only exclusion of exon)
        filter((sign(psi__median_diff)>0 & sign(mean)<0)) %>%
        pull(GENE) %>%
        unique()
    
    ## by cancer subtypes
    events_oi_subtypes = tcga_diff_result_subtypes %>%
        # differentially spliced
        filter(psi__is_significant) %>%
        # targetable (only exclusion of exon)
        filter((sign(psi__median_diff)>0 & sign(mean)<0)) %>%
        pull(event_gene) %>%
        unique()
    
    genes_oi_subtypes = tcga_diff_result_subtypes %>%
        # differentially spliced
        filter(psi__is_significant) %>%
        # targetable (only exclusion of exon)
        filter((sign(psi__median_diff)>0 & sign(mean)<0)) %>%
        pull(GENE) %>%
        unique()
    
    ## combine
    events_oi = union(events_oi_types, events_oi_subtypes)
    genes_oi = union(genes_oi_types, genes_oi_subtypes)
    
    # plot
    plts = make_plots(ccle_stats, genes_oi, events_oi,
                      ccle_harm_stats, ccle_harm, ccle_splicing, ccle_genexpr,
                      events_genes, protein_impact, available_cells)
    
    # make figdata
    figdata = make_figdata(events_oi, annot, event_info, protein_impact, ccle_harm_stats,
                           cancer_events, cosmic_genes, encore_ontology)

    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}