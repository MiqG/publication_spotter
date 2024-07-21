require(optparse)
require(ensembldb)
require(EnsDb.Hsapiens.v86)
require(BSgenome.Hsapiens.UCSC.hg38)
require(Biostrings)
require(tidyverse)
require(ggpubr)
require(cowplot)
require(ggrepel)
require(extrafont)

# variables
SELECTED_CELL_LINES = c(
    "A549_LUNG",
    "MDAMB231_BREAST",
    "HT29_LARGE_INTESTINE"
)

SELECTED_EXONS = c(
    "HsaEX0034998_KRAS",
    "HsaEX0070392_VLDLR",
    "HsaEX0050345_PRPF18",
    "HsaEX0049558_PPP1R12A",
    "HsaEX0026116_FNBP1",
    "HsaEX0071941_YAP1",
    "HsaEX0052877_RCC1",
    "HsaEX0044398_NUP85"
)

# formatting
PAL_SINGLE_ACCENT = "orange"
PAL_SINGLE_LIGHT = "#efb300ff"#"#6AC2BF"
PAL_SINGLE_DARK = "#007e67ff"
PAL_SINGLE_NEUTRAL = "#716454"
PAL_DUAL = c(PAL_SINGLE_DARK, PAL_SINGLE_LIGHT)
PAL_FDR_DARK = "#005AB5"
PAL_FDR_LIGHT = "#DC3220"
PAL_CELLS = setNames(
    get_palette("Accent", length(SELECTED_CELL_LINES)),
    SELECTED_CELL_LINES
)
PAL_REPLICATES = setNames(get_palette("simpsons",3), 1:3)

LINE_SIZE = 0.25

FONT_SIZE = 2 # for additional labels
FONT_FAMILY = "Arial"

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# SUPPORT_DIR = file.path(ROOT,'support')
# RESULTS_DIR = file.path(ROOT,'results','experimental_validation')

# annotation_file = file.path(RAW_DIR,'VastDB','event_annotation-Hs2.tsv.gz')
# event_info_file = file.path(RAW_DIR,'VastDB','EVENT_INFO-hg38_noseqs.tsv')

# splicing_file = file.path(PREP_DIR,"event_psi","inhouse-EX.tsv.gz")
# genexpr_file = file.path(PREP_DIR,"genexpr_tpm","inhouse.tsv.gz")
# spldep_file = file.path(RESULTS_DIR,"files","splicing_dependency-EX","mean.tsv.gz")

# sso_optimization_file = file.path(RAW_DIR,"experiments","validation_therapeutic_potential","20220928-psi-aso-optimization.tsv")
# validation_splicing_file = file.path(RAW_DIR,"experiments","validation_therapeutic_potential","20230215-psi-aso.tsv")
# validation_od_file = file.path(RAW_DIR,"experiments","validation_therapeutic_potential","clonogenic_assay-od-merged.tsv")

# prolif_file = file.path(RAW_DIR,"experiments","validation_therapeutic_potential","resazurin_cell_proliferation-merged.tsv")

# exon_info_file = file.path(ROOT,"support","exon_mapping_validated_exons.tsv")
# domains_file = file.path(RESULTS_DIR,'files','proteoform_info',"domains.tsv.gz")

# sso_seqs_file = file.path(SUPPORT_DIR,"sso_sequences.tsv")

# figs_dir = file.path(RESULTS_DIR,'figures','validation')


##### FUNCTIONS #####
plot_validation = function(validation_clonogenic, validation_harm_scores){
    
    plts = list()
    
    # Proliferation changes upon treatment
    plts[["validation-od_raw"]] = validation_clonogenic %>%
        #mutate(event_gene = fct_reorder(event_gene, -od, mean)) %>%
        ggplot(aes(x=event_gene, y=od)) +
        geom_boxplot(width=0.5, outlier.shape=NA) +
        geom_jitter(aes(color=as.factor(replicate_biological)), size=0.5, width=0.1) +
        color_palette(PAL_REPLICATES) +
        labs(x="Condition", y="Cell Proliferation (OD570)", color="Replicate") +
        theme_pubr(x.text.angle = 70) +
        stat_compare_means(ref.group="CONTROL_NEG", method="t.test", label="p.format", 
                           size=FONT_SIZE, family=FONT_FAMILY) + 
        facet_wrap(~CCLE_Name, ncol=1) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    plts[["validation-od_norm"]] = validation_clonogenic %>%
        #mutate(event_gene = fct_reorder(event_gene, -od_norm, mean)) %>%
        ggplot(aes(x=event_gene, y=od_norm)) +
        geom_boxplot(width=0.5, outlier.shape=NA) +
        geom_jitter(aes(color=as.factor(replicate_biological)), size=0.5, width=0.1) +
        color_palette(PAL_REPLICATES) +
        labs(x="Condition", y="Cell Proliferation (Norm. OD570)", color="Replicate") +
        theme_pubr(x.text.angle = 70) +
        stat_compare_means(ref.group="CONTROL_NEG", method="t.test", label="p.format", 
                           size=FONT_SIZE, family=FONT_FAMILY) + 
        geom_hline(yintercept=1, linetype='dashed', size=LINE_SIZE) +
        facet_wrap(~CCLE_Name, ncol=1) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    plts[["validation-od_norm_averaged"]] = validation_clonogenic %>%
        distinct(CCLE_Name, replicate_biological, event_gene, od_norm_avg) %>%
        #mutate(event_gene = fct_reorder(event_gene, -od_norm_avg, mean)) %>%
        ggplot(aes(x=event_gene, y=od_norm_avg)) +
        geom_boxplot(width=0.5, outlier.shape=NA) +
        geom_jitter(aes(color=as.factor(replicate_biological)), size=0.5, width=0.1) +
        color_palette(PAL_REPLICATES) +
        labs(x="Condition", y="Cell Proliferation (Norm. OD570)", color="Replicate") +
        theme_pubr(x.text.angle = 70) +
        stat_compare_means(ref.group="CONTROL_NEG", method="t.test", label="p.format", 
                           size=FONT_SIZE, family=FONT_FAMILY) + 
        geom_hline(yintercept=1, linetype='dashed', size=LINE_SIZE) +
        facet_wrap(~CCLE_Name, ncol=1) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    # comparison of all conditions FCs across cell lines
    X = validation_clonogenic %>%
        distinct(CCLE_Name, replicate_biological, event_gene, od_norm_avg) %>%
        group_by(CCLE_Name, event_gene) %>%
        summarize(
            od_norm_avg_avg = mean(od_norm_avg)
        ) %>%
        ungroup()
    
    plts[["validation-od_norm_averaged-a549_vs_mdamb231"]] = X %>%
        pivot_wider(names_from=CCLE_Name, values_from=od_norm_avg_avg) %>%
        ggplot(aes(x=A549_LUNG, y=MDAMB231_BREAST)) +
        geom_hline(yintercept = 1, linetype="dashed", size=LINE_SIZE, color="black") +
        geom_vline(xintercept = 1, linetype="dashed", size=LINE_SIZE, color="black") +
        geom_smooth(method="lm", linetype="dashed", color="black", size=LINE_SIZE, alpha=0.2) +
        geom_point(aes(color=event_gene), size=1) +
        geom_text_repel(
            aes(label=event_gene, color=event_gene), 
            size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1,
            max.overlaps=50
        ) +
        stat_cor(method="pearson", size=FONT_SIZE, family=FONT_FAMILY) +
        color_palette(get_palette("uchicago",12)) +
        theme_pubr() +
        theme(aspect.ratio=1) +
        labs(color="Condition")
    
    plts[["validation-od_norm_averaged-a549_vs_ht29"]] = X %>%
        pivot_wider(names_from=CCLE_Name, values_from=od_norm_avg_avg) %>%
        ggplot(aes(x=A549_LUNG, y=HT29_LARGE_INTESTINE)) +
        geom_hline(yintercept = 1, linetype="dashed", size=LINE_SIZE, color="black") +
        geom_vline(xintercept = 1, linetype="dashed", size=LINE_SIZE, color="black") +
        geom_smooth(method="lm", linetype="dashed", color="black", size=LINE_SIZE, alpha=0.2) +
        geom_point(aes(color=event_gene), size=1) +
        geom_text_repel(
            aes(label=event_gene, color=event_gene), 
            size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1,
            max.overlaps=50
        ) +
        stat_cor(method="pearson", size=FONT_SIZE, family=FONT_FAMILY) +
        color_palette(get_palette("uchicago",12)) +
        theme_pubr() +
        theme(aspect.ratio=1) +
        labs(color="Condition")
    
    plts[["validation-od_norm_averaged-mdamb231_vs_ht29"]] = X %>%
        pivot_wider(names_from=CCLE_Name, values_from=od_norm_avg_avg) %>%
        ggplot(aes(x=MDAMB231_BREAST, y=HT29_LARGE_INTESTINE)) +
        geom_hline(yintercept = 1, linetype="dashed", size=LINE_SIZE, color="black") +
        geom_vline(xintercept = 1, linetype="dashed", size=LINE_SIZE, color="black") +
        geom_smooth(method="lm", linetype="dashed", color="black", size=LINE_SIZE, alpha=0.2) +
        geom_point(aes(color=event_gene), size=1) +
        geom_text_repel(
            aes(label=event_gene, color=event_gene), 
            size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1,
            max.overlaps=50
        ) +
        stat_cor(method="pearson", size=FONT_SIZE, family=FONT_FAMILY) +
        color_palette(get_palette("uchicago",12)) +
        theme_pubr() +
        theme(aspect.ratio=1) +
        labs(color="Condition")
    
    # PSI changes upon treatment
    X = validation_harm_scores %>%
        pivot_longer(
            c(psi_treated, psi_scrambled, psi_untreated), names_to="condition", values_to="psi"
        ) %>%
        mutate(
            condition_lab = case_when(
                condition=="psi_treated" ~ "+SSOs",
                condition=="psi_scrambled" ~ "CONTROL_NEG",
                condition=="psi_untreated" ~ "WATER"
            )
        ) %>%
        mutate(
            condition_lab = factor(
                condition_lab, levels=c("CONTROL_NEG", "WATER", "+SSOs")
            )
        ) %>%
        drop_na(psi) 
    
    plts[["validation-psi_effects"]] = X %>%
        group_by(CCLE_Name, event_gene, condition_lab) %>%
        mutate(avg_psi = mean(psi, na.rm=TRUE)) %>%
        ungroup() %>%
        ggplot(aes(x=event_gene)) +
        geom_bar(
            aes(y=avg_psi, fill=condition_lab), 
            . %>% distinct(CCLE_Name, event_gene, avg_psi, condition_lab),
            stat="identity", 
            color=NA, position=position_dodge(0.9)
        ) +
        geom_point(aes(y=psi, fill=condition_lab), size=0.5, color="black", 
                   position=position_jitterdodge(jitter.width=0.2, dodge.width=0.9)) +
        fill_palette("nejm") + 
        theme_pubr(x.text.angle=70) +
        labs(x="Event & Gene", y="PSI", fill="Condition") +
        facet_wrap(~CCLE_Name, ncol=1) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY))
    
    # PSI RT-PCR vs PSI RNAseq
    X = validation_harm_scores %>%
        distinct(CCLE_Name, EVENT, psi_untreated, psi_untreated_rnaseq)
    
    plts[["validation-psi_rtpcr_vs_psi_rnaseq"]] = X %>%
        ggplot(aes(x=psi_untreated, y=psi_untreated_rnaseq)) +
        geom_smooth(method="lm", linetype="dashed", color="black", size=LINE_SIZE, alpha=0.2) +
        geom_point(aes(color=CCLE_Name)) +
        color_palette(PAL_CELLS) +
        stat_cor(method="pearson", size=FONT_SIZE, family=FONT_FAMILY) +
        theme_pubr() +
        facet_wrap(~CCLE_Name, nrow=1) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        guides(color="none") +
        labs(x="PSI RT-PCR", y="PSI RNAseq")
    
    # spotter predictions
    od_ctl = validation_clonogenic %>%
        filter(event_gene=="CONTROL_NEG") %>%
        group_by(CCLE_Name, replicate_biological) %>%
        summarize(od_ctl = mean(od)) %>%
        ungroup()
    
    x = validation_clonogenic %>%
        left_join(od_ctl, by=c("CCLE_Name","replicate_biological")) %>%
        mutate(od_fc = log2(od / od_ctl)) %>%
        group_by(CCLE_Name, event_gene, replicate_biological) %>%
        summarize(
            od_fc = mean(od_fc),
            replicate_biological = as.factor(replicate_biological)
        ) %>%
        ungroup() %>%
        distinct() %>%
        left_join(
            validation_harm_scores %>%
                group_by(CCLE_Name, event_gene) %>%
                summarize(
                    harm_score = mean(harm_score, na.rm=TRUE),
                    delta_psi = mean(delta_psi, na.rm=TRUE)
                ) %>%
                ungroup(), 
            by=c("CCLE_Name","event_gene")
        ) %>%
        drop_na(harm_score, od_fc)
    
    plts[["validation-od_vs_harm"]] = x %>%
        ggplot(aes(x=harm_score, y=od_fc)) +
        geom_smooth(method="lm", linetype="dashed", color="black", size=LINE_SIZE) +
        geom_point(aes(color=event_gene)) +
        #geom_text_repel(aes(label=event_gene), size=FONT_SIZE, family=FONT_FAMILY, segment.size=0.1) +
        stat_cor(method="pearson", size=FONT_SIZE, family=FONT_FAMILY) +
        facet_wrap(~CCLE_Name+replicate_biological, ncol=2, scales="free_y") +
        color_palette(get_palette("simpsons",8)) +
        theme_pubr() +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY), aspect.ratio=1) +
        labs(x="Harm Score", y="Obs. Cell Prolif. Effect", color="Replicate")
    
    plts[["validation-od_vs_harm_combined"]] = x %>%
        ggplot(aes(x=harm_score, y=od_fc)) +
        geom_smooth(aes(color=replicate_biological, fill=replicate_biological), method="lm", linetype="dashed", size=LINE_SIZE, se=FALSE) +
        geom_point(aes(color=replicate_biological)) +
        stat_cor(
            aes(color=replicate_biological), 
            method="pearson", size=FONT_SIZE, family=FONT_FAMILY) +
        facet_wrap(~CCLE_Name, ncol=2, scales="free_y") +
        color_palette(PAL_REPLICATES) +
        fill_palette(PAL_REPLICATES) +
        theme_pubr() +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY), aspect.ratio=1) +
        labs(x="Harm Score", y="Obs. Cell Prolif. Effect", color="Event & Gene")
    
    plts[["validation-od_vs_harm_combined_minmax"]] = x %>%
        group_by(CCLE_Name, event_gene, harm_score) %>%
        summarize(od_fc = median(od_fc, na.rm=TRUE)) %>%
        ungroup() %>%
        group_by(CCLE_Name) %>%
        mutate(
            od_fc = (od_fc - min(od_fc))/(max(od_fc) - min(od_fc))
        ) %>%
        ungroup() %>%
        ggplot(aes(x=harm_score, y=od_fc)) +
        geom_smooth(
            method="lm", color="black", fill="lightgray", 
            linetype="dashed", size=LINE_SIZE, se=TRUE
        ) +
        geom_point(aes(color=CCLE_Name)) +
        color_palette(PAL_CELLS) +
        fill_palette(PAL_CELLS) +
        stat_cor(
            method="pearson", size=FONT_SIZE, family=FONT_FAMILY
        ) + 
        theme_pubr() +
        theme(aspect.ratio=1) +
        labs(x="Harm Score", y="Scaled Obs. Cell Prolif. Effect", color="Cell Line")
    
    plts[["validation-od_vs_harm_combined_centered"]] = x %>%
        group_by(CCLE_Name, event_gene, harm_score) %>%
        summarize(od_fc = median(od_fc, na.rm=TRUE)) %>%
        ungroup() %>%
        group_by(CCLE_Name) %>%
        mutate(
            od_fc = od_fc - median(od_fc, na.rm=TRUE)
        ) %>%
        ungroup() %>%
        ggplot(aes(x=harm_score, y=od_fc)) +
        geom_smooth(
            method="lm", color="black", fill="lightgray", 
            linetype="dashed", size=LINE_SIZE, se=TRUE
        ) +
        geom_point(aes(color=CCLE_Name)) +
        color_palette(PAL_CELLS) +
        fill_palette(PAL_CELLS) +
        stat_cor(
            method="pearson", size=FONT_SIZE, family=FONT_FAMILY
        ) + 
        theme_pubr() +
        theme(aspect.ratio=1) +
        labs(x="Harm Score", y="Norm. Obs. Cell Prolif. Effect", color="Cell Line")

    plts[["validation-od_vs_harm_combined_centered-high_prolif"]] = x %>%
        group_by(CCLE_Name, event_gene, harm_score) %>%
        summarize(od_fc = median(od_fc, na.rm=TRUE)) %>%
        ungroup() %>%
        group_by(CCLE_Name) %>%
        mutate(
            od_fc = od_fc - median(od_fc, na.rm=TRUE)
        ) %>%
        ungroup() %>%
        filter(CCLE_Name != "HT29_LARGE_INTESTINE") %>%
        ggplot(aes(x=harm_score, y=od_fc)) +
        geom_smooth(
            method="lm", color="black", fill="lightgray", 
            linetype="dashed", size=LINE_SIZE, se=TRUE
        ) +
        geom_point(aes(color=CCLE_Name)) +
        color_palette(PAL_CELLS) +
        fill_palette(PAL_CELLS) +
        stat_cor(
            method="pearson", size=FONT_SIZE, family=FONT_FAMILY
        ) + 
        theme_pubr() +
        theme(aspect.ratio=1) +
        labs(x="Harm Score", y="Norm. Obs. Cell Prolif. Effect", color="Cell Line")
    
    # SSO effect vs deltaPSI across cell types
    plts[["validation-od_fc_vs_delta_psi"]] = x %>% 
        ggscatter(x="delta_psi", y="od_fc", color="CCLE_Name", palette=PAL_CELLS, size=1) +
        geom_smooth(
            method="lm", color="black", fill="lightgray", 
            linetype="dashed", size=LINE_SIZE, se=TRUE
        ) +
        stat_cor(
            method="pearson", size=FONT_SIZE, family=FONT_FAMILY
        ) + 
        theme_pubr() +
        facet_wrap(~event_gene, scales="free_x", nrow=2) +
        theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="Delta PSI", y="Obs. Cell Prolif. Effect", color="Cell Line")
    
    return(plts)
}


plot_prolif = function(prolif){
    plts = list()
    
    X = prolif %>% 
        mutate(timepoint = factor(timepoint, levels=c(0,24,48,72,96)))
    
    plts[["prolif-time_vs_od-scatter"]] = X %>%
        ggplot() +
        geom_point(aes(x=timepoint, y=od_norm, color=CCLE_Name), size=0.5) +
        geom_smooth(
            aes(x=as.numeric(timepoint), y=od_norm, color=CCLE_Name, fill=CCLE_Name), 
            linetype="dashed", size=LINE_SIZE, alpha=0.2, span=0.5
        ) +
        theme_pubr() +
        theme(aspect.ratio = 1) + 
        color_palette(PAL_CELLS) +
        fill_palette(PAL_CELLS) +
        labs(x="Time (h)", y="Relative FLuorescence", color="Cell Line", fill="Cell Line")
    
    return(plts)
}


plot_protein_domains = function(exon_info){
    
    plts = list()
    X = exon_info
    
    palette = get_palette("rickandmorty", length(unique(X[["description"]])))
    plts[["protein_domains-sequence-rect"]] = X %>%
        ggplot() +
        # main proteoform
        geom_rect(aes(xmin=1, xmax=main_inclusion_proteoform_length, ymin=-0.25, ymax=0.25), fill="darkgrey") +
        # mapped exon
        geom_rect(aes(xmin=exon_proteoform_start, xmax=exon_proteoform_end, ymin=-0.25, ymax=0.25), fill="orange") +
        # mapped domains
        geom_rect(aes(xmin=start, xmax=end, ymin=(-0.45*index)-0.15, ymax=-0.45*index, fill=description)) +
        # add domain labels 
        geom_text(
            aes(x=main_inclusion_proteoform_length+25, y=(-0.45*index)-0.075, label=description),
            size=FONT_SIZE, family=FONT_FAMILY, hjust=0
        ) +
        fill_palette(palette) +
        facet_wrap(~proteoform_lab, ncol=1) +
        theme_pubr() +
        theme(
            strip.text.x = element_text(size=6, family=FONT_FAMILY),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank()
        ) +
        xlim(1,2000) +
        ylim(NA,2) +
        guides(fill="none") +
        labs(x="Protein Position", y="")
    
    return(plts)
}


plot_sso_seqs = function(sso_seqs){
    
    EDB = EnsDb.Hsapiens.v86
    dna = BSgenome.Hsapiens.UCSC.hg38
    
    df = sso_seqs %>%
        mutate(
            chromosome = as.integer(gsub("chr","",gsub(":.*","",COORD_o))),
            event_start = as.integer(gsub(".*:","",gsub("-.*","",COORD_o))),
            event_end = as.integer(gsub(".*-","",COORD_o)),
            event_strand = gsub(".*:","",REF_CO)
        ) %>% as.data.frame()
    
    plts_tmp = list()
    for (i in 1:nrow(df)) {
        # unpack
        gene_oi = df[i, "GENE"]
        event_gene = df[i, "event_gene"]
        aso_id = df[i, "aso_id"]
        sso_seq = df[i, "sequence"]
        event_start = df[i, "event_start"]
        event_end = df[i, "event_end"]
        event_strand = df[i, "event_strand"]
        margin = 25
        
        # get sequence
        gene_position = genes(EDB, filter = GeneNameFilter(gene_oi))
        gene_position = GRanges(
            seqnames = paste0("chr", as.character(seqnames(gene_position))), 
            ranges = ranges(gene_position)
        )
        gene_sequence = getSeq(dna, gene_position)[[1]]

        if (event_strand=="+"){
            sequence_to_find = reverseComplement(RNAString(sso_seq))
        }else{
            sequence_to_find = RNAString(sso_seq)
        }

        matches = matchPattern(pattern=sequence_to_find, subject=gene_sequence)
        target = subseq(gene_sequence, start=start(matches)-margin, end=end(matches)+margin)
        
        # plot
        X = data.frame(
                nt = target %>% as.character() %>% strsplit(split="") %>% unlist()
            ) %>%
            mutate(
                position = row_number(),
                genomic_position = position + start(matches) + start(gene_position) - margin,
                sso_target = (position >= margin) & (position <= (length(target)-margin)),
                feature = ifelse(
                    (genomic_position >= event_start) & (genomic_position <= event_end),
                    "exon",
                    "intron"
                )
            )
        
        title = sprintf("%s-%s", event_gene, aso_id)
        plts_tmp[[sprintf("sso_seqs-target_region-%s-%s", event_gene, aso_id)]] = X %>%
            ggplot(aes(x=genomic_position)) +
            geom_rect(aes(xmin=min(genomic_position)-0.5, xmax=max(genomic_position)+0.5, 
                          ymin=-1, ymax=1),
                      . %>% filter(sso_target), 
                      fill="lightgray", alpha=0.5) + 
            geom_text(aes(y=0, label=nt, color=feature), size=FONT_SIZE, family=FONT_FAMILY) +
            geom_text_repel(
                aes(y=0, label=genomic_position),
                . %>% filter(sso_target) %>% slice_min(genomic_position),
                min.segment.length=0, nudge_y=4, segment.size=0.1,
                size=FONT_SIZE, family=FONT_FAMILY
            ) +
            geom_text_repel(
                aes(y=0, label=genomic_position),
                . %>% filter(sso_target) %>% slice_max(genomic_position),
                min.segment.length=0, nudge_y=4, segment.size=0.1,
                size=FONT_SIZE, family=FONT_FAMILY
            ) +
            color_palette("Dark2") + 
            theme_pubr() +
            theme(axis.line.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
            ylim(-5,5) +
            labs(x="Genomic Position", y="", subtitle=title)        
    }
    
    plts = list(
        "sso_seqs-target_regions" = ggarrange(plotlist=plts_tmp, ncol=2, nrow=5, common.legend=TRUE)
    )
    
    return(plts)
}


plot_sso_optimization = function(sso_optimization){
    plts = list()
    
    X = sso_optimization 

    plts[["sso_optimization-ssos_vs_psi-bar"]] = X %>%
        ggbarplot(x="aso_id", y="psi", fill="condition", palette="nejm", color=NA, position=position_dodge(0.9)) +
        geom_text(
            aes(label="#", group=condition), 
            . %>% filter(selected_treatment),
            position=position_dodge(0.9)
        ) +
        facet_wrap(~event_gene, nrow=2) +
        theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
        labs(x="SSO ID", y="PSI", fill="Condition")
    
    return(plts)
}

make_plots = function(validation_clonogenic, validation_harm_scores, prolif, exon_info, sso_optimization, sso_seqs){
    plts = list(
        plot_validation(validation_clonogenic, validation_harm_scores),
        plot_prolif(prolif),
        plot_protein_domains(exon_info),
        plot_sso_seqs(sso_seqs),
        plot_sso_optimization(sso_optimization)
    )
    plts = do.call(c,plts)
    return(plts)
}


make_figdata = function(validation_clonogenic, validation_harm_scores, prolif, exon_info, sso_optimization, sso_seqs){

    figdata = list(
        "experiments" = list(
            "validation_clonogenic" = validation_clonogenic,
            "validation_harm_scores" = validation_harm_scores,
            "proliferation_assay" = prolif,
            "exon_info" = exon_info,
            "sso_optimization" = sso_optimization,
            "sso_seqs" = sso_seqs
        )
    )
}

make_source_data = function(plts){
    
    source_data = list(
        # FIGURE 3
        ## Fig. 3c
        "fig03c" = plts[["validation-od_norm_averaged"]][["data"]] %>%
            filter(CCLE_Name!="HT29_LARGE_INTESTINE"),
        
        # SUPPLEMENTARY FIGURE 8
        ## Sup. Fig. 8b
        "supfig08b" = plts[["validation-psi_effects"]][["data"]] %>%
            filter(CCLE_Name!="HT29_LARGE_INTESTINE"),
        
        ## Sup. Fig. 8c
        "supfig08c" = plts[["validation-od_vs_harm_combined_centered-high_prolif"]][["data"]],
        
        # SUPPLEMENTARY FIGURE 10
        ## Sup. Fig. 10a
        "supfig10a" = plts[["prolif-time_vs_od-scatter"]][["data"]],
        
        ## Sup. Fig. 10c
        "supfig10c" = plts[["validation-od_norm_averaged"]][["data"]] %>%
            filter(CCLE_Name=="HT29_LARGE_INTESTINE"),
        
        # SUPPLEMENTARY FIGURE 11
        ## Sup. Fig. 11b
        "supfig11b" = plts[["validation-psi_effects"]][["data"]] %>%
            filter(CCLE_Name=="HT29_LARGE_INTESTINE"),
        
        # SUPPLEMENTARY FIGURE 13
        "supfig13" = plts[["validation-od_norm_averaged-a549_vs_mdamb231"]][["data"]]
    )
    
    return(source_data)
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
    save_plt(plts, "validation-od_raw", '.pdf', figs_dir, width=5, height=16)
    save_plt(plts, "validation-od_norm", '.pdf', figs_dir, width=5, height=16)
    save_plt(plts, "validation-od_norm_averaged", '.pdf', figs_dir, width=5, height=16)
    
    save_plt(plts, "validation-od_norm_averaged-a549_vs_mdamb231", '.pdf', figs_dir, width=9, height=10)
    save_plt(plts, "validation-od_norm_averaged-a549_vs_ht29", '.pdf', figs_dir, width=9, height=10)
    save_plt(plts, "validation-od_norm_averaged-mdamb231_vs_ht29", '.pdf', figs_dir, width=9, height=10)
    
    save_plt(plts, "validation-psi_effects", '.pdf', figs_dir, width=4.5, height=16)
    save_plt(plts, "validation-psi_rtpcr_vs_psi_rnaseq", '.pdf', figs_dir, width=8, height=4.5)
    save_plt(plts, "validation-od_vs_harm", '.pdf', figs_dir, width=12, height=19)
    save_plt(plts, "validation-od_vs_harm_combined", '.pdf', figs_dir, width=8, height=9)
    save_plt(plts, "validation-od_vs_harm_combined_minmax", '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, "validation-od_vs_harm_combined_centered", '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, "validation-od_vs_harm_combined_centered-high_prolif", '.pdf', figs_dir, width=5, height=6)
    save_plt(plts, "validation-od_fc_vs_delta_psi", '.pdf', figs_dir, width=9, height=9)
    save_plt(plts, "prolif-time_vs_od-scatter", '.pdf', figs_dir, width=5, height=6)
    
    save_plt(plts, "protein_domains-sequence-rect", '.pdf', figs_dir, width=10, height=65)
    save_plt(plts, "sso_seqs-target_regions", '.pdf', figs_dir, width=25, height=20)
    
    save_plt(plts, "sso_optimization-ssos_vs_psi-bar", '.pdf', figs_dir, width=11.5, height=9)
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

save_source_data = function(source_data, dir){
    d = file.path(dir,"figdata",'source_data')
    dir.create(d, recursive=TRUE)
    lapply(names(source_data), function(nm){
        df = source_data[[nm]]
        filename = file.path(d, paste0(nm,'.tsv.gz'))
        write_tsv(df, filename)
        print(filename)
    })
}

parseargs = function(){
    
    option_list = list( 
        make_option("--annotation_file", type="character"),
        make_option("--event_info_file", type="character"),
        make_option("--validation_splicing_file", type="character"),
        make_option("--splicing_file", type="character"),
        make_option("--genexpr_file", type="character"),
        make_option("--spldep_file", type="character"),
        make_option("--validation_spldep_file", type="character"),
        make_option("--validation_od_file", type="character"),
        make_option("--prolif_file", type="character"),
        make_option("--exon_info_file", type="character"),
        make_option("--domains_file", type="character"),
        make_option("--sso_seqs_file", type="character"),
        make_option("--sso_optimization_file", type="character"),
        make_option("--figs_dir", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    annotation_file = args[["annotation_file"]]
    event_info_file = args[["event_info_file"]]
    validation_splicing_file = args[["validation_splicing_file"]]
    splicing_file = args[["splicing_file"]]
    genexpr_file = args[["genexpr_file"]]
    spldep_file = args[["spldep_file"]]
    validation_spldep_file = args[["validation_splicing_file"]]
    validation_od_file = args[["validation_od_file"]]
    prolif_file = args[["prolif_file"]]
    exon_info_file = args[["exon_info_file"]]
    domains_file = args[["domains_file"]]
    sso_seqs_file = args[["sso_seqs_file"]]
    sso_optimization_file = args[["sso_optimization_file"]]
    figs_dir = args[["figs_dir"]]
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    event_info = read_tsv(event_info_file)
    annot = read_tsv(annotation_file) %>%
        mutate(event_gene = paste0(EVENT,"_",GENE))
    events_genes = annot %>%
        dplyr::select(EVENT,GENE,ENSEMBL,event_gene)
    
    splicing = read_tsv(splicing_file)
    genexpr = read_tsv(genexpr_file)
    spldep = read_tsv(spldep_file) %>% dplyr::rename(EVENT=index)
    
    validation_psi = read_tsv(validation_splicing_file)
    validation_clonogenic = read_tsv(validation_od_file) %>%
        group_by(DepMap_ID, CCLE_Name, event_gene, replicate_technical, replicate_biological) %>%
        summarize(od = mean(od)) %>% # summarize OD replicates
        ungroup() %>%
        separate(event_gene, c("EVENT","GENE"), remove=FALSE)
    prolif = read_tsv(prolif_file)

    exon_info = read_tsv(exon_info_file)
    domains = read_tsv(domains_file)
    
    sso_seqs = read_tsv(sso_seqs_file) %>%
        left_join(exon_info, by=c("EVENT"="exon_vastdb","GENE"="gene")) %>%
        left_join(event_info %>% distinct(EVENT,COORD_o,REF_CO), by="EVENT")
    
    sso_optimization = read_tsv(sso_optimization_file)
    
    events_val = validation_psi %>% pull(EVENT) %>% unique()
    
    # compute harm scores
    ## using experimentally measured (RT-PCR) delta PSIs
    validation_harm_scores = spldep %>%
        filter(EVENT %in% events_val) %>%
        pivot_longer(-EVENT, names_to="DepMap_ID", values_to="spldep") %>%
        left_join(validation_psi, by=c("DepMap_ID","EVENT")) %>% 
        # compute harm score
        mutate(
            delta_psi = psi_treated - psi_untreated,
            harm_score = (-1) * delta_psi * spldep
        ) %>%
        left_join(
            splicing %>%
                pivot_longer(-EVENT, names_to="DepMap_ID", values_to="psi_untreated_rnaseq"),
            by = c("EVENT","DepMap_ID")
        )
    
    # normalize OD by biological replicate
    od_ctl_neg = validation_clonogenic %>%
        filter(event_gene=="CONTROL_NEG") %>%
        group_by(CCLE_Name, replicate_biological) %>%
        summarize(od_ctl_neg = mean(od)) %>%
        ungroup()
    
    validation_clonogenic = validation_clonogenic %>%
        left_join(od_ctl_neg, by=c("CCLE_Name","replicate_biological")) %>%
        mutate(
            od_norm = od / od_ctl_neg,
            replicate_biological = factor(replicate_biological, levels=1:3)
        ) %>%
        group_by(CCLE_Name, replicate_biological, event_gene) %>%
        mutate(od_norm_avg = mean(od_norm)) %>%
        ungroup()
    
    # subset
    extras = c("CONTROL_NEG","EMPTY","WATER","CONTROL_POS")
    validation_clonogenic = validation_clonogenic %>%
        filter(CCLE_Name %in% SELECTED_CELL_LINES) %>%
        filter(event_gene %in% c(extras,SELECTED_EXONS)) %>%
        mutate(
            CCLE_Name = factor(CCLE_Name, levels=SELECTED_CELL_LINES),
            event_gene = factor(event_gene, levels=c(extras,SELECTED_EXONS)) 
        )

    validation_harm_scores = validation_harm_scores %>%
        filter(CCLE_Name %in% SELECTED_CELL_LINES) %>%
        filter(event_gene %in% c(extras,SELECTED_EXONS)) %>%
        mutate(
            CCLE_Name = factor(CCLE_Name, levels=SELECTED_CELL_LINES),
            event_gene = factor(event_gene, levels=c(extras,SELECTED_EXONS)) 
        )
    
    # average and normalize cell proliferation OD
    prolif = prolif %>%
        group_by(CCLE_Name, DepMap_ID, replicate_technical, replicate_biological, timepoint) %>%
        # average measurements
        summarize(
            od = mean(od)
        ) %>%
        ungroup() %>%
        group_by(CCLE_Name, DepMap_ID, replicate_biological, timepoint) %>%
        # average technical replicates
        mutate(
            od_avg = mean(od)
        ) %>%
        ungroup()
    
    prolif_ctl = prolif %>%
        filter(timepoint == 0) %>%
        rename(od_avg_ctl = od_avg) %>%
        dplyr::select(-c(od, timepoint, replicate_technical)) %>%
        distinct()
    
    prolif = prolif %>%
        left_join(
            prolif_ctl, 
            by=c("CCLE_Name","DepMap_ID","replicate_biological")
        ) %>%
        mutate(
            od_norm = od / od_avg_ctl
        ) %>%
        group_by(CCLE_Name, replicate_biological, timepoint) %>%
        mutate(
            od_norm_avg = mean(od_norm)
        ) %>%
        ungroup()
    
    ## prepare annotated proteoform domains
    domains = domains %>%
        drop_na(description) %>% 
        mutate(
            width = end - start,
            description = sprintf("%s (%s)", description, type) 
        ) %>%
        filter(width>0) %>%
        group_by(seq_region_name, description) %>%
        slice_max(width, n=1) %>%
        ungroup() %>%
        distinct(seq_region_name, start, end, description)
    
    exon_info = exon_info %>%
        mutate(
            proteoform_lab = sprintf("%s_%s | %s", exon_vastdb, gene, main_inclusion_proteoform_ensembl)
        ) %>%
        left_join(
            domains,
            by=c("main_inclusion_proteoform_ensembl"="seq_region_name")
        ) %>% 
        drop_na(description) %>%
        group_by(proteoform_lab) %>%
        mutate(index = row_number()) %>%
        ungroup()
    
    # prep sso optimization
    sso_optimization = sso_optimization %>%
        distinct(event_gene, psi_treated, psi_scrambled, psi_untreated, aso_id, selected_treatment) %>%
        pivot_longer(-c(event_gene,aso_id,selected_treatment), names_to="condition", values_to="psi") %>%
        mutate(
            aso_id = factor(aso_id, levels=c("A1","A2","A1+A2","A3")),
            event_gene = factor(event_gene, levels=SELECTED_EXONS),
            condition = case_when(
                condition=="psi_treated" ~ "+SSOs",
                condition=="psi_untreated" ~ "EMPTY",
                condition=="psi_scrambled" ~ "CONTROL_NEG"
            ),
            condition = factor(condition, levels=c("CONTROL_NEG","EMPTY","+SSOs"))
        )
    
    # plot
    plts = make_plots(validation_clonogenic, validation_harm_scores, prolif, exon_info, sso_optimization, sso_seqs)
    
    # make figdata
    figdata = make_figdata(validation_clonogenic, validation_harm_scores, prolif, exon_info, sso_optimization, sso_seqs)

    # make source data
    source_data = make_source_data(plts)
    
    # save
    save_plots(plts, figs_dir)
    save_figdata(figdata, figs_dir)
    save_source_data(source_data, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}