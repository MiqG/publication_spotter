#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# 
# Outline
# -------
# - prediction of fitness effect of exon CRISPR screen Thomas 2020: if we run the model with less data, do we get the same predictions? (How robust are the predictions?)

require(tidyverse)
require(ggrepel)
require(ggpubr)

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

RANDOM_SEED = 1234
THRESH_LR_PVALUE = 0.001

# Development
# -----------
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# RESULTS_DIR = file.path(ROOT,'results','model_splicing_dependency')
# protein_impact_file = file.path(RAW_DIR,'VastDB','PROT_IMPACT-hg38-v3.tab.gz')
# annotation_file = file.path(RAW_DIR,'VastDB','EVENT_INFO-hg38_noseqs.tsv')
# psi_file = file.path(RAW_DIR,'articles','Thomas2020','vast_out','PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz')
# genexpr_file = file.path(RAW_DIR,'articles','Thomas2020','vast_out','TPM-hg38-2.tab.gz')
# crispr_screen_file = file.path(PREP_DIR,'Thomas2020','crispr_screen.tsv.gz')
# models_file = file.path(RESULTS_DIR,'files','models_gene_dependency-EX.tsv.gz')
# total=100
# figs_dir = file.path(RESULTS_DIR,'figures','validation_crispr')



##### FUNCTIONS #####
sample_splicing_dependencies = function(models, psi, genexpr, total, random_seed=NULL){
    set.seed(random_seed)
    
    samples = setdiff(colnames(psi), 'EVENT')
    
    models = models %>% column_to_rownames('EVENT')
    psi = psi %>% column_to_rownames('EVENT')
    genexpr = genexpr %>% dplyr::select(-NAME) %>% column_to_rownames('ID')
    genexpr = log2(genexpr+1)
    
    result = sapply(rownames(models), function(event){
        gene = models[event,] %>% pull(ENSEMBL)
        
        # sample coefficients
        intercept_coefs = rnorm(n=total, mean=pull(models[event,],intercept_coefficient), sd=pull(models[event,],intercept_stderr))
        event_coefs = rnorm(n=total, mean=pull(models[event,],event_coefficient), sd=pull(models[event,],event_stderr))
        gene_coefs = rnorm(n=total, mean=pull(models[event,],gene_coefficient), sd=pull(models[event,],gene_stderr))
        interaction_coefs = rnorm(n=total, mean=pull(models[event,],interaction_coefficient), sd=pull(models[event,],interaction_stderr))
        
        # make predictions for each combination of coefficients
        psi_mean = models[event,] %>% pull(event_mean)
        psi_std = models[event,] %>% pull(event_std)
        genexpr_mean = models[event,] %>% pull(gene_mean)
        genexpr_std = models[event,] %>% pull(gene_std)
        
        preds = sapply(samples, function(sample_oi){
            x_psi = (psi[event,sample_oi] - psi_mean)/psi_std
            x_genexpr = (genexpr[gene,sample_oi] - genexpr_mean)/genexpr_std
            x_interaction = x_psi * x_genexpr
            y_pred = intercept_coefs + event_coefs*x_psi + gene_coefs*x_genexpr + interaction_coefs*x_interaction
            return(y_pred)
        }, simplify=FALSE)
        
        return(preds)
    }, simplify=FALSE)
    
    # switch order
    result = sapply(samples, function(sample_oi){ 
        r = sapply(result, function(r){r[[sample_oi]]},simplify=FALSE) 
        r = do.call(rbind,r)
        return(r)
    }, simplify=FALSE)
    
    return(result)
}


plot_crispr_validation = function(crispr_screen, spldeps, models){
    events_considered = intersect(crispr_screen[['EVENT']],rownames(spldeps[[1]]))
    
    # compute distribution of correlations between predictions and real values
    comparisons = crispr_screen %>% pull(comparison) %>% unique()
    corrs = sapply(names(spldeps), function(sample_oi){
        preds = spldeps[[sample_oi]]
        
        corrs_sample = sapply(comparisons, function(comparison_oi){
            true = crispr_screen %>% 
                filter(cell_line %in% sample_oi & 
                       comparison %in% comparison_oi) %>% 
                column_to_rownames('EVENT')
        
            # measure correlations with sampled splicing dependencies
            corrs_comparison = apply(preds, 2, function(pred){ 
                x = pred
                y = true[rownames(preds),'fitness_score']
                
                df = tryCatch({
                    test = cor.test(x,y,method='spearman',use='pairwise.complete.obs')
                    df = data.frame(rho=test$estimate[['rho']], pvalue=test$p.value)
                    return(df)
                }, error=function(e){
                    df = data.frame(rho=NA, pvalue=NA)
                    df
                })
                return(df)
            })
            corrs_comparison = do.call(rbind,corrs_comparison)
            corrs_comparison[['comparison']] = comparison_oi
            return(corrs_comparison)
        }, simplify=FALSE)
        corrs_sample = do.call(rbind,corrs_sample)
        corrs_sample[['cell_line']] = sample_oi

        return(corrs_sample)
    },simplify=FALSE)
    corrs = do.call(rbind,corrs)
   
    
    plts = list()
    # CRISPR effects on selected exons
    plts[['crispr_selected-scatter']] = crispr_screen %>%
        left_join(models, by='EVENT') %>%
        mutate(log10_pvalue=-log10(pvalue)) %>%
        filter(is_selected) %>%
        ggscatter(x='event_gene', y='fitness_score', size='log10_pvalue',
                  color='cell_line', position=position_dodge(0.2)) +
        facet_wrap(~comparison, ncol=2) +
        labs(x='Event & Gene', y='Fitness Score', 
             color='Cell Line', size='-log10(p-value)') +
        theme_pubr(x.text.angle = 70, border=TRUE) +
        geom_hline(yintercept=0, linetype='dashed')
    
    # distributions of correlations with predictions (sampling size)
    plts[['predictions-corrs_violins']] = corrs %>% 
        #filter(comparison %in% c('0d-vs-14d','0d-vs-8d')) %>% 
        ggviolin(x='comparison', y='rho', facet.by='cell_line', 
                 fill='cell_line', color=FALSE) + 
        geom_boxplot(width=0.5) +
        guides(fill='none') +
        labs(x='Comparison', y='Spearman Correlation') +
        stat_compare_means(method='wilcox.test') + 
        theme_pubr(x.text.angle = 70)
    
    # distributions of pvalues of correlations with predictions
    plts[['predictions-corrs_pvalues']] = corrs %>% 
        # filter(comparison %in% c('0d-vs-14d','0d-vs-8d')) %>% 
        gghistogram(x='pvalue', facet.by=c('comparison','cell_line'), bins=50)
    
    # scatter of the average predictions
    spldeps_long = lapply(names(spldeps), function(sample_oi){
        spldeps[[sample_oi]] %>% 
        as.data.frame() %>% 
        rownames_to_column('EVENT') %>% 
        pivot_longer(-EVENT, names_to='iteration', values_to='pred_fitness') %>%
        mutate(cell_line=sample_oi)
    })
    spldeps_long = do.call(rbind, spldeps_long)
    
    X = spldeps_long %>% 
        group_by(cell_line,EVENT) %>%
        summarize(pred_mean=mean(pred_fitness,na.rm=TRUE),
                  pred_max=quantile(pred_fitness,0.975,na.rm=TRUE),
                  pred_min=quantile(pred_fitness,0.025,na.rm=TRUE)) %>% 
        left_join(crispr_screen, by=c('EVENT','cell_line')) %>% 
        left_join(models, by='EVENT') %>%
        # filter(comparison %in% c('0d-vs-14d','0d-vs-8d')) %>%
        mutate(diff=fitness_score - pred_mean) %>%
        ungroup()
    
    
    plts[['predictions-means_scatter']] = X %>%
        ggplot(aes(x=fitness_score, y=pred_mean)) +
        facet_wrap(comparison~cell_line, ncol=2) +
        geom_point(aes(color=is_selected), size=1) +
        geom_errorbar(aes(ymin=pred_min, ymax=pred_max, color=is_selected), 
                      width=0.05) +
        theme_pubr(border=TRUE) +    
        ggpubr::color_palette('jco') +
        stat_cor(aes(color=is_selected),method='spearman') +
        geom_abline(intercept = 0, slope = 1, linetype='dashed') + 
        geom_text_repel(
            aes(label=event_gene),
            X %>%
            group_by(comparison,cell_line) %>%
            slice_max(order_by = abs(diff), n=3)
        ) +
        labs(x='True Splicing Dependency', y='Predicted Splicing Dependency')

    plts[['predictions-diff_vs_std']] = X %>%
        ggscatter(x='event_std', y='diff', 
                  facet.by=c('comparison','cell_line')) +
        labs(x='Event PSI Std. across CCLE', y='Diff. Splicing Dependency') +
        geom_text_repel(
            aes(label=event_gene), 
            X %>% 
            group_by(comparison,cell_line) %>%
            slice_max(order_by=abs(diff), n=3)
        )
    
    # what are the events that we could map and compare?
    plts[['predictions-protein_impact']] = X %>% 
        distinct(EVENT,target_type,ONTO) %>% 
        count(target_type,ONTO) %>% 
        ggbarplot(x='ONTO', y='n', fill='target_type', position=position_dodge(0.8),
                  color=FALSE, label=TRUE, palette='jco') + 
        theme_pubr(x.text.angle=45, border=TRUE) +
        labs(x='Predicted Protein Impact', y='Count', fill='CRISPR Target') +
        ylim(0,50)
    
    return(plts)
}


make_plots = function(crispr_screen, spldeps, models){
    plts = list(
        plot_crispr_validation(crispr_screen, spldeps, models)
    )
    plts = do.call(c,plts)
    return(plts)
}


save_plt = function(plts, plt_name, extension='.pdf', 
                      directory='', dpi=350, 
                      width = par("din")[1], height = par("din")[2]){
        filename = file.path(directory,paste0(plt_name,extension))
        cowplot::save_plot(filename, 
                  plts[[plt_name]], 
                  base_width=width, base_height=height, dpi=dpi)
}


save_plots = function(plts, figs_dir){
    save_plt(plts, 'crispr_selected-scatter', '.pdf', figs_dir, width=7, height=7)
    save_plt(plts, 'predictions-corrs_violins', '.pdf', figs_dir, width=10, height=6)
    save_plt(plts, 'predictions-corrs_pvalues', '.pdf', figs_dir, width=10, height=10)
    save_plt(plts, 'predictions-means_scatter', '.pdf', figs_dir, width=10, height=10)
    save_plt(plts, 'predictions-diff_vs_std', '.pdf', figs_dir, width=10, height=10)
    save_plt(plts, 'predictions-protein_impact', '.pdf', figs_dir, width=7, height=5)
}


main = function(){
    args = getParsedArgs()

    crispr_screen_file = args$crispr_screen_file
    psi_file = args$psi_file
    genexpr_file = args$genexpr_file
    models_file = args$models_file
    annotation_file = args$annotation_file
    protein_impact_file = args$protein_impact_file
    total = args$total
    figs_dir = args$figs_dir
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    ## Spotter
    models = read_tsv(models_file) %>%
        mutate(
            is_selected = pearson_correlation>0 & lr_pvalue<THRESH_LR_PVALUE,
            event_gene = paste0(EVENT,'_',GENE),
            event_type = ifelse(grepl('EX',EVENT),'EX','ALT')
        )
    ## VastDB
    annot = read_tsv(annotation_file)
    protein_impact = read_tsv(protein_impact_file) %>%
        dplyr::rename(EVENT=EventID)
    ## Thomas 2020
    crispr_screen = read_tsv(crispr_screen_file) %>%
        left_join(protein_impact, by="EVENT")
    psi = read_tsv(psi_file) %>% 
        dplyr::rename(HeLa=SRR7946515_1, PC9=SRR7946516_1)
    genexpr = read_tsv(genexpr_file) %>% 
        dplyr::rename(HeLa=SRR7946515_1, PC9=SRR7946516_1)
    
    tmp = crispr_screen %>% 
        left_join(models, by='EVENT') %>%
        left_join(
            psi %>% 
            pivot_longer(c(HeLa,PC9), 
                         names_to = 'cell_line', 
                         values_to = 'psi'), by=c('EVENT','cell_line')) %>% 
        left_join(
            genexpr %>% 
            dplyr::rename(ENSEMBL=ID) %>%
            pivot_longer(c(HeLa,PC9), 
                         names_to = 'cell_line', 
                         values_to = 'genexpr'), by=c('ENSEMBL','cell_line')) %>%
        mutate(genexpr=log2(genexpr+1))
    
    # filter available data
    common_events = Reduce(intersect, list(crispr_screen[['EVENT']], psi[['EVENT']], models[['EVENT']]))
    models = models %>% filter(EVENT %in% common_events)
    psi = psi %>% filter(EVENT %in% common_events)
    genexpr = genexpr %>% filter(ID %in% models[['ENSEMBL']])
    
    # sample splicing dependencies
    spldeps = sample_splicing_dependencies(models, psi, genexpr, total, RANDOM_SEED)
    
    # plot
    plts = make_plots(crispr_screen, spldeps, models)

    # save
    save_plots(plts, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}