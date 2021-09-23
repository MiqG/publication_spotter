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
require(writexl)
require(ggpubr)
require(cowplot)
require(scattermore)
require(ggrepel)

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

# variables
MIN_OBS = 50
THRESH_ZSCORE = 1.96
THRESH_PVALUE = 0.05

PALETTE_VARS = 'uchicago'

METADATA_OI = c(
    'sex',
    'source',
    'culture_type',
    'culture_medium',
    'sample_collection_site',
    'primary_or_metastasis',
    'age',
    'cancer_type',
    'mitotic_index',
    'stemness'
)

METADATA_OI_CONT = c(
    'age',
    'mitotic_index',
    'stemness'
)

METADATA_PALETTES = list(
    'sex' = 'lancet',
    'source' = 'default',
    'culture_type' = 'Set2',
    'culture_medium' = 'Dark2',
    'sample_collection_site' = 'simpsons',
    'primary_or_metastasis' = 'npg',
    'age' = 'Greens',
    'cancer_type' = 'Paired',
    'mitotic_index' = 'Reds',
    'stemness' = 'Blues'
)

# Development
# -----------
# PREP_DIR = file.path(ROOT,'data','prep')
# RESULTS_DIR = file.path(ROOT,'results','model_splicing_dependency')
# psi_file = file.path(PREP_DIR,'event_psi','CCLE-EX.tsv.gz')
# genexpr_file = file.path(PREP_DIR,'genexpr_tpm','CCLE.tsv.gz')
# rnai_file = file.path(PREP_DIR,'demeter2','CCLE.tsv.gz')
# splicing_dependency_file = file.path(RESULTS_DIR,'files','splicing_dependency_mean-EX.tsv.gz')
# embedded_dependency_file = file.path(RESULTS_DIR,'files','embedded_splicing_dependency_mean-EX.tsv.gz')
# metadata_file = file.path(PREP_DIR,'metadata','CCLE.tsv.gz')
# indices_file = file.path(PREP_DIR,'transcriptome_indices','CCLE.tsv.gz')
# genes_splicing_factors_file = file.path(ROOT,'support','splicing_factors-symbols.txt')
# figs_dir = file.path(RESULTS_DIR,'figures','splicing_dependency_embeddings')

##### FUNCTIONS #####
diff_splicing_dependency = function(spldep, embedding){
    X = spldep %>% column_to_rownames('index')
    clusters = unique(embedding[['leiden_labels']])
    
    results = lapply(clusters, function(cluster){
        samples_in = embedding %>% 
            filter(leiden_labels %in% cluster) %>% 
            pull(index)
        samples_out = embedding %>% 
            filter(!(leiden_labels %in% cluster)) %>% 
            pull(index)
        res = lapply(rownames(X), function(event){
            x_in = as.numeric(X[event,samples_in])
            x_out = as.numeric(X[event,samples_out])
            
            df = tryCatch({
                    test = wilcox.test(x_in, x_out)
                    data.frame(
                        EVENT = event,
                        pvalue = test$p.value,
                        median_diff = median(x_in, na.rm=TRUE) - median(x_out, na.rm=TRUE),
                        median_in = median(x_in, na.rm=TRUE),
                        median_out = median(x_out, na.rm=TRUE),
                        std_in = sd(x_in, na.rm=TRUE),
                        std_out = sd(x_out, na.rm=TRUE)
                    )
                }, error = function(e){
                    data.frame(
                        EVENT = event,
                        pvalue = NA,
                        median_diff = median(x_in, na.rm=TRUE) - median(x_out, na.rm=TRUE),
                        median_in = median(x_in, na.rm=TRUE),
                        median_out = median(x_out, na.rm=TRUE),
                        std_in = sd(x_in, na.rm=TRUE),
                        std_out = sd(x_out, na.rm=TRUE)
                    )
                })
            
            return(df)
        })
        res = do.call(rbind, res)
        res[['comparison']] = paste0(cluster,'_vs_rest')
        return(res)
    })
    results = do.call(rbind, results)
    return(results)
}


plot_embeddings = function(embedding, metadata, pattern, figtitle){
    X = embedding %>% left_join(metadata, by='index')
    
    plts = list()
    plts[['pca']] = X %>% 
        ggplot(aes_string(x='PC0', y='PC1')) + 
        geom_scattermore(pixels=c(1000,1000), pointsize = 2, alpha=0.8) +
        theme_pubr() +
        labs(title = figtitle)
    
    plts[['umap']] = X %>% 
        ggplot(aes_string(x='UMAP0', y='UMAP1')) + 
        geom_scattermore(pixels=c(1000,1000), pointsize = 2, alpha=0.8) +
        theme_pubr() +
        labs(title = figtitle)
    
    plts[['umap_clusters']] = X %>%
        ggplot(aes_string(x='UMAP0', y='UMAP1', color='leiden_labels')) + 
        geom_scattermore(pixels=c(1000,1000), pointsize = 2, alpha=0.8) +
        theme_pubr() +
        labs(title = figtitle) +
        ggpubr::color_palette(
            palette=get_palette('jco', length(unique(X[['leiden_labels']])))
        )
    
    
    # with metadata information    
    tmp = sapply(METADATA_OI, function(col_oi){
        # make plot
        plt = X %>%
            drop_na(col_oi) %>%
            ggplot(aes_string(x='PC0', y='PC1', color=col_oi)) + 
            geom_scattermore(pixels = c(1000,1000), pointsize = 2, alpha=0.8) +
            theme_pubr() +
            labs(title = figtitle)
        # set palette
        if (col_oi %in% METADATA_OI_CONT){
            palette = METADATA_PALETTES[[col_oi]]
            plt = plt + gradient_color(palette)
        }else{
            palette = get_palette(METADATA_PALETTES[[col_oi]], 
                              length(unique(metadata[[col_oi]])))
            plt = set_palette(plt, palette=palette)
        }
        return(plt)
    }, simplify=FALSE)
    names(tmp) = paste0('pca-',names(tmp))
    plts = c(plts,tmp)
    
    tmp = sapply(METADATA_OI, function(col_oi){
        # make plot
        plt = X %>%
            drop_na(col_oi) %>%
            ggplot(aes_string(x='UMAP0', y='UMAP1', color=col_oi)) + 
            geom_scattermore(pixels = c(500,500), pointsize = 2, alpha=0.8) +
            theme_pubr() +
            labs(title = figtitle)
        # set palette
        if (col_oi %in% METADATA_OI_CONT){
            palette = METADATA_PALETTES[[col_oi]]
            plt = plt + gradient_color(palette)
        }else{
            palette = get_palette(METADATA_PALETTES[[col_oi]], 
                              length(unique(metadata[[col_oi]])))
            plt = set_palette(plt, palette=palette)
        }
        return(plt)
    }, simplify=FALSE)
    names(tmp) = paste0('umap-',names(tmp))
    plts = c(plts,tmp)
    
    # prepend pattern name
    names(plts) = paste0(pattern,'-',names(plts))
    
    return(plts)
}


plot_clusters = function(embedding, metadata, spldep){
    X = merge(metadata, embedding[,c('index','leiden_labels')], by='index', all = FALSE)
    
    col_oi = 'cancer_type'
    
    plts = list()
    # overview
    plts[['overview']] = X %>% 
        count(leiden_labels) %>% 
        ggbarplot(x='leiden_labels', y='n', label=TRUE, 
                  fill='leiden_labels', color=NA, 
                  palette=get_palette('jco', length(unique(X[['leiden_labels']])))) + 
        guides(fill="none") + 
        labs(x='Cluster', y='Count')
    
    for (col_oi in METADATA_OI){
        if (col_oi %in% METADATA_OI_CONT){
            # continous metadata
            palette = METADATA_PALETTES[[col_oi]]
            plt = X %>% 
                ggviolin(x='leiden_labels', y=col_oi, 
                         fill='leiden_labels', color=NA,
                         palette=get_palette('jco', length(unique(X[['leiden_labels']])))) +
                geom_boxplot(width=0.5) + 
                guides(fill="none") + 
                labs(x='Cluster')
            plts[[col_oi]] = plt
        }else{
            # for discrete metadata
            palette = get_palette(METADATA_PALETTES[[col_oi]], 
                              length(unique(metadata[[col_oi]])))
            plt = X %>% 
                group_by(leiden_labels,get(col_oi)) %>% 
                summarize(n=n()) %>% 
                mutate(freq=n/sum(n)) %>% 
                ggbarplot(x='leiden_labels', y='freq', fill='get(col_oi)', color=NA, 
                          palette=palette) +
                labs(x='Cluster', y='Proportion', fill=col_oi)
            plts[[col_oi]] = plt
        }
    }
    
    names(plts) = paste0('clusters-',names(plts))
    return(plts)
}



make_plots = function(embedding, metadata){
    plts = list(
        plot_embeddings(embedding, metadata, 'splicing_dependency', 'Splicing Dependency'),
        plot_clusters(embedding, metadata)
    )
    plts = do.call(c,plts)
    return(plts)
}


save_plt = function(plts, plt_name, extension='.pdf', 
                      directory='', dpi=350, 
                      width = par("din")[1], height = par("din")[2]){
        filename = file.path(directory,paste0(plt_name,extension))
        save_plot(filename, 
                  plts[[plt_name]], 
                  base_width=width, base_height=height, dpi=dpi)
}


save_plots = function(plts, figs_dir){
    # splicing dependencies
    plt_names = grep('splicing_dependency', names(plts), value = TRUE)
    sapply(plt_names, function(plt_name){
        save_plt(plts, plt_name, '.pdf', figs_dir, width=5, height=5)
    })
    
    # clusters
    plt_names = grep('clusters', names(plts), value = TRUE)
    sapply(plt_names, function(plt_name){
        save_plt(plts, plt_name, '.pdf', figs_dir, width=5, height=5)
    })
}


main = function(){
    args = getParsedArgs()
    metadata_file = args$metadata_file
    indices_file = args$indices_file
    embedded_dependency_file = args$embedded_dependency_file
    # splicing_dependency_file = args$splicing_dependency_file
    figs_dir = args$figs_dir
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    indices = read_tsv(indices_file)
    metadata = read_tsv(metadata_file) %>%
        dplyr::rename(index=DepMap_ID) %>%
        left_join(indices, by='index')
    embedding = read_tsv(embedded_dependency_file) %>%
        mutate(leiden_labels=as.factor(leiden_labels))
    # spldep = read_tsv(splicing_dependency_file)
    
    # plot
    plts = make_plots(embedding, metadata)

    # save
    save_plots(plts, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}