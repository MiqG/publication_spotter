#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# EDA of datasets used.
# 
# Outline
# -------
# - gene dependencies across cell lines
# - mRNA levels across cell lines
# - imputed PSI across cell lines
# - consider transcriptome indices: stemness

require(tidyverse)
require(ggpubr)
require(cowplot)
require(scattermore)

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

# variables
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
# RESULTS_DIR = file.path(ROOT,'results')
# PREP_DIR = file.path(ROOT,'data','prep')
# demeter_file = file.path(RESULTS_DIR,'eda_data','files','decompose_and_cluster','demeter2.tsv.gz')
# genexpr_file = file.path(RESULTS_DIR,'eda_data','files','decompose_and_cluster','genexpr.tsv.gz')
# psi_EX_file = file.path(RESULTS_DIR,'eda_data','files','decompose_and_cluster','psi-EX.tsv.gz')
# metadata_file = file.path(PREP_DIR,'metadata','CCLE.tsv')
# indices_file = file.path(PREP_DIR,'transcriptome_indices','CCLE.tsv.gz')

# figs_dir = file.path(ROOT,'results','eda_data','figures','eda')


##### FUNCTIONS #####
plot_metadata = function(metadata){
    X = metadata
    
    # discrete
    plts = sapply(setdiff(METADATA_OI,METADATA_OI_CONT), simplify=FALSE,
        function(col_oi){
            
            palette = get_palette(METADATA_PALETTES[[col_oi]], 
                                  length(unique(metadata[[col_oi]])))
            
            X %>%
            group_by_at(col_oi) %>%
            summarize(n=n()) %>%
            ggbarplot(x=col_oi, y='n', label = TRUE,
                      fill=col_oi, color=NA,
                      palette=palette) +
            guides(fill=FALSE) +
            labs(y='Count') +
            theme_pubr(x.text.angle = 70)
        }
    )
    
    # continuous
    plts[['age']] = X %>%
        gghistogram(x='age', fill='darkgreen', color=NA) +
        labs(y='Count')
    plts[['mitotic_index']] = X %>%
        gghistogram(x='mitotic_index', fill='darkred', color=NA) +
        labs(y='Count')
    plts[['stemness']] = X %>%
        gghistogram(x='stemness', fill='darkblue', color=NA) +
        labs(y='Count')
    
    # prepend plot name
    names(plts) = paste0('metadata-',names(plts))
    
    return(plts)
}


plot_dimreds = function(df, metadata, pattern, figtitle){
    X = df %>% left_join(metadata, by='index')
    
    plts = list()
    plts[['pca']] = 
    X %>% ggplot(aes_string(x='PC0', y='PC1')) + 
        geom_scattermore(pointsize = 2, alpha=0.8) +
        theme_pubr() +
        labs(title = figtitle)
    
    plts[['umap']] = X %>% 
        ggplot(aes_string(x='UMAP0', y='UMAP1')) + 
        geom_scattermore(pointsize = 2, alpha=0.8) +
        theme_pubr() +
        labs(title = figtitle)
    
    plts[['umap_clusters']] = X %>%
        ggplot(aes_string(x='UMAP0', y='UMAP1', color='leiden_labels')) + 
        geom_scattermore(pointsize = 2, alpha=0.8) +
        theme_pubr() +
        labs(title = figtitle)
    plts[['umap_clusters']] = set_palette(plts[['umap_clusters']], palette = 'Paired')
    
    
    # with metadata information    
    tmp = sapply(METADATA_OI, function(col_oi){
        # make plot
        plt = X %>%
            ggplot(aes_string(x='PC0', y='PC1', color=col_oi)) + 
            geom_scattermore(pointsize = 2, alpha=0.8) +
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
            ggplot(aes_string(x='UMAP0', y='UMAP1', color=col_oi)) + 
            geom_scattermore(pointsize = 2, alpha=0.8) +
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


make_plots = function(rnai, genexpr, psi, metadata){
    plts = list(
        plot_metadata(metadata),
        plot_dimreds(rnai, metadata, 'demeter2', 'Demeter2'),
        plot_dimreds(genexpr, metadata, 'genexpr', 'Gene Expression'),
        plot_dimreds(psi, metadata, 'psi_EX', 'PSI')
    )
    plts = do.call(c,plts)
    return(plts)
}


save_plt = function(plts, plt_name, extension='.pdf', 
                      directory='', dpi=350, 
                      width = par("din")[1], height = par("din")[2]){
        filename = file.path(directory,paste0(plt_name,extension))
        save_plot(filename, 
                  plts[[plt_name]], units='cm',
                  base_asp = 1,
                  base_width=width, base_height=height, dpi=dpi)
}


save_plots = function(plts, figs_dir){
    lapply(names(plts), function(plt_oi){
        save_plt(plts, plt_oi, '.pdf', figs_dir, 
                 width=12, height=12, dpi=350)
    })
}


main = function(){
    args = getParsedArgs()
    demeter_file = args$demeter_file
    genexpr_file = args$genexpr_file
    psi_EX_file = args$psi_EX_file
    metadata_file = args$metadata_file
    indices_file = args$indices_file
    figs_dir = args$figs_dir
    
    dir.create(figs_dir, recursive = TRUE)
    
    # load
    rnai = read_tsv(demeter_file) %>% mutate(leiden_labels=as.factor(leiden_labels))
    genexpr = read_tsv(genexpr_file) %>% mutate(leiden_labels=as.factor(leiden_labels))
    psi = read_tsv(psi_EX_file) %>% mutate(leiden_labels=as.factor(leiden_labels))
    indices = read_tsv(indices_file)
    metadata = read_tsv(metadata_file) %>%
        dplyr::rename(index=DepMap_ID) %>% 
        left_join(indices, by='index')
    
    # plot
    plts = make_plots(rnai, genexpr, psi, metadata)

    # save
    save_plots(plts, figs_dir)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}