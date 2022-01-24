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

require(tidyverse)
require(ggpubr)
require(cowplot)
require(extrafont)
require(scattermore)

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

# Development
# -----------
# PREP_DIR = file.path(ROOT,'data','prep')
# RESULTS_DIR = file.path(ROOT,'results','splicing_dependency_tcga')
# embedding_file = file.path(RESULTS_DIR,'files','PANCAN','embedded_splicing_dependency_mean-EX.tsv.gz')
# metadata_file = file.path(PREP_DIR,'metadata','PANCAN.tsv.gz')
# figs_dir = file.path(RESULTS_DIR,'figures','personalization')

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

make_plots = function(embedding, metadata){
    plts = list(
        plot_patient_clusters(embedding, metadata)
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
    
    # plot
    plts = make_plots(embedding, metadata)

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
