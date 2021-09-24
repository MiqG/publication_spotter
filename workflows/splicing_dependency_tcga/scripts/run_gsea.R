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
require(clusterProfiler)

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

# variables
THRESH_FDR = 0.05

# Development
# -----------
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# RESULTS_DIR = file.path(ROOT,'results','splicing_dependency_tcga')
# diff_result_file = file.path(RESULTS_DIR,'files','HNSC','mannwhitneyu-spldep-RESPONDER_vs_NONRESPONDER-EX.tsv.gz')
# annotation_file = file.path(RAW_DIR,'VastDB','EVENT_INFO-hg38_noseqs.tsv')
# msigdb_dir = file.path(RAW_DIR,'MSigDB','msigdb_v7.4','msigdb_v7.4_files_to_download_locally','msigdb_v7.4_GMTs')
# protein_impact_file = file.path(RAW_DIR,'VastDB','PROT_IMPACT-hg38-v3.tab.gz')


##### FUNCTIONS #####
get_genes_list = function(diff_result){
    df = diff_result 
    genes_list = df %>%
        filter(padj < THRESH_FDR) %>%
        pull(GENE) %>%
        unique()
        
    return(genes_list)
}


get_events_list = function(diff_result){
    df = diff_result 
    events_list = df %>%
        filter(padj < THRESH_FDR) %>%
        pull(EVENT) %>%
        unique()
        
    return(events_list)
}


get_universe = function(diff_result){
    df = diff_result
    universe = df[,c('EVENT','GENE')] %>% apply(., 2, unique)
    names(universe) = c('events','genes')
    return(universe)
}


run_enrichment = function(genes, events, universe, ontologies){
    enrichments = list()
    if(length(genes)>0){
        enrichments[['hallmarks']] = enricher(genes, TERM2GENE=ontologies[['hallmarks']], universe=universe[['genes']])
    enrichments[['oncogenic_signatures']] = enricher(genes, TERM2GENE=ontologies[['oncogenic_signatures']], universe=universe[['genes']])
    enrichments[['GO_BP']] = enricher(genes, TERM2GENE=ontologies[['GO_BP']], universe=universe[['genes']])
    }
    
    if(length(events)>0){
        enrichments[['protein_impact']] = enricher(events, TERM2GENE=ontologies[['protein_impact']], universe=universe[['events']])
    }
    return(enrichments)
}


prep_enrichments = function(enrichments){
    if (length(enrichments)>0){
        ontos = names(enrichments)
        results = lapply(ontos, function(onto){
            result = enrichments[[onto]]
            if(!is.null(result)){
                result = result@result
                result[['ontology']] = onto
            }
            return(result)
        })
        results[sapply(results, is.null)] = NULL
        results = do.call(rbind,results)
        ## filter by p.adjusted
        results = results %>% filter(p.adjust<THRESH_FDR)
    }else{
        results = data.frame(
            ID=NA, Description=NA, GeneRatio=NA, BgRatio=NA, pvalue=NA,
            p.adjust=NA, qvalue=NA, geneID=NA, Count=NA, ontology=NA
        ) %>% drop_na()
    }
    
    return(results)
}


main = function(){
    args = getParsedArgs()
    diff_result_file = args$diff_result_file
    annotation_file = args$annotation_file
    msigdb_dir = args$msigdb_dir
    protein_impact_file = args$protein_impact_file
    output_file = args$output_file
    
    # load
    diff_result = read_tsv(diff_result_file) %>% 
        rename_all(recode, index = "EVENT")
    annot = read_tsv(annotation_file)
    ontologies = list(
        "hallmarks" = read.gmt(file.path(msigdb_dir,'h.all.v7.4.symbols.gmt')),
        "oncogenic_signatures" = read.gmt(file.path(msigdb_dir,'c6.all.v7.4.symbols.gmt')),
        "GO_BP" = read.gmt(file.path(msigdb_dir,'c5.go.bp.v7.4.symbols.gmt')),
        "protein_impact" = read_tsv(protein_impact_file) %>%
                            dplyr::rename(EVENT=EventID, term=ONTO) %>%
                            dplyr::select(term,EVENT)
    )
    
    # prepare diff_result
    diff_result = diff_result %>% 
        left_join(annot[,c('EVENT','GENE')], by='EVENT') %>%
        filter(str_detect(EVENT,'EX')) # only EX events
    
    # run enrichments
    genes_list = get_genes_list(diff_result)
    events_list = get_events_list(diff_result)
    universe = get_universe(diff_result)
    enrichments = run_enrichment(genes_list, events_list, universe, ontologies)
    results_enrich = prep_enrichments(enrichments)
    
    # save
    write_tsv(results_enrich, output_file)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}