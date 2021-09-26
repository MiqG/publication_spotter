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
require(org.Hs.eg.db)

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

# variables
THRESH_FDR = 0.05
THRESH_MEDIAN_DIFF = 5
EVENT_TYPE = 'EX' #!
ORGDB = org.Hs.eg.db

# Development
# -----------
RAW_DIR = file.path(ROOT,'data','raw')
PREP_DIR = file.path(ROOT,'data','prep')
RESULTS_DIR = file.path(ROOT,'results','splicing_dependency_tcga')
diff_psi_file = file.path(RESULTS_DIR,'files','READ','mannwhitneyu-psi-PrimaryTumor_vs_SolidTissueNormal.tsv.gz')
diff_spldep_file = file.path(RESULTS_DIR,'files','READ','mannwhitneyu-spldep-PrimaryTumor_vs_SolidTissueNormal-EX.tsv.gz')
annotation_file = file.path(RAW_DIR,'VastDB','EVENT_INFO-hg38_noseqs.tsv')
msigdb_dir = file.path(RAW_DIR,'MSigDB','msigdb_v7.4','msigdb_v7.4_files_to_download_locally','msigdb_v7.4_GMTs')
protein_impact_file = file.path(RAW_DIR,'VastDB','PROT_IMPACT-hg38-v3.tab.gz')


##### FUNCTIONS #####
load_diff_result = function(filename, pattern){
    df = read_tsv(filename) %>% 
        rename_all(recode, index = "EVENT") %>%
        column_to_rownames('EVENT') %>%
        setNames(paste0(pattern, names(.))) %>%
        rownames_to_column('EVENT') %>%
        filter(str_detect(EVENT,EVENT_TYPE))
    return(df)
}


get_genes_lists = function(diff_result){
    df = diff_result
    cols_oi = c('psi__is_significant','spldep__is_significant')
    
    genes_lists = df[,c('GENE',cols_oi)] %>% 
        pivot_longer(all_of(cols_oi), names_to='dataset', values_to='is_significant') %>%
        mutate(dataset = gsub('__is_significant','',dataset)) %>%
        filter(is_significant) %>%
        drop_na() %>%
        distinct(GENE, dataset) %>%
        with(., split(GENE, dataset))
        
    return(genes_lists)
}


get_events_lists = function(diff_result){
    df = diff_result 
    cols_oi = c('psi__is_significant','spldep__is_significant')
    
    events_lists = df[,c('EVENT',cols_oi)] %>% 
        pivot_longer(all_of(cols_oi), names_to='dataset', values_to='is_significant') %>%
        mutate(dataset = gsub('__is_significant','',dataset)) %>%
        filter(is_significant) %>%
        drop_na() %>%
        distinct(EVENT, dataset) %>%
        with(., split(EVENT, dataset))
        
    return(events_lists)
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
        enrichments[['GO_BP']] = enrichGO(genes, OrgDb=ORGDB, keyType='SYMBOL', universe=universe[['genes']])
    }
    
    if(length(events)>0){
        enrichments[['protein_impact']] = enricher(events, TERM2GENE=ontologies[['protein_impact']], universe=universe[['events']])
    }
    return(enrichments)
}


run_enrichments = function(genes_lists, events_lists, universe, ontologies){
    enrichments = sapply(names(genes_lists), function(cond){
            genes = genes_lists[[cond]]
            events = events_lists[[cond]]
            run_enrichment(genes, events, universe, ontologies)
        }, simplify=FALSE)
    return(enrichments)
}


get_enrichment_result = function(enrich_list, thresh){
    ## datasets are extracted from names
    if (length(enrich_list)>0){
        ontos = names(enrich_list[[1]])
        datasets = names(enrich_list)
        results = sapply(
            ontos, function(onto){
            result = lapply(datasets, function(dataset){
                result = enrich_list[[dataset]][[onto]]
                if(!is.null(result)){
                    result = result@result
                    result[['dataset']] = dataset
                    result[['ontology']] = onto
                }
                return(result)
            })
            result[sapply(result, is.null)] = NULL
            result = do.call(rbind,result)
            ## filter by p.adjusted
            result = result %>% filter(p.adjust<thresh)
        }, simplify=FALSE)
        results = do.call(rbind,results)
    }else{
        results = data.frame(
            ID=NA, Description=NA, GeneRatio=NA, BgRatio=NA, pvalue=NA,
            p.adjust=NA, qvalue=NA, geneID=NA, Count=NA, dataset=NA, ontology=NA
        ) %>% drop_na()
    }
    
    return(results)
}


# prep_enrichments = function(enrichments){
#     if (length(enrichments)>0){
#         ontos = names(enrichments)
#         results = lapply(ontos, function(onto){
#             result = enrichments[[onto]]
#             if(!is.null(result)){
#                 result = result@result
#                 result[['ontology']] = onto
#             }
#             return(result)
#         })
#         results[sapply(results, is.null)] = NULL
#         results = do.call(rbind,results)
#         ## filter by p.adjusted
#         results = results %>% filter(p.adjust<THRESH_FDR)
#     }else{
#         results = data.frame(
#             ID=NA, Description=NA, GeneRatio=NA, BgRatio=NA, pvalue=NA,
#             p.adjust=NA, qvalue=NA, geneID=NA, Count=NA, ontology=NA
#         ) %>% drop_na()
#     }
    
#     return(results)
# }


main = function(){
    args = getParsedArgs()
    diff_psi_file = args$diff_psi_file
    diff_spldep_file = args$diff_spldep_file
    annotation_file = args$annotation_file
    msigdb_dir = args$msigdb_dir
    protein_impact_file = args$protein_impact_file
    output_file = args$output_file
    
    # load
    diff_psi = load_diff_result(diff_psi_file,'psi__')
    diff_spldep = load_diff_result(diff_spldep_file,'spldep__')
    annot = read_tsv(annotation_file)
    ontologies = list(
        "hallmarks" = read.gmt(file.path(msigdb_dir,'h.all.v7.4.symbols.gmt')),
        "oncogenic_signatures" = read.gmt(file.path(msigdb_dir,'c6.all.v7.4.symbols.gmt')),
        "protein_impact" = read_tsv(protein_impact_file) %>%
            dplyr::rename(EVENT=EventID, term=ONTO) %>%
            dplyr::select(term,EVENT)
    )
    
    # prepare diff_result
    diff_result = merge(diff_psi, diff_spldep, all=TRUE, by='EVENT') %>%
        left_join(annot[,c('EVENT','GENE')], by='EVENT') %>%
        mutate(psi__is_significant = psi__padj<THRESH_FDR & 
                                     abs(psi__median_diff)>THRESH_MEDIAN_DIFF,
               spldep__is_significant = spldep__padj<THRESH_FDR &
                                        psi__is_significant)

    # run enrichments
    genes_lists = get_genes_lists(diff_result)
    events_lists = get_events_lists(diff_result)
    universe = get_universe(diff_result)
    enrichments = run_enrichments(genes_lists, events_lists, universe, ontologies)
    results_enrich = get_enrichment_result(enrichments, THRESH_FDR)
    
    # save
    write_tsv(results_enrich, output_file)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}