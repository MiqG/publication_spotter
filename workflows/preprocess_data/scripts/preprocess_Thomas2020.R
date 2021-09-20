#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Prepare metadata from thomas 2020

require(readxl)
require(tidyverse)

ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))

# Development
# -----------
# RAW_DIR = file.path(ROOT,'data','raw')
# thomas_crispr_screen_file = file.path(RAW_DIR,'articles','Thomas2020','crispr_screen.xlsx')
# thomas_event_mapping_file = file.path(RAW_DIR,'articles','Thomas2020','event_mapping_vastdb.tsv')


##### FUNCTIONS #####
prep_thomas = function(X){
    cols = grep('\\|', colnames(X), value = TRUE)
    conditions = do.call(rbind,strsplit(cols, '\\|'))
    stats = unique(conditions[,1])
    conditions = unique(conditions[,2])
    
    screen = lapply(conditions, function(cond){
        screen = X %>% dplyr::select(contains(cond))
        colnames(screen) = stats
        cond = unlist(strsplit(cond, '\\.'))
        screen = screen %>% 
            mutate(cell_line=cond[1], 
                   comparison=cond[2], 
                   replicate=cond[3])
        cols_oi = c('event','EVENT','target_type','geneName')
        screen = cbind(X[,cols_oi],screen)
        return(screen)
    })
    screen = do.call(rbind,screen)
    return(screen)
}


load_crispr_screen = function(
    thomas_crispr_screen_file,
    thomas_event_mapping_file
){
    thomas_crispr_screen = read_excel(thomas_crispr_screen_file, skip=1, sheet=5)
    thomas_mapping = read_tsv(thomas_event_mapping_file)

    thomas_crispr_screen = thomas_crispr_screen %>% 
        left_join(thomas_mapping, by=c('event','target_type')) %>%
        filter(!is.na(EVENT))
    thomas_crispr_screen = prep_thomas(thomas_crispr_screen)
    thomas_crispr_screen = thomas_crispr_screen %>% mutate(
        fitness_score = log2(fc_norm_mean),
        study = 'Thomas (2020)',
        pvalue = pval,
        hit = fdr < 0.05
    )
    
    return(thomas_crispr_screen)
}
    

main = function(){
    args = getParsedArgs()
    thomas_crispr_screen_file = args$thomas_crispr_screen_file
    thomas_event_mapping_file = args$thomas_event_mapping_file
    output_file = args$output_file
    
    crispr_screen = load_crispr_screen(
        thomas_crispr_screen_file, 
        thomas_event_mapping_file
    )

    write_tsv(crispr_screen, output_file)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}