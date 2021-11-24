# Script purpose
# --------------
# Map mutations to exons. Results in a maf file that has new columns starting 
# with "EVENT*". Note that mutations could map to multiple exons or to none.
# We used a MARGIN of 2 to consider mutations on splice donor or acceptor sites,
# but used the original coordinates for the result.

require(tidyverse)
require(GenomicRanges)
ROOT = here::here()
source(file.path(ROOT,'src','R','utils.R'))


# Development
# -----------
# RAW_DIR = file.path(ROOT,'data','raw')
# PREP_DIR = file.path(ROOT,'data','prep')
# input_file = file.path(RAW_DIR,'DepMap','achilles_ccle','CCLE_mutations.csv')
# annotation_file = file.path(PREP_DIR,'references','EVENT_INFO-hg38_noseqs-wmargins.tsv')
# event_type = 'EX'

##### FUNCTIONS #####
createGranges = function(chr, start, end, names){
    # 'df' with chr, start, end
    GRanges(
            seqnames = chr,
            ranges = IRanges(
                start = start,
                end = end,
                names = names
            )
        )
}

merge_exons_and_mutations = function(maf, exon_annotation, include_margins=FALSE){
    # create GRanges objects
    ## maf, (create unique index)
    maf$mutation_id = 1:nrow(maf)
    maf_ranges = createGranges(
        maf$chr, 
        maf$Start_position, 
        maf$End_position, 
        maf$mutation_id
    )
    
    ## exon_annotation
    if (include_margins){
        exon_ranges = createGranges(
            exon_annotation$chr, 
            exon_annotation$start_wmargin, 
            exon_annotation$end_wmargin, 
            exon_annotation$EVENT
        )
    }else{
        exon_ranges = createGranges(
            exon_annotation$chr, 
            exon_annotation$start, 
            exon_annotation$end, 
            exon_annotation$EVENT
        )
    }
    
    # join
    overlaps = findOverlaps(exon_ranges, maf_ranges, ignore.strand=TRUE)
    
    # prepare outputs
    ## concatenate overlaps
    mutations = maf[subjectHits(overlaps),]
    exons = exon_annotation[queryHits(overlaps),c('EVENT','chr','start','end','EVENT_type')]
    colnames(exons) = c('EVENT','EVENT_chr','EVENT_start','EVENT_end','EVENT_type')
    result = cbind(mutations, exons)
    ## add not-overlapped mutations
    not_overlapped = setdiff(maf$mutation_id, mutations$mutation_id)
    result = bind_rows(result,maf[not_overlapped,])
    
    return(result)
}


main = function(){
    args = getParsedArgs()
    
    input_file = args$input_file
    annotation_file = args$annotation_file
    output_file = args$output_file
    event_type = args$event_type
    
    # create output directory
    dir.create(dirname(output_file))
    
    # load data
    print('Loading data...')
    maf = read_csv(input_file)
    exon_annotation = read_tsv(annotation_file)
    
    print('Subsetting by event type of interest...')
    maf = maf %>% mutate(chr = paste0('chr',Chromosome))
    exon_annotation = exon_annotation %>% filter(str_detect(EVENT,event_type))
    
    # merge annotations and exons
    print('Merging event annotations with mutation coordinates...')
    merged = merge_exons_and_mutations(maf, exon_annotation) 
    
    # save
    print('Saving...')
    write_tsv(merged, output_file)
}

# #### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
