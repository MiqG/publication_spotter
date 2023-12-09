#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Prepare metadata from thomas 2020

require(optparse)
require(tidyverse)
require(ensembldb)
require(Biostrings)
require(foreach)
require(doParallel)

CHR_OI = c(as.character(1:22),c("X","Y","MT"))

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,"data","raw")
# PREP_DIR = file.path(ROOT,"data","prep")
# mapping_file = file.path(RAW_DIR,"DepMap","demeter2","shRNA-mapping.csv")
# alignment_file = file.path(PREP_DIR,"demeter2","shRNA_to_gencode.v44.transcripts.bed")
# genome_fasta_file = file.path(RAW_DIR,"GENCODE","GRCh38.p13.genome.fa.gz")
# genome_annot_file = file.path(RAW_DIR,"GENCODE","gencode.v39.annotation.gtf.gz")
# genome_annot_file = "https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz"
# genome_annot_file = file.path(RAW_DIR,"ENSEMBL","Homo_sapiens.GRCh38.110.sqlite")
# genome_annot_file = file.path(RAW_DIR,"GENCODE","Homo_sapiens.GRCh38.gencode_v44.sqlite")
# gene_oi = "RPS6KA1"
# seq_oi = "AAAAATGGCATCAACCACCAT"

##### FUNCTIONS #####
parseargs = function(){
    
    option_list = list( 
        make_option("--alignment_file", type="character"),
        make_option("--genome_annot_file", type="character"),
        make_option("--chunk_size", type="integer"),
        make_option("--n_jobs", type="integer"),
        make_option("--output_file", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    alignment_file = args[["alignment_file"]]
    genome_annot_file = args[["genome_annot_file"]]
    chunk_size = args[["chunk_size"]]
    n_jobs = args[["n_jobs"]]
    output_file = args[["output_file"]]
    
    # load
    alignment = read_tsv(alignment_file, col_names = c("transcript_id","start","end","shrna_id","rgb","strand"))
    genome_annot = EnsDb(genome_annot_file)

    # prep
    ## alignment
    alignment = alignment %>% 
        separate(
            transcript_id, 
            into=c("transcript_id","gene_id","havana_transcript","havana_gene",
                   "transcript_name","gene_name","entrez_id","gene_type","empty"),
            sep="\\|"
        ) %>%
        separate(shrna_id, into=c("barcode_sequence","gene_name_prevmap"), sep="\\|") %>%
        mutate(
            tx_id = gsub("\\..*","",transcript_id),
            tx_id = transcript_id,
            alignment_index = row_number()
        ) %>%
        dplyr::select(-strand) 
    
    ## add transcript info
    transcripts_oi = alignment %>% pull(tx_id) %>% unique()
    transcripts_annot = filter(genome_annot, filter = TxIdFilter(transcripts_oi))
    transcripts_info = transcripts(transcripts_annot) %>% as.data.frame()

    alignment = alignment %>%
        left_join(
            transcripts_info %>% distinct(tx_id, seqnames, strand),
            by="tx_id"
        )

    alignment_ranges = GRanges(
        seqnames = alignment[["seqnames"]],
        ranges = IRanges(
            start = alignment[["start"]],
            end = alignment[["end"]],
            name = alignment[["tx_id"]]
        ),
        strand = alignment[["strand"]],
        tx_id = alignment[["tx_id"]],
        shrna_barcode_sequence = alignment[["barcode_sequence"]],
        shrna_gene_name = alignment[["gene_name"]],
        alignment_index = alignment[["alignment_index"]]
    )
    
    # map alignments to genome annotations
    print("Mapping alignments to genomic features...")
    start_time = Sys.time()
    
    registerDoParallel(cores=n_jobs)
    
    #chunks = split(1:length(alignment_ranges[1:100]), ceiling(1:length(alignment_ranges[1:100]) / chunk_size))
    chunks = split(1:length(alignment_ranges), ceiling(1:length(alignment_ranges) / chunk_size))
    
    result = foreach(chunk = chunks) %dopar% {
        # open connection
        genome_annot = EnsDb(genome_annot_file)
        transcripts_annot = filter(genome_annot, filter = TxIdFilter(transcripts_oi))

        # mapping
        mapping = transcriptToGenome(ranges(alignment_ranges[chunk]), transcripts_annot)
        mapping = lapply(1:length(mapping), function(i){
                mapping[[i]] %>% data.frame() %>% mutate(alignment_index = chunk[i])
            }) %>% 
            bind_rows() %>%
            dplyr::rename(genomic_start = start, genomic_end = end) %>%
            left_join(alignment, by=c("alignment_index", "tx_id", "seqnames", "strand"))

        # close connection
        on.exit(RSQLite::dbDisconnect(dbconn(genome_annot)))
        gc()
        return(mapping)
    } %>% bind_rows()
    
    end_time = Sys.time()
    execution_time <- end_time - start_time
    print(execution_time)
    
    # add exon info
    exons_info = exons(transcripts_annot) %>% data.frame() %>% distinct()
    result = result %>%
        left_join(exons_info, by=c("exon_id","seqnames","strand"), suffix=c("","_exon"))
    
    # save
    write_tsv(result, output_file)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}