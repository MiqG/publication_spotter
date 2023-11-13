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

CHR_OI = c(as.character(1:22),c("X","Y","MT"))

# Development
# -----------
# ROOT = here::here()
# RAW_DIR = file.path(ROOT,"data","raw")
# mapping_file = file.path(RAW_DIR,"DepMap","demeter2","shRNA-mapping.csv")
# genome_fasta_file = file.path(RAW_DIR,"GENCODE","GRCh38.p13.genome.fa.gz")
# genome_annot_file = file.path(RAW_DIR,"GENCODE","gencode.v39.annotation.gtf.gz")
# genome_annot_file = "https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz"
# genome_annot_file = file.path(RAW_DIR,"ENSEMBL","Homo_sapiens.GRCh38.110.sqlite")
# gene_oi = "RPS6KA1"
# seq_oi = "AAAAATGGCATCAACCACCAT"

##### FUNCTIONS #####
get_seq_genomic_coordinates = function(seq_oi, gene_oi, genome_annot, genome_seq){
    
    # get the exon sequences of each transcript of the gene
    ## subset gene annotations
    gene_annot = filter(genome_annot, filter=GeneNameFilter(gene_oi))
    
    ## subset chromosome sequence
    gene_chr = genes(gene_annot) %>% seqnames() %>% as.character() %>% unique() %>% sprintf("chr%s",.)
    gene_chr_seq = genome_seq[names(genome_seq)==gene_chr,][[1]]
    
    ## subset gene sequence
    gene_start = genes(gene_annot) %>% start()
    gene_end = genes(gene_annot) %>% end()
    gene_seq = gene_chr_seq[gene_start:gene_end]
    
    ## reverse complement if negative strand
    gene_strand = genes(gene_annot) %>% strand() %>% as.character() %>% unique()
    if (gene_strand=="-"){
        gene_seq = reverseComplement(gene_seq)       
    }    
    
    ## find sequence overlaps on transcripts
    seq_coords = matchPattern(seq_oi, gene_seq) %>% ranges()
    names(seq_coords) = seq_oi
    seq_coords = shift(seq_coords, shift=gene_start)
    seq_coords = GRanges(
        seqnames = gene_chr %>% gsub("chr","",.) %>% as.numeric(),
        ranges = seq_coords,
        strand = gene_strand,
        seq_oi = seq_oi,
        gene_oi = gene_oi
    )
    
    ## find overlaps with exons coordinates
    exons_annot = exons(gene_annot)
    overlaps = findOverlaps(exons_annot, seq_coords)    
    idx_exons = queryHits(overlaps)
    exons_annot = exons_annot[idx_exons,] %>% data.frame()
    
    if (length(seq_coords)>0){
        ## add corresponding genome coordinates
        seq_coords = seq_coords %>% 
            data.frame() %>%
            left_join(exons_annot, by=c("seqnames","strand"), suffix=c("","_exon_genomic"))
    } else {
        seq_coords = data.frame(
            seq_oi = seq_oi,
            gene_oi = gene_oi            
        )
    }
    
    return(seq_coords)
}


get_seq_coordinates_from_transcripts = function(seq_oi, gene_oi, genome_annot, genome_seq){
    
    # get the exon sequences of each transcript of the gene
    ## subset gene annotations
    gene_annot = filter(genome_annot, filter=GeneNameFilter(gene_oi))
    
    ## subset chromosome sequence
    gene_chr = genes(gene_annot) %>% seqnames() %>% as.character() %>% unique() %>% sprintf("chr%s",.)
    gene_chr_seq = genome_seq[names(genome_seq)==gene_chr,][[1]]
    
    transcripts_oi = exonsBy(gene_annot, by="tx", filter=~gene_name==gene_oi)
    transcripts_seqs = extractTranscriptSeqs(gene_chr_seq, ranges(transcripts_oi))
    
    ## reverse complement if negative strand
    gene_strand = genes(gene_annot) %>% strand() %>% as.character() %>% unique()
    if (gene_strand=="-"){
        transcripts_seqs = reverseComplement(transcripts_seqs)       
    }
    
    ## find sequence overlaps on transcripts
    seq_coords_transcripts = vmatchPattern(
        seq_oi, transcripts_seqs
    ) %>% unlist()
    
    ## exons coordinates
    exons_annot = exons(gene_annot) %>% data.frame()
    
    if (length(seq_coords_transcripts)>0){
        ## add corresponding genome coordinates
        seq_coords = transcriptToGenome(seq_coords_transcripts, gene_annot) %>% 
            unlist() %>%
            data.frame() %>%
            mutate(
                seq_oi = seq_oi,
                gene_oi = gene_oi
            ) %>%
            left_join(exons_annot, by=c("exon_id","seqnames","strand"), suffix=c("","_exon_genomic"))
    } else {
        seq_coords = data.frame(
            seq_oi = seq_oi,
            gene_oi = gene_oi            
        )
    }
    
    return(seq_coords)
}

parseargs = function(){
    
    option_list = list( 
        make_option("--mapping_file", type="character"),
        make_option("--genome_annot_file", type="character"),
        make_option("--genome_fasta_file", type="character"),
        make_option("--output_file", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}

main = function(){
    args = parseargs()
    
    mapping_file = args[["mapping_file"]]
    genome_annot_file = args[["genome_annot_file"]]
    genome_fasta_file = args[["genome_fasta_file"]]
    output_file = args[["output_file"]]
    
    # load
    mapping = read_csv(mapping_file)
    genome_seq = readDNAStringSet(genome_fasta_file)
    genome_annot = EnsDb(genome_annot_file)
    
    # prep
    ## correct chromosome numbers
    names(genome_seq) = gsub(" .*", "", names(genome_seq))
    ## subset genome annotations
    genome_annot = filter(genome_annot, filter = SeqNameFilter(CHR_OI))
    
    ## subset
    common_genes = intersect(genes(genome_annot)$gene_name, mapping[["Gene Symbol"]])
    genome_annot = filter(genome_annot, filter = GeneNameFilter(common_genes))
    mapping = mapping %>% dplyr::filter(`Gene Symbol` %in% common_genes)
    
    # get coordinates
    coordinates_transcripts = lapply(
        #1:nrow(mapping), function(idx){
        1:10, function(idx){
            
        seq_oi = mapping[idx,"Barcode Sequence"][[1]]
        gene_oi = mapping[idx,"Gene Symbol"][[1]]
        
        seq_coords = tryCatch({
            
            seq_coords = get_seq_coordinates_from_transcripts(
                seq_oi,
                gene_oi,
                genome_annot,
                genome_seq
            )
            return(seq_coords)
            
        }, error= function(cond) {
            
            seq_coords = data.frame(
                seq_oi = seq_oi,
                gene_oi = gene_oi            
            )
            print("Error with:"); print(seq_coords)
            return(seq_coords)
            
        })
        
        return(seq_coords)
        
    }) %>% bind_rows()
    
    coordinates = lapply(
        #1:nrow(mapping), function(idx){
        1:10, function(idx){
            
        seq_oi = mapping[idx,"Barcode Sequence"][[1]]
        gene_oi = mapping[idx,"Gene Symbol"][[1]]
        
        seq_coords = tryCatch({
            
            seq_coords = get_seq_genomic_coordinates(
                seq_oi,
                gene_oi,
                genome_annot,
                genome_seq
            )
            return(seq_coords)
            
        }, error= function(cond) {
            
            seq_coords = data.frame(
                seq_oi = seq_oi,
                gene_oi = gene_oi            
            )
            print("Error with:"); print(seq_coords)
            return(seq_coords)
            
        })
        
        return(seq_coords)
        
    }) %>% bind_rows()
    
    # save
    write_tsv(coordinates, output_file)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}