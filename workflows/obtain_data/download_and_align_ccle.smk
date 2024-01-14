"""
Workflow purpose
----------------
Download fastq files and quantify splicing with vast-tools for CCLE.
"""

import os
import pandas as pd
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
DATA_DIR = os.path.join(ROOT,'data','raw')
SUPPORT_DIR = os.path.join(ROOT,'support')
VASTDB_DIR = os.path.join(DATA_DIR,"VastDB")
CCLE_DIR = os.path.join(DATA_DIR,'CCLE')
VASTTOOLS_DIR = "~/repositories/vast-tools/"

# prep metadata
SAMPLES_CCLE = pd.read_table(os.path.join(SUPPORT_DIR,'SraAccList-CCLE.txt'),header=None)[0]
metadata = pd.read_table(os.path.join(SUPPORT_DIR,'ENA_filereport-PRJNA523380-CCLE.tsv'))
metadata = metadata.loc[metadata["library_source"]=="TRANSCRIPTOMIC"]
metadata = metadata.loc[metadata['run_accession'].isin(SAMPLES_CCLE)]
metadata["size"] = metadata["fastq_bytes"].astype(str).str.split(";").apply(lambda x: max(np.array(x, dtype=int)))
metadata = metadata.sort_values("size", ascending=False)

# URLs to download
URLS_CCLE = metadata["fastq_ftp"].str.split(';').str[0].apply(os.path.dirname).to_list()
URLS_CCLE = {os.path.basename(url): url for url in URLS_CCLE}

SAMPLES_CCLE = list(URLS_CCLE.keys())

# read ends
ENDS = ["1","2"]

# fastq sizes
SIZES = metadata.set_index("run_accession")["size"].to_dict()
SIZE_THRESH = 5e9

rule all:
    input:        
        # download
        expand(os.path.join(CCLE_DIR,'fastqs','{sample}_{end}.fastq.gz'), sample=SAMPLES_CCLE, end=ENDS),
        
        # quantify splicing and gene expression
        ## Download .fastq files and Quantify splicing and expression
        expand(os.path.join(CCLE_DIR,'fastqs','.done','{sample}_{end}'), end=ENDS, sample=SAMPLES_CCLE),
        expand(os.path.join(CCLE_DIR,'vast_out','.done','{sample}'), sample=SAMPLES_CCLE),
        ## Combine into single tables
        os.path.join(CCLE_DIR,'vast_out','.done/vasttools_combine'),
        
        ## get raw read counts
        os.path.join(CCLE_DIR,'vast_out','inclusion_reads.tsv.gz'),
        os.path.join(CCLE_DIR,'vast_out','exclusion_reads.tsv.gz'),
        os.path.join(CCLE_DIR,'vast_out','total_reads.tsv.gz'),
        
        ## Tidy PSI
        '.done/CCLE.done',
        os.path.join(CCLE_DIR,'vast_out','PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz'),
        ## Copy annotations
        os.path.join(CCLE_DIR,'ENA_filereport-PRJNA523380-CCLE.tsv'),
        
        
##### CCLE #####
rule copy_sample_identifiers:
    input:
        os.path.join('support','ENA_filereport-PRJNA523380-CCLE.tsv')
    output:
        os.path.join(CCLE_DIR,'ENA_filereport-PRJNA523380-CCLE.tsv')
    shell:
        """
        cp {input} {output}
        """
        
        
rule download:
    params:
        sample = "{sample}",
        end = "{end}",
        url = lambda wildcards: URLS_CCLE[wildcards.sample],
        fastqs_dir = os.path.join(CCLE_DIR,"fastqs"),
        bin_dir=VASTTOOLS_DIR
    output:
        os.path.join(CCLE_DIR,'fastqs','{sample}_{end}.fastq.gz'),
        download_done = touch(os.path.join(CCLE_DIR,'fastqs','.done','{sample}_{end}'))
    threads: 3
    resources:
        runtime = 3600*2, # 2 h
        memory = 1 # GB
    shell:
        """
        set -eo pipefail

        # download
        echo "Downloading {params.sample}_{params.end}..."
        nice axel --num-connections={threads} \
             --output={params.fastqs_dir}/{params.sample}_{params.end}.fastq.gz \
             {params.url}/{params.sample}_{params.end}.fastq.gz
        
        echo "Finished downloading {params.sample}_{params.end}."
        echo $(date)
        
        echo "Done!"
        """
        
        
rule quantify_psi_vasttools:
    params:
        sample = '{sample}',
        fastqs_dir = os.path.join(CCLE_DIR,'fastqs'),
        bin_dir=VASTTOOLS_DIR,
        vast_out = directory(os.path.join(CCLE_DIR,'vast_out','{sample}'))
    input:
        dbDir = os.path.join(VASTDB_DIR,'assemblies'),
        download_done = [os.path.join(CCLE_DIR,'fastqs','.done','{sample}_{end}').format(end=end, sample='{sample}') for end in ENDS]
    output:
        align_done = touch(os.path.join(CCLE_DIR,'vast_out','.done','{sample}'))
    threads: 16
    resources:
        runtime = 86400, # most 6h is enough; some needed 24h (more reads).
        memory = 15
    group: "CCLE"
    shell:
        """
        set -eo pipefail

        # align paired reads
        echo "Aligning {params.sample}..."
        {params.bin_dir}/vast-tools align \
                    {params.fastqs_dir}/{params.sample}_1.fastq.gz \
                    {params.fastqs_dir}/{params.sample}_2.fastq.gz \
                    --sp Hs2 \
                    --dbDir {input.dbDir} \
                    --expr \
                    --EEJ_counts \
                    --cores {threads} \
                    --output {params.vast_out}
                    
        echo "Finished aligning {params.sample}."
        echo $(date)
        
        echo "Done!"
        """

        
rule vasttools_combine_CCLE:
    input:
        [os.path.join(CCLE_DIR,'vast_out','{sample}').format(sample=sample) for sample in SAMPLES_CCLE],
        dbDir = os.path.join(VASTDB_DIR,'assemblies')
    output:
        touch(os.path.join(CCLE_DIR,'vast_out','.done','vasttools_combine')),
        tpm = os.path.join(CCLE_DIR,'vast_out','TPM-hg38-1019.tab.gz'),
        psi = os.path.join(CCLE_DIR,'vast_out','INCLUSION_LEVELS_FULL-hg38-1019.tab.gz')
    params:
        bin_dir=VASTTOOLS_DIR,
        folder = os.path.join(CCLE_DIR,'vast_out')
    threads: 56
    resources:
        runtime = 172800, # 48h
        memory = 200
    shell:
        """
        set -eo pipefail
        
        # group results
        echo "Grouping results..."
        mkdir -p {params.folder}/to_combine
        ln -s {params.folder}/*/to_combine/* {params.folder}/to_combine/
        mkdir -p {params.folder}/expr_out
        ln -s {params.folder}/*/expr_out/* {params.folder}/expr_out/
        
        # combine runs
        echo "Combining runs..."
        {params.bin_dir}/vast-tools combine \
                    --cores {threads} \
                    --sp Hs2 \
                    --dbDir {input.dbDir} \
                    --keep_raw_reads \
                    --keep_raw_incl \
                    --output {params.folder} \
                    --TPM \
                    -C
        
        # compress outputs
        echo "Compressing outputs..."
        gzip -f {params.folder}/raw_incl/*
        gzip -f {params.folder}/raw_reads/*
        gzip -f {params.folder}/*.tab
        
        # remove grouped results
        echo "Removing grouped results..."
        rm -r {params.folder}/to_combine
        rm -r {params.folder}/expr_out
        
        echo "Done!"
        """
    
    
rule vasttools_tidy_CCLE:
    input:
        os.path.join(CCLE_DIR,'vast_out','INCLUSION_LEVELS_FULL-hg38-1019.tab.gz') ## combined table
    output:
        '.done/CCLE.done',
        tidy = os.path.join(CCLE_DIR,'vast_out','PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz')
    params:
        bin_dir=VASTTOOLS_DIR
    threads: 1
    resources:
        runtime = 43200, # 12h
        memory = 60
    shell:
        """
        echo "Tidying up..."
        {params.bin_dir}/vast-tools tidy <(zcat {input}) -min_N 1 -min_SD 0 --min_ALT_use 25 --noVLOW --log -outFile {output.tidy}
        gzip --force {output.tidy}
        mv {output.tidy}.gz {output.tidy}
        
        touch .done/CCLE.done
        
        echo "Done!"
        """

        
rule get_raw_reads:
    input:
        inclusion_levels = os.path.join(CCLE_DIR,'vast_out','INCLUSION_LEVELS_FULL-hg38-1019.tab.gz')
    output:
        inclusion_reads = os.path.join(CCLE_DIR,'vast_out','inclusion_reads.tsv.gz'),
        exclusion_reads = os.path.join(CCLE_DIR,'vast_out','exclusion_reads.tsv.gz'),
        total_reads = os.path.join(CCLE_DIR,'vast_out','total_reads.tsv.gz')
    params:
        chunksize = 10000
    threads: 24
    resources:
        runtime = 6*3600, # 6h
        memory = 20 # GB
    shell:
        """
        nice python scripts/get_raw_reads.py \
                    --inclusion_levels_file={input.inclusion_levels} \
                    --n_jobs={threads} \
                    --chunksize={params.chunksize} \
                    --inclusion_reads_file={output.inclusion_reads} \
                    --exclusion_reads_file={output.exclusion_reads} \
                    --total_reads_file={output.total_reads}
        """
