"""
Workflow purpose
----------------
Download fastq files and quantify splicing with vast-tools for Hart 2015.
"""

import os
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
SUPPORT_DIR = os.path.join(ROOT,'support')
DATA_DIR = os.path.join(ROOT,'data','raw')
ARTICLES_DIR = os.path.join(DATA_DIR,'articles')
ARTICLE_DIR = os.path.join(ARTICLES_DIR,'Hart2015')
VASTDB_DIR = os.path.join(DATA_DIR,'VastDB')
VASTTOOLS_DIR = "~/repositories/vast-tools/"

# load metadata
metadata = pd.read_table(os.path.join(SUPPORT_DIR,'ENA_filereport-PRJNA302748-Hart2015.tsv'))
metadata = metadata.loc[metadata["library_source"]=="TRANSCRIPTOMIC"]

# prepare
## URLs to download
URLS = metadata['fastq_ftp'].str.split(';').str[0].apply(os.path.dirname).to_list()
URLS = {os.path.basename(url): url for url in URLS}
SAMPLES = list(URLS.keys())

## fastq sizes
SIZES = metadata.set_index("run_accession")["fastq_bytes"].to_dict()
SIZE_THRESH = 5e9

rule all:
    input:
        # Download .fastq files 
        expand(os.path.join(ARTICLE_DIR,'fastqs','.done','{sample}'), sample=SAMPLES),
        
        # Quantify splicing and expression
        expand(os.path.join(ARTICLE_DIR,'vast_out','.done','{sample}'), sample=SAMPLES),

        # Combine into single tables
        os.path.join(ARTICLE_DIR,'vast_out','TPM-hg38-12.tab.gz'),
        os.path.join(ARTICLE_DIR,'vast_out','INCLUSION_LEVELS_FULL-hg38-12.tab.gz'),
        
        # Tidy PSI
        os.path.join(ARTICLE_DIR,'vast_out','PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz')

        
rule download:
    params:
        sample = '{sample}',
        url = lambda wildcards: URLS[wildcards.sample],
        fastqs_dir = os.path.join(ARTICLE_DIR,'fastqs'),
        bin_dir=VASTTOOLS_DIR
    output:
        download_done = os.path.join(ARTICLE_DIR,'fastqs','.done','{sample}')
    threads: 10
    resources:
        runtime = 7200, # 2h
        memory = 2
    group: "Hart2015"
    shell:
        """
        set -eo pipefail
        
        # download
        echo "Downloading {params.sample}..."
             
        nice axel --num-connections={threads} \
             {params.url}/{params.sample}.fastq.gz \
             --output={params.fastqs_dir}/{params.sample}.fastq.gz
        
        touch {output.download_done}
        echo "Finished downloading {params.sample}."
        echo $(date)
        
        echo "Done!"
        """
        

rule align:
    params:
        sample = '{sample}',
        fastqs_dir = os.path.join(ARTICLE_DIR,'fastqs'),
        bin_dir=VASTTOOLS_DIR,
        vast_out = directory(os.path.join(ARTICLE_DIR,'vast_out','{sample}'))
    input:
        dbDir = os.path.join(VASTDB_DIR,'assemblies'),
        download_done = os.path.join(ARTICLE_DIR,'fastqs','.done','{sample}')
    output:
        align_done = touch(os.path.join(ARTICLE_DIR,'vast_out','.done','{sample}'))
    threads: 16
    resources:
        runtime = lambda wildcards: 86400 if SIZES[wildcards.sample]>SIZE_THRESH else 21600, # most 6h is enough; some needed 24h (more reads).
        memory = 15
    group: "Hart2015"
    shell:
        """
        set -eo pipefail
        
        # align paired reads
        echo "Aligning {params.sample}..."
        echo $(date)
        
        {params.bin_dir}/vast-tools align \
                    {params.fastqs_dir}/{params.sample}.fastq.gz \
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

        
rule vasttools_combine:
    input:
        [os.path.join(ARTICLE_DIR,'vast_out','.done','{sample}').format(sample=sample) for sample in SAMPLES],
        dbDir = os.path.join(VASTDB_DIR,'assemblies')
    output:
        touch(os.path.join(ARTICLE_DIR,'vast_out','.done','vasttools_combine')),
        tpm = os.path.join(ARTICLE_DIR,'vast_out','TPM-hg38-{n_samples}.tab.gz').format(n_samples=len(SAMPLES)),
        psi = os.path.join(ARTICLE_DIR,'vast_out','INCLUSION_LEVELS_FULL-hg38-{n_samples}.tab.gz').format(n_samples=len(SAMPLES))
    params:
        bin_dir=VASTTOOLS_DIR,
        folder = os.path.join(ARTICLE_DIR,'vast_out')
    threads: 16
    resources:
        runtime = 24*3600, # 24h
        memory = 60
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
                    --TPM

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
    
    
rule vasttools_tidy:
    input:
        os.path.join(ARTICLE_DIR,'vast_out','INCLUSION_LEVELS_FULL-hg38-{n_samples}.tab.gz').format(n_samples=len(SAMPLES)) ## combined table
    output:
        touch('.done/Hart2015.done'),
        tidy = os.path.join(ARTICLE_DIR,'vast_out','PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz')
    params:
        bin_dir=VASTTOOLS_DIR
    threads: 1
    resources:
        runtime = 12*3600, # 12h
        memory = 20
    shell:
        """
        set -eo pipefail

        echo "Tidying up..."
        {params.bin_dir}/vast-tools tidy <(zcat {input}) -min_N 1 -min_SD 0 --min_ALT_use 25 --noVLOW --log -outFile {output.tidy}
        gzip --force {output.tidy}
        mv {output.tidy}.gz {output.tidy}
        
        echo "Done!"
        """