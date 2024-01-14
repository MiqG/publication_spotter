"""
Workflow purpose
----------------
Download fastq files and quantify splicing with vast-tools for Thomas 2020.
"""

import os

ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
DATA_DIR = os.path.join(ROOT,'data','raw')
SUPPORT_DIR = os.path.join(ROOT,'support')
ARTICLES_DIR = os.path.join(DATA_DIR,'articles')
SAMPLES_THOMAS2020 = ['SRR7946516','SRR7946515'] # PC9 and HeLa
VASTDB_DIR = os.path.join(DATA_DIR,'VastDB')
VASTTOOLS_DIR = "~/repositories/vast-tools/"

rule all:
    input:
        # download fastq
        os.path.join(ARTICLES_DIR,'Thomas2020','fastqs'),
        
        # align
        expand(os.path.join(ARTICLES_DIR,'Thomas2020','vast_out','{sample}'), sample=SAMPLES_THOMAS2020),
        
        # finish
        os.path.join(ARTICLES_DIR,'Thomas2020','vast_out','TPM-hg38-2.tab.gz'),
        os.path.join(ARTICLES_DIR,'Thomas2020','vast_out','INCLUSION_LEVELS_FULL-hg38-2.tab.gz'),
        os.path.join(ARTICLES_DIR,'Thomas2020','vast_out','PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz'),
        
        
        
rule download:
    input:
        sra = os.path.join(SUPPORT_DIR,'SraAccList-Thomas2020.txt')
    output:
        folder = directory(os.path.join(ARTICLES_DIR,'Thomas2020','fastqs'))
    threads: 10
    shell:
        """
        # download fastq files
        cat {input.sra} | parallel -j {threads} "fastq-dump --split-files --gzip --outdir {output.folder}"
        
        echo "Done!"
        """
        
        
rule vasttools_align:
    input:
        os.path.join(ARTICLES_DIR,'Thomas2020','fastqs'),
        dbDir = os.path.join(VASTDB_DIR,'assemblies')
    output:
        folder = directory(os.path.join(ARTICLES_DIR,'Thomas2020','vast_out','{sample}'))
    params:
        r1 = os.path.join(ARTICLES_DIR,'Thomas2020','fastqs','{sample}_1.fastq.gz'),
        bin_dir=VASTTOOLS_DIR
    threads: 10
    shell:
        """
        {params.bin_dir}/vast-tools align \
                    {params.r1} \
                    --sp Hs2 \
                    --dbDir {input.dbDir} \
                    --expr \
                    --EEJ_counts \
                    --cores {threads} \
                    --output {output}
        """
        
        
rule vasttools_combine:
    input:
        [os.path.join(ARTICLES_DIR,'Thomas2020','vast_out','{sample}').format(sample=sample) for sample in SAMPLES_THOMAS2020],
        dbDir = os.path.join(VASTDB_DIR,'assemblies')
    output:
        touch(os.path.join(ARTICLES_DIR,'Thomas2020','.done','vasttools_combine')),
        tpm = os.path.join(ARTICLES_DIR,'Thomas2020','vast_out','TPM-hg38-2.tab.gz'),
        psi = os.path.join(ARTICLES_DIR,'Thomas2020','vast_out','INCLUSION_LEVELS_FULL-hg38-2.tab.gz')
    params:
        bin_dir=VASTTOOLS_DIR,
        folder = os.path.join(ARTICLES_DIR,'Thomas2020','vast_out')
    threads: 6
    resources:
        runtime = 21600, # 6h
        memory = 10 # G
    shell:
        """
        # group results
        mkdir -p {params.folder}/to_combine
        ln -sfn {params.folder}/*/to_combine/* {params.folder}/to_combine/
        mkdir -p {params.folder}/expr_out
        ln -sfn {params.folder}/*/expr_out/* {params.folder}/expr_out/
        
        # combine runs
        {params.bin_dir}/vast-tools combine \
                    --cores {threads} \
                    --sp Hs2 \
                    --dbDir {input.dbDir} \
                    --keep_raw_reads \
                    --keep_raw_incl \
                    --output {params.folder} \
                    --TPM
        
        # compress outputs
        gzip -f {params.folder}/raw_incl/*
        gzip -f {params.folder}/raw_reads/*
        gzip -f {params.folder}/*.tab
        
        # remove grouped results
        rm -r {params.folder}/to_combine
        rm -r {params.folder}/expr_out
        
        echo "Done!"
        """

        
rule vasttools_tidy:
    input:
        os.path.join(ARTICLES_DIR,'Thomas2020','vast_out','INCLUSION_LEVELS_FULL-hg38-2.tab.gz') ## combined table
    output:
        tidy = os.path.join(ARTICLES_DIR,'Thomas2020','vast_out','PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz')
    params:
        bin_dir=VASTTOOLS_DIR
    resources:
        runtime = 21600, # 6h
        memory = 10 # G
    shell:
        """
        echo "Tidying up..."
        {params.bin_dir}/vast-tools tidy <(zcat {input}) -min_N 1 -min_SD 0 --min_ALT_use 25 --noVLOW --log -outFile {output.tidy}
        gzip --force {output.tidy}
        mv {output.tidy}.gz {output.tidy}
        
        echo "Done!"
        """
