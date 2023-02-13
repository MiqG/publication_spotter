"""
Workflow purpose
----------------
Download fastq files and quantify splicing with vast-tools for TCGA.
Also, download gene expression counts from GDC.
"""

import os
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
SUPPORT_DIR = os.path.join(ROOT,'support')
DATA_DIR = os.path.join(ROOT,'data','raw')
ARTICLES_DIR = os.path.join(DATA_DIR,'articles')
ARTICLE_DIR = os.path.join(ARTICLES_DIR,'Pietilla2021')
VASTDB_DIR = os.path.join(DATA_DIR,'VastDB')
TOKEN_FILE = os.path.join(SUPPORT_DIR,".private","ega_aspera_pass.txt")
DECRYPTION_KEY_FILE = os.path.join(SUPPORT_DIR,".private","ega_decryption_key")
DECRYPTION_PASS_FILE = os.path.join(SUPPORT_DIR,".private","ega_decryption_pass")
METADATA_FILE = os.path.join(SUPPORT_DIR,".private","ega_files-EGAD00001006456.txt")

# load metadata
metadata = pd.read_table(METADATA_FILE, header=None).rename(columns={0:"file_id"})
metadata["sampleID"] = metadata["file_id"].apply(os.path.basename).str.replace("_[1-2].fq.gz.c4gh", "", regex=True)

# prepare
## URLs to download
SAMPLES = list(metadata["sampleID"].unique())
URLS = {s: metadata.loc[metadata["sampleID"]==s,"file_id"].apply(os.path.dirname).values[0] for s in SAMPLES}
## paired ends
ENDS = ["1","2"]

rule all:
    input:
        # Download .fastq files and Quantify splicing and expression
        expand(os.path.join(ARTICLE_DIR,'fastqs','.done','{sample}_{end}'), end=ENDS, sample=SAMPLES),
        expand(os.path.join(ARTICLE_DIR,'fastqs','.done_decryption','{sample}_{end}'), end=ENDS, sample=SAMPLES),
        expand(os.path.join(ARTICLE_DIR,'vast_out','.done','{sample}'), sample=SAMPLES),

        # Combine into single tables
        expand(os.path.join(ARTICLE_DIR,"vast_out",".done_combine-{n_samples}"), n_samples=len(SAMPLES)),
        
        # Tidy PSI
        '.done/Pietilla2021.done',
        os.path.join(ARTICLE_DIR,'vast_out','PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz'),

        
rule download:
    params:
        sample = "{sample}",
        url = lambda wildcards: URLS[wildcards.sample],
        end = "{end}",
        fastqs_dir = os.path.join(ARTICLE_DIR,'fastqs'),
        bin_dir="~/repositories/aspera/connect/bin"
    input:
        token = TOKEN_FILE
    output:
        download_done = os.path.join(ARTICLE_DIR,'fastqs','.done','{sample}_{end}')
    threads: 1
    resources:
        runtime = 7200, # 2h
        memory = 2
    group: "Pietilla2021"
    shell:
        """
        set -eo pipefail
       
        echo "Downloading {params.sample}_{params.end}..."
        mkdir -p {params.fastqs_dir}
        
        export ASPERA_SCP_PASS=$(cat {input.token})
        
        nice {params.bin_dir}/ascp \
            --ignore-host-key \
            -k 1 \
            --partial-file-suffix=PART \
            -QT \
            dbox3021@xfer.ega-archive.org:{params.url}/{params.sample}_{params.end}.fq.gz.c4gh \
            {params.fastqs_dir}/{params.sample}_{params.end}.fq.gz.c4gh 
        
        touch {output.download_done}
        echo "Finished downloading {params.sample}."
        echo $(date)
        
        echo "Done!"        
        """
        
        
rule decrypt:
    params:
        sample = "{sample}",
        end = "{end}",
        fastqs_dir = os.path.join(ARTICLE_DIR,'fastqs')
    input:
        decryption_key = DECRYPTION_KEY_FILE,
        decryption_pass = DECRYPTION_PASS_FILE
    output:
        decryption_done = os.path.join(ARTICLE_DIR,'fastqs','.done_decryption','{sample}_{end}')
    threads: 1
    resources:
        runtime = 2*3600, # 2h
        memory = 2
    shell:
        """
        set -eo pipefail
        
        echo "Decrypting {params.sample}_{params.end}..."
        echo $(date)
        
        export C4GH_SECRET_KEY={input.decryption_key}
        export C4GH_PASSPHRASE=$(cat {input.decryption_pass})
        
        nice cat {params.fastqs_dir}/{params.sample}_{params.end}.fq.gz.c4gh |\
            crypt4gh decrypt > \
            {params.fastqs_dir}/{params.sample}_{params.end}.fq.gz
        
        touch {output.decryption_done}
        echo "Finished downloading {params.sample}."
        echo $(date)
        
        echo "Done!"
        """
        
        
rule align:
    params:
        sample = '{sample}',
        fastqs_dir = os.path.join(ARTICLE_DIR,'fastqs'),
        bin_dir="~/repositories/vast-tools/",
        vast_out = directory(os.path.join(ARTICLE_DIR,'vast_out','{sample}'))
    input:
        dbDir = os.path.join(VASTDB_DIR,'assemblies'),
        download_done = [os.path.join(ARTICLE_DIR,'fastqs','.done','{sample}_{end}').format(end=end, sample='{sample}') for end in ENDS],
        decryption_done = [os.path.join(ARTICLE_DIR,'fastqs','.done_decryption','{sample}_{end}').format(end=end, sample='{sample}') for end in ENDS]
    output:
        align_done = touch(os.path.join(ARTICLE_DIR,'vast_out','.done','{sample}'))
    threads: 16
    resources:
        runtime = 6*3600, # 6h
        memory = 15
    group: "Pietilla2021"
    shell:
        """
        set -eo pipefail
        
        echo "Aligning {params.sample}..."
        echo $(date)

        nice {params.bin_dir}/vast-tools align \
                    {params.fastqs_dir}/{params.sample}_1.fq.gz \
                    {params.fastqs_dir}/{params.sample}_2.fq.gz \
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
        touch(os.path.join(ARTICLE_DIR,'vast_out','.done_combine-{n_samples}')),
        tpm = os.path.join(ARTICLE_DIR,'vast_out','TPM-hg38-{n_samples}.tab.gz'),
        psi = os.path.join(ARTICLE_DIR,'vast_out','INCLUSION_LEVELS_FULL-hg38-{n_samples}.tab.gz')
    params:
        bin_dir="~/repositories/vast-tools/",
        folder = os.path.join(ARTICLE_DIR,'vast_out')
    threads: 16
    resources:
        runtime = 24*3600, # 24h
        memory = 90
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
        touch('.done/Pietilla2021.done'),
        tidy = os.path.join(ARTICLE_DIR,'vast_out','PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz')
    params:
        bin_dir="~/repositories/vast-tools/"
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
