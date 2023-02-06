"""
Workflow purpose
----------------
Download fastq files and quantify splicing with vast-tools for ENCORE.
"""

import os
import pandas as pd
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
DATA_DIR = os.path.join(ROOT,'data','raw')
SUPPORT_DIR = os.path.join(ROOT,'support')
ENCODE_DIR = os.path.join(DATA_DIR,'ENCODE') # parent
ENCORE_DIR = os.path.join(ENCODE_DIR,'ENCORE') # child
VASTDB_DIR = os.path.join(DATA_DIR,'VastDB')
VASTTOOLS_DIR = "~/repositories/vast-tools/"

ENDS = ["1","2"]

# load metadata (previously downloaded)
metadata = pd.read_table(os.path.join(SUPPORT_DIR,'ENCORE-metadata.tsv'))
metadata["sample"] = metadata["dbxrefs"].str.replace("SRA:","")
metadata["end"] = metadata["Paired end"].astype('str')

# dictionary of paired end fastq files
SAMPLES_ENCORE = {}
for sample_oi in metadata["sample"].unique():
    SAMPLES_ENCORE[sample_oi] = metadata.loc[
        metadata["sample"]==sample_oi,
        ["end","File accession"]]\
        .set_index("end")\
        .to_dict()["File accession"]
    
URLS = metadata.set_index("File accession")["File download URL"].to_dict()
SIZES = metadata.set_index("File accession")["Size"].to_dict()
SIZE_THRESH = 5e9

rule all:
    input:
        # download metadata
        os.path.join(ENCORE_DIR,'metadata','ENCORE.tsv'),
        
        # Download .fastq files
        expand(os.path.join(ENCORE_DIR,'fastqs','.done','{sample}_{end}'), end=ENDS, sample=SAMPLES_ENCORE.keys()),
        
        # quantify splicing
        expand(os.path.join(ENCORE_DIR,'vast_out','.done','{sample}'), sample=SAMPLES_ENCORE.keys()),

        # Combine into single tables
        os.path.join(ENCORE_DIR,'vast_out','.done/vasttools_combine'),
        
        # Tidy PSI
        ".done/ENCORE.done",
        os.path.join(ENCORE_DIR,'vast_out','PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz')
        
        
rule download_metadata:
    params:
        metadata = "'https://www.encodeproject.org/metadata/?assay_title=shRNA+RNA-seq&files.file_type=fastq&status=released&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&searchTerm=shrna&type=Experiment'"
    output:
        metadata = os.path.join(ENCORE_DIR,'metadata','ENCORE.tsv'),
        readme = os.path.join(ENCORE_DIR,'metadata','README.md')
    shell:
        """
        # metadata
        wget --user-agent="Chrome" --no-clobber --no-check-certificate {params.metadata} -O {output.metadata}
        
        # readme
        echo "Downloaded on $(date)." > {output.readme}
        
        echo Done!
        """
        
        
rule download:
    params:
        sample = '{sample}',
        end = '{end}',
        url = lambda wildcards: URLS[SAMPLES_ENCORE[wildcards.sample][wildcards.end]],
        fastqs_dir = os.path.join(ENCORE_DIR,'fastqs'),
    output:
        download_done = os.path.join(ENCORE_DIR,'fastqs','.done','{sample}_{end}')
    threads: 1
    resources:
        runtime = 5400,
        memory = 2
    group: "ENCORE"
    shell:
        """
        # download
        echo "Downloading {params.sample}..."
        
        wget --user-agent="Chrome" --no-check-certificate {params.url} -O {params.fastqs_dir}/{params.sample}_{params.end}.fastq.gz
        
        touch {output.download_done}
        echo "Finished downloading {params.sample}."
        echo $(date)
        
        echo "Done!"
        """
        
        
rule vasttools_align:
    params:
        sample = '{sample}',
        fastqs_dir = os.path.join(ENCORE_DIR,'fastqs'),
        bin_dir=VASTTOOLS_DIR,
        vast_out = directory(os.path.join(ENCORE_DIR,'vast_out','{sample}'))
    input:
        dbDir = os.path.join(VASTDB_DIR,'assemblies'),
        download_done = [os.path.join(ENCORE_DIR,'fastqs','.done','{sample}_{end}').format(end=end, sample='{sample}') for end in ENDS]
    output:
        align_done = touch(os.path.join(ENCORE_DIR,'vast_out','.done','{sample}'))
    threads: 16
    resources:
        runtime = lambda wildcards: 86400 if any([SIZES[SAMPLES_ENCORE[wildcards.sample]["1"]]>SIZE_THRESH,SIZES[SAMPLES_ENCORE[wildcards.sample]["2"]]>SIZE_THRESH]) else 21600, # most 6h is enough; some needed 24h (more reads).
        memory = 15
    group: "ENCORE"
    shell:
        """
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

        
rule vasttools_combine:
    input:
        [os.path.join(ENCORE_DIR,'vast_out','{sample}').format(sample=sample) for sample in SAMPLES_ENCORE.keys()],
        dbDir = os.path.join(VASTDB_DIR,'assemblies')
    output:
        touch(os.path.join(ENCORE_DIR,'vast_out','.done','vasttools_combine')),
        tpm = os.path.join(ENCORE_DIR,'vast_out','TPM-hg38-1097.tab.gz'),
        psi = os.path.join(ENCORE_DIR,'vast_out','INCLUSION_LEVELS_FULL-hg38-1097.tab.gz')
    params:
        bin_dir=VASTTOOLS_DIR,
        folder = os.path.join(ENCORE_DIR,'vast_out')
    threads: 16
    resources:
        runtime = 172800, # 48h
        memory = 60
    shell:
        """
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
        os.path.join(ENCORE_DIR,'vast_out','INCLUSION_LEVELS_FULL-hg38-1097.tab.gz') ## combined table
    output:
        touch('.done/ENCORE.done'),
        tidy = os.path.join(ENCORE_DIR,'vast_out','PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz')
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
        
        echo "Done!"
        """