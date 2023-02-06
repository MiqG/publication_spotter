"""
Workflow purpose
----------------
Download fastq files and quantify splicing with vast-tools for TCGA.
Also, download gene expression counts from GDC.
"""

import os
import pandas as pd
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
DATA_DIR = os.path.join(ROOT,'data','raw')
SUPPORT_DIR = os.path.join(ROOT,'support')
TCGA_DIR = os.path.join(DATA_DIR,'TCGA') 
GDC_DIR = os.path.join(DATA_DIR,"GDC")
VASTDB_DIR = os.path.join(DATA_DIR,'VastDB')
VASTTOOLS_DIR = "~/repositories/vast-tools/"
TOKEN_FILE = os.path.join(SUPPORT_DIR,".private","gdc-user-token.txt")

# variables
CANCER_TYPES = [
    "BLCA",
    "BRCA",
    "COAD",
    "HNSC",
    "KICH",
    "KIRC",
    "KIRP",
    "LIHC",
    "LUAD",
    "LUSC",
    "PRAD",
    "THCA",
    "UCEC"
]

# load metadata
metadata_json = pd.concat([
    pd.read_json(os.path.join(SUPPORT_DIR,"TCGA-metadata.cart.2022-07-22-brca.json")).assign(sequencer="IlluminaHiSeq"),
    pd.read_json(os.path.join(SUPPORT_DIR,"TCGA-metadata.cart.2022-07-22-rest.json")).assign(sequencer="IlluminaHiSeq"),
    pd.read_json(os.path.join(SUPPORT_DIR,"TCGA-metadata.cart.2022-12-27-rest_illumina_ga.json")).assign(sequencer="IlluminaGA")
])
metadata = []
for idx, row in metadata_json.iterrows():
        metadata.append({
            "file_name": row["file_name"],
            "file_id": row["file_id"],
            "file_size": row["file_size"],
            "project_id": row["cases"][0]["project"]["project_id"],
            "sequencer": row["sequencer"]
        })
metadata = pd.DataFrame(metadata).drop_duplicates()       
metadata["cancer_type"] = metadata["project_id"].str.replace("TCGA-","")
## sort by file size
metadata = metadata.sort_values("file_size", ascending=False)
metadata = metadata.loc[metadata["cancer_type"].isin(CANCER_TYPES)]

# prepare
## get types of cancer
CANCER_TYPES = metadata["cancer_type"].to_list()
## UUIDs to download
SAMPLES = metadata["file_id"].to_list()
## get sequencer
SEQUENCER = metadata.set_index("file_id")["sequencer"].to_dict()
## tar sizes
TAR_SIZES = metadata.set_index("file_id")["file_size"].to_dict()
## tar filenames
FILENAMES = metadata.set_index("file_id")["file_name"].to_dict()
## number of samples per cancer type
N_SAMPLES = metadata.groupby("cancer_type").size().reset_index(name="n").set_index("cancer_type")["n"].to_dict()

# fastq ends
ENDS = ["1","2"]

SIZE_THRESH = 5e9

rule all:
    input:
        # create metadata
        os.path.join(TCGA_DIR,"metadata"),
        
        # Download .fastq
        expand(os.path.join(TCGA_DIR,"manifests","{cancer_type}.txt"), cancer_type=set(CANCER_TYPES)),
        expand(os.path.join(TCGA_DIR,"{cancer_type}","fastqs",".done_download"), cancer_type=set(CANCER_TYPES)),
        expand(os.path.join(TCGA_DIR,"{cancer_type}","fastqs",".done","{sample}"), zip, sample=SAMPLES, cancer_type=CANCER_TYPES),

        # Quantify splicing with vast-tools
        expand(os.path.join(TCGA_DIR,"{cancer_type}","vast_out",".done","{sample}"), zip, sample=SAMPLES, cancer_type=CANCER_TYPES),

        # Combine into single tables
        expand(os.path.join(TCGA_DIR,"{cancer_type}","vast_out",".done_combine-{n_samples}"), zip, cancer_type=N_SAMPLES.keys(), n_samples=N_SAMPLES.values()),

        # Tidy PSI
        expand(os.path.join(".done","{cancer_type}.done"), cancer_type=set(CANCER_TYPES)),
        expand(os.path.join(TCGA_DIR,"{cancer_type}","vast_out","PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz"), cancer_type=set(CANCER_TYPES)),
        
        # Download counts from GDC
        ## Download
        expand('.done/download_GDC_tcga_gene_counts_{cancer}.done', cancer=CANCER_TYPES),
        ## Prepare
        expand(os.path.join(GDC_DIR,'gene_counts','{cancer}.tsv.gz'), cancer=CANCER_TYPES)
        
        
rule create_manifests:
    output:
        os.path.join(TCGA_DIR,"manifests","{cancer_type}.txt")
    params:
        cancer_type = "{cancer_type}",
        metadata = metadata
    threads: 1
    resources:
        runtime = int(0.25*3600), # 15 min
        memory = 1
    run:
        import pandas as pd
        
        manifest = pd.DataFrame({
            "id": metadata.loc[metadata["cancer_type"]==params.cancer_type, "file_id"].values
        })
        
        manifest.to_csv(output[0], sep="\t", index=False)
            

rule create_metadata_tables:
    input:
        json1 = os.path.join(SUPPORT_DIR,"TCGA-metadata.cart.2022-07-22-brca.json"),
        json2 = os.path.join(SUPPORT_DIR,"TCGA-metadata.cart.2022-07-22-rest.json"),
        json3 = os.path.join(SUPPORT_DIR,"TCGA-metadata.cart.2022-12-27-rest_illumina_ga.json")
    output:
        folder = directory(os.path.join(TCGA_DIR,"metadata"))
    run:
        import os
        import pandas as pd
        
        # load
        metadata_json = pd.concat([
            pd.read_json(input.json1).assign(sequencer="IlluminaHiSeq"), 
            pd.read_json(input.json2).assign(sequencer="IlluminaHiSeq"), 
            pd.read_json(input.json3).assign(sequencer="IlluminaGA")
        ])
        
        # generate metadata table
        metadata = []
        for idx, row in metadata_json.iterrows():
                metadata.append({
                    "file_name": row["file_name"],
                    "file_id": row["file_id"],
                    "file_size": row["file_size"],
                    "project_id": row["cases"][0]["project"]["project_id"],
                    "sample_type": row["cases"][0]["samples"][0]["sample_type"],
                    "submitter_id": row["cases"][0]["samples"][0]["submitter_id"],
                    "is_ffpe": row["cases"][0]["samples"][0]["is_ffpe"],
                    "sequencer": row["sequencer"]
                })
        metadata = pd.DataFrame(metadata).drop_duplicates()
        metadata["file_id_vasttools"] = metadata["file_id"]+"_1"
        metadata["sampleID"] = metadata["submitter_id"].str.slice(0,15)
        metadata["cancer_type"] = metadata["project_id"].str.replace("TCGA-","")
        
        # save
        os.makedirs(output.folder, exist_ok=True)
        filename = os.path.join(output.folder,"PANCAN.tsv.gz")
        metadata.to_csv(filename, sep="\t", compression="gzip", index=False)
        
        for cancer_oi in set(metadata["cancer_type"]):
            print(cancer_oi)
            
            filename = os.path.join(output.folder,cancer_oi+".tsv.gz")
            metadata.loc[
                metadata["cancer_type"]==cancer_oi
            ].to_csv(filename, sep="\t", compression="gzip", index=False)
        
        print("Done!")
        
            
rule download_fastq:
    input:
        manifest = os.path.join(TCGA_DIR,"manifests","{cancer_type}.txt")
    output:
        download_done = os.path.join(TCGA_DIR,"{cancer_type}","fastqs",".done_download")
    params:
        token = TOKEN_FILE,
        cancer_type = "{cancer_type}",
        fastqs_dir = os.path.join(TCGA_DIR,"{cancer_type}","fastqs"),
        bin_dir = os.path.join("~/repositories")
    threads: 16
    resources:
        runtime = 3600*120,
        memory = 2
    shell:
        """
        set -eo pipefail
        
        date "+%A %W %Y %X"
        
        # we download a tar file
        echo "Downloading..."
        {params.bin_dir}/gdc-client download \
                    --token-file={params.token} \
                    --n-processes={threads} \
                    --dir={params.fastqs_dir} \
                    --no-verify \
                    --manifest={input.manifest}
        
        touch {output.download_done}
        echo "Finished downloading {params.cancer_type}."                    

        date "+%A %W %Y %X"
        
        echo "Done!"
        """

        
rule prepare_fastq:
    input:
        download_done = os.path.join(TCGA_DIR,"{cancer_type}","fastqs",".done_download")
    output:
        prepare_done = os.path.join(TCGA_DIR,"{cancer_type}","fastqs",".done","{sample}")
    params:
        cancer_type = "{cancer_type}",
        sample = "{sample}",
        tar_dir = os.path.join(TCGA_DIR,"{cancer_type}","fastqs","{sample}"),
        filename = lambda wildcards: FILENAMES[wildcards.sample],
        sequencer = lambda wildcards: SEQUENCER[wildcards.sample],
        fastqs_dir = os.path.join(TCGA_DIR,"{cancer_type}","fastqs")
    threads: 8
    resources:
        runtime = 2*3600, # 2h
        memory = 2
    shell:
        """
        set -eo pipefail

        date "+%A %W %Y %X"

        # we unpack the tar file
        echo "Unpacking..."
        tar -xvf \
            {params.tar_dir}/{params.filename} \
            --directory {params.tar_dir}/

        # we compress the fastq files
        echo "Compressing..."
        nice pigz \
            --force \
            --processes {threads} \
            {params.tar_dir}/*.fastq

        # normalize names of fastq files
        if [[ "{params.sequencer}" == "IlluminaHiSeq" ]]; then
            mv {params.tar_dir}/*_1.fastq.gz {params.fastqs_dir}/{params.sample}_1.fastq.gz
            mv {params.tar_dir}/*_2.fastq.gz {params.fastqs_dir}/{params.sample}_2.fastq.gz

        elif [[ "{params.sequencer}" == "IlluminaGA" ]]; then
            mv {params.tar_dir}/*.fastq.gz {params.fastqs_dir}/{params.sample}_1.fastq.gz
        fi

        # remove unnecessary stuff within the tar file
        echo "Cleaning up..."
        rm -r {params.tar_dir}/

        touch {output.prepare_done}
        echo "Finished preparing {params.sample}."                    

        date "+%A %W %Y %X"

        echo "Done!"
        """        
    
        
# NOTE: different stepSize=24 and trimLen=48, and bypass strand guessing with --nc
rule align:
    input:
        dbDir = os.path.join(VASTDB_DIR,"assemblies"),
        download_done = os.path.join(TCGA_DIR,"{cancer_type}","fastqs",".done","{sample}")
    output:
        align_done = touch(os.path.join(TCGA_DIR,"{cancer_type}","vast_out",".done","{sample}"))
    params:
        sample = "{sample}",
        fastqs_dir = os.path.join(TCGA_DIR,"{cancer_type}","fastqs"),
        bin_dir="~/repositories/vast-tools",
        sequencer = lambda wildcards: SEQUENCER[wildcards.sample],
        vast_out = directory(os.path.join(TCGA_DIR,"{cancer_type}","vast_out","{sample}"))
    threads: 16
    resources:
        runtime = lambda wildcards: 86400 if (TAR_SIZES[wildcards.sample]/2)>SIZE_THRESH else 21600, # most 6h is enough; some needed 24h (more reads).
        memory = 15
    group: "TCGA"
    shell:
        """
        set -eo pipefail

        date "+%A %W %Y %X"

        # align paired reads
        echo "Aligning {params.sample} {params.sequencer}..."
        if [[ "{params.sequencer}" == "IlluminaHiSeq" ]]; then
            nice {params.bin_dir}/vast-tools align \
                        {params.fastqs_dir}/{params.sample}_1.fastq.gz \
                        {params.fastqs_dir}/{params.sample}_2.fastq.gz \
                        --ns \
                        --stepSize 24 \
                        --trimLen 48 \
                        --sp Hs2 \
                        --dbDir {input.dbDir} \
                        --EEJ_counts \
                        --cores {threads} \
                        --output {params.vast_out}

        elif [[ "{params.sequencer}" == "IlluminaGA" ]]; then
            nice {params.bin_dir}/vast-tools align \
                        {params.fastqs_dir}/{params.sample}_1.fastq.gz \
                        --ns \
                        --stepSize 24 \
                        --trimLen 48 \
                        --sp Hs2 \
                        --dbDir {input.dbDir} \
                        --EEJ_counts \
                        --cores {threads} \
                        --output {params.vast_out}
        fi

        echo "Finished aligning {params.sample}."

        date "+%A %W %Y %X"

        echo "Done!"
        """

        
import datetime
t = datetime.datetime.now()
t = t.strftime('%Y%m%d%H%M%S')

rule vasttools_combine:
    input:
        done = lambda wildcards: [os.path.join(TCGA_DIR,"{cancer_type}","vast_out",".done","{sample}").format(cancer_type=wildcards.cancer_type, sample=sample) for sample in  metadata.loc[metadata["cancer_type"]==wildcards.cancer_type,"file_id"]],
        dbDir = os.path.join(VASTDB_DIR,"assemblies")
    output:
        touch(os.path.join(TCGA_DIR,"{cancer_type}","vast_out",".done_combine-{n_samples}")),
        psi = os.path.join(TCGA_DIR,"{cancer_type}","vast_out","INCLUSION_LEVELS_FULL-hg38-{n_samples}.tab.gz")
    params:
        bin_dir="~/repositories/vast-tools/",
        vast_out = os.path.join(TCGA_DIR,"{cancer_type}","vast_out"),
        folder = os.path.join(TCGA_DIR,"{cancer_type}","vast_out",t)
    threads: 16
    resources:
        runtime = 3600*24, # 24h
        memory = 90
    shell:
        """
        set -eo pipefail

        # group results
        echo "Grouping results..."
        mkdir -p {params.folder}/to_combine
        ln -s {params.vast_out}/*/to_combine/* {params.folder}/to_combine/
        mkdir -p {params.folder}/expr_out
        ln -s {params.vast_out}/*/expr_out/* {params.folder}/expr_out/

        # combine runs
        echo "Combining runs..."
        nice {params.bin_dir}/vast-tools combine \
                    --cores {threads} \
                    --sp Hs2 \
                    --dbDir {input.dbDir} \
                    --keep_raw_reads \
                    --keep_raw_incl \
                    --output {params.folder}

        # compress outputs
        echo "Compressing outputs..."
        gzip -f {params.folder}/raw_incl/*
        gzip -f {params.folder}/raw_reads/*
        gzip -f {params.folder}/*.tab

        # move results from tmp to vast_out
        mv {params.folder}/raw_incl {params.vast_out}/
        mv {params.folder}/raw_reads {params.vast_out}/
        mv {params.folder}/*.tab.gz {params.vast_out}/

        # remove tmp directory
        echo "Removing grouped results..."
        rm -r {params.folder}


        echo "Done!"
        """
    
    
rule vasttools_tidy:
    input:
        lambda wildcards: os.path.join(TCGA_DIR,"{cancer_type}","vast_out","INCLUSION_LEVELS_FULL-hg38-{n_samples}.tab.gz").format(cancer_type=wildcards.cancer_type, n_samples=N_SAMPLES[wildcards.cancer_type]) ## combined table
    output:
        touch(os.path.join(".done","{cancer_type}.done")),
        tidy = os.path.join(TCGA_DIR,"{cancer_type}","vast_out","PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz")
    params:
        bin_dir="~/repositories/vast-tools/"
    threads: 1
    resources:
        runtime = 3600*12, # 12h
        memory = 40
    shell:
        """
        set -eo pipefail

        echo "Tidying up..."
        {params.bin_dir}/vast-tools tidy <(zcat {input}) -min_N 1 -min_SD 0 --min_ALT_use 25 --noVLOW --log -outFile {output.tidy}
        gzip --force {output.tidy}
        mv {output.tidy}.gz {output.tidy}
        
        echo "Done!"
        """
        
        
rule download_tcga_genexpr_counts:
    output:
        touch('.done/download_GDC_tcga_gene_counts_{cancer}.done')
    params:
        download_dir = GDC_DIR
    shell:
        """
        nice Rscript scripts/tcga__download_from_GDC.R \
                    --cancer_type={wildcards.cancer} \
                    --download_dir={params.download_dir} \
                    --data_category='Transcriptome Profiling' \
                    --data_type='Gene Expression Quantification' \
                    --workflow_type='STAR - Counts'
        """
        
        
rule make_gene_counts_matrices_from_GDC:
    input:
        '.done/download_GDC_tcga_gene_counts_{cancer}.done'
    output:
        expr_mat = os.path.join(GDC_DIR,'gene_counts','{cancer}.tsv.gz')
    params:
        download_dir = GDC_DIR
    shell:
        """
        nice Rscript scripts/tcga__make_omic_matrix_from_GDC.R \
                    --cancer_type={wildcards.cancer} \
                    --output_file={output.expr_mat} \
                    --download_dir={params.download_dir} \
                    --data_category='Transcriptome Profiling' \
                    --data_type='Gene Expression Quantification' \
                    --workflow_type='STAR - Counts'
        """