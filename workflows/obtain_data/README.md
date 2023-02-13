# Data download workflow

This workflow will download the files required to complete all the analyses.

## Important remarks
Note that two datasets require you to have private tokens:

- `vast-tools` human genome assembly
Make sure to place it at the following path:
```
publication_spotter/data/raw/VastDB/assemblies/Hs2
```

- COSMIC
download manually from https://cancer.sanger.ac.uk/census. Expected output directory at:
```
publication_spotter/data/raw/COSMIC/cancer_gene_census.tsv
```

- AS Cancer Atlas
Download manually from https://ngdc.cncb.ac.cn/ascancer/download. Expected output directory at:
```
publication_spotter/data/raw/ASCancerAtlas/CASE_all.csv.gz
```

- TCGA
Download from [GDC](https://portal.gdc.cancer.gov/) with dbGaP access. Private token expected at the following path:
```
publication_spotter/support/.private/gdc-user-token.txt
```

- Pietilla 2021
Download from [EGA](https://ega-archive.org/studies/EGAS00001004714). The list of fastq files, token and decryption key and pass files are expected to be found at:
```
# list of fastq files
publication_spotter/support/.private/ega_files-EGAD00001006456.txt

# token
publication_spotter/support/.private/ega_aspera_pass.txt

# decryption key
publication_spotter/support/.private/ega_decryption_key

# decryption pass
publication_spotter/support/.private/ega_decryption_pass
```

## Recomendations
Due to their computational burden, we recommend running all `download_and_align*` workflows separately in a cluster before running the main `snakefile` workflow. See [snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cluster.html) for further information on how to run it on your specific computer cluster configuration.
