# Data preprocessing workflow

Once obtained the models of splicing dependency from the `exon_drug_interactions` workflow, this workflow uses the transcriptomes from [Pietila2021](https://www.nature.com/articles/s41467-021-24009-8) from untreated patient tumors to predict their sensitivity to drug treatment.

## Important remarks

Make sure to have all packages required to run the scripts.

## Recommendations
Run the workflow in a machine with at least 20GB of RAM using
```
snakemake --cores 6
```
In case you want to run the rules on your cluster, refer to [snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cluster.html).