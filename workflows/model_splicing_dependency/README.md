# Inference of splicing dependency models

Once preprocessd all files by running the `preprocess_data` workflow, this workflow will generate the models of splicing dependency.

## Important remarks

Make sure to have all dependencies installed.

## Recommendations
Run the workflow in a machine with at least 20GB of RAM using
```
snakemake --cores 6
```
In case you want to run the rules on your cluster, refer to [snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cluster.html).
