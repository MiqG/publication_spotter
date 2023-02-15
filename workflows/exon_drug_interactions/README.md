# Data preprocessing workflow

Once obtained the models of splicing dependency from the `model_splicing_dependency` workflow, this workflow associates their predictions with drug sensitivities.

## Important remarks

Make sure to have all packages required to run the scripts.

## Recommendations
Run the workflow in a machine with at least 20GB of RAM using
```
snakemake --cores 6
```
In case you want to run the rules on your cluster, refer to [snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cluster.html).