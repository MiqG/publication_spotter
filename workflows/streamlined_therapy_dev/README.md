# Application of splicing dependency models

Once fitted and selected all splicing dependency models by running the `model_splicing_dependency` workflow, this workflow will use these models to prioritize clinically-relevant cancer-driver exons with therapeutic potential.

## Important remarks

Make sure to have all dependencies installed.

## Recommendations
Run the workflow in a machine with at least 20GB of RAM using
```
snakemake --cores 6
```
In case you want to run the rules on your cluster, refer to [snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cluster.html).
