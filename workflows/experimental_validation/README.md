# Experimental validation of selected cancer-driver exons with therapeutic potential

Once selected clinically-relevant exon targets with therapeutic poential by running the `stremlined_therapy_dev` workflow, this workflow will prioritize the vulnerabilities of our laboratory cancer cells to validate experimentally the interesting targets.

## Important remarks

Make sure to have all dependencies installed.

## Recommendations
Run the workflow in a machine with at least 20GB of RAM using
```
snakemake --cores 6
```
In case you want to run the rules on your cluster, refer to [snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cluster.html).
