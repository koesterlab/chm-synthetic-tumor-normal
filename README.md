# Snakemake workflow: generating a synthetic tumor/normal dataset from the haploid cell lines CHM1 and CHM13

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.7.0-brightgreen.svg)](https://snakemake.readthedocs.io)

This workflow generates a synthetic tumor/normal dataset based on the [syndip benchmark](https://doi.org/10.1038/s41592-018-0054-7).

## Authors

* Johannes Köster (@johanneskoester)

## Usage

In any case, if you use this workflow in a paper, don't forget to give credits to the authors (and the authors of the syndip benchmark) by citing the URL of this (original) repository.

Perform a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --use-conda --cluster qsub --jobs 100

or

    snakemake --use-conda --drmaa --jobs 100

See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.
