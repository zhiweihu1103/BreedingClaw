# Input Contract

## Required files

- `PREFIX.tped`
- `PREFIX.tfam`
- A kinship matrix
- `phenotype.tsv`

## Optional files

- `covariates.tsv`

## `phenotype.tsv` columns

- The first column must be `sample_id`
- Columns from the second onward are trait names
- Trait names should use letters, numbers, and underscores when possible

## Recommended result columns

- `trait`
- `chrom`
- `pos`
- `snp_id`
- `pvalue`
- `beta`
- `stderr`
