# Input Contract

## Required files

- A cohort VCF or an already filtered PLINK dataset
- `phenotype.tsv`

## `phenotype.tsv` columns

- The first column must be `sample_id`
- Columns from the second onward represent one or more traits
- Use `NA` as the preferred missing-value marker

## Sample consistency requirements

- Genotype and phenotype files must use the same sample naming scheme
- Remove duplicate samples before calculating PCA and kinship
- If population-structure covariates already exist, provide them as a separate covariate table

## EMMAX output requirements

- `PREFIX.tped`
- `PREFIX.tfam`
- A kinship matrix
- Per-trait phenotype vectors that can be extracted independently
