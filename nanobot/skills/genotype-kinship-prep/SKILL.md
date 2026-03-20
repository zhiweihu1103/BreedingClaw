---
name: genotype-kinship-prep
description: Use when you need to transform filtered genotype data plus phenotype and covariate tables into EMMAX-ready inputs including TPED/TFAM, kinship, PCA covariates, and sample-intersected trait matrices.
---

# genotype-kinship-prep

## When to use

- You already have a filtered VCF or PLINK dataset
- You need to reconcile sample intersections, covariates, PCA, and kinship
- The downstream association engine is fixed to `emmax-intel64`

## Quick start

- Prepare genotype and phenotype matrices as described in `references/input_contract.md`
- Copy and edit `assets/genotype_kinship_prep.env.template`
- Use `emmax-intel64` and `emmax-kin-intel64` as the default binary names; do not assume a plain `emmax` executable exists in `PATH`
- Run `bash scripts/run_genotype_kinship_prep.sh --config assets/genotype_kinship_prep.env.template`

## Default workflow

1. Validate genotype, phenotype, and sample identifiers
2. Convert the VCF into EMMAX-compatible `tped/tfam`
3. Build PLINK binaries and compute PCA covariates
4. Generate the kinship matrix
5. Extract trait names and write a trait manifest for batch GWAS

## Primary outputs

- `trait_list.txt`
- `workflow_plan.txt`
- `emmax/`
- `covariates/`
- `intermediate/`

## Read only when needed

- Input matrix rules: `references/input_contract.md`
- `emmax-intel64` preparation details: `references/workflow.md`
