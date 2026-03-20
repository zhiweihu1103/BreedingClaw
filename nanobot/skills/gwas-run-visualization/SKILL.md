---
name: gwas-run-visualization
description: Use when EMMAX-ready genotype, kinship, phenotype, and optional covariate files are already prepared and you need batch GWAS execution, result collection, and visualization inputs such as Manhattan, QQ, and lead SNP tables.
---

# gwas-run-visualization

## When to use

- `tped/tfam`, kinship, and phenotype inputs are already prepared
- You need to run `emmax-intel64` across many traits
- You want a consistent set of Manhattan, QQ, and lead SNP outputs

## Quick start

- Prepare EMMAX inputs and phenotype matrices described in `references/input_contract.md`
- Copy and edit `assets/gwas_run_visualization.env.template`
- Use `emmax-intel64` as the default association binary; do not probe plain `emmax` first unless your environment explicitly uses that name
- Run `bash scripts/run_gwas_visualization.sh --config assets/gwas_run_visualization.env.template`

## Default workflow

1. Validate `tped/tfam`, kinship, phenotype, and optional covariate inputs
2. Build a trait list from the phenotype matrix
3. Generate per-trait phenotype vectors and EMMAX commands
4. Collect association results, significant hits, and plotting tables
5. Export a standard directory structure for Manhattan, QQ, and summary figures

## Primary outputs

- `trait_list.txt`
- `traits/`
- `assoc/`
- `run_emmax_commands.sh`
- `plot_plan.txt`

## Read only when needed

- Input file and column requirements: `references/input_contract.md`
- Result directory organization guidance: `references/workflow.md`
- Plot colors, panel layout, and threshold style: `references/plot_style.md`
