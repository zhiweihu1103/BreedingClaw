---
name: qtl-candidate-gene
description: Use when GWAS summary statistics are available and you need to collapse significant SNPs into QTL regions, refine intervals with LD or fixed windows, annotate nearby genes, and export candidate gene tables for breeding interpretation.
---

# qtl-candidate-gene

## When to use

- You already have significant GWAS hits or full summary statistics
- You need to move from SNP-level signals to QTL intervals
- You want to collect candidate genes using annotation, LD, or fixed windows

## Quick start

- Prepare GWAS results and gene annotation as described in `references/input_contract.md`
- Copy and edit `assets/qtl_candidate_gene.env.template`
- Run `bash scripts/run_qtl_candidate_gene.sh --config assets/qtl_candidate_gene.env.template`

## Default workflow

1. Filter significant hits using a project-level p-value threshold
2. Group lead SNPs by trait and chromosome
3. Build LD-supported QTL intervals from genotype data and fall back to fixed windows only when LD is unavailable
4. Intersect refined QTL intervals with gene annotation to collect candidate genes
5. Export QTL summaries, candidate gene tables, and plotting-ready QTL context tables

## Primary outputs

- `significant_hits.tsv`
- `workflow_plan.txt`
- `qtl_regions/`
- `candidates/`
- `plots/`

## Read only when needed

- Input summary-table requirements: `references/input_contract.md`
- QTL and candidate-gene strategy notes: `references/workflow.md`
- Plotting table requirements: `references/plot_input_contract.md`
- Plot colors, panel layout, and highlight style: `references/plot_style.md`
