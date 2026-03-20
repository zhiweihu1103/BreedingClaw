---
name: variant-call-filter
description: Use when you already have aligned BAM files and need cohort-level variant calling, normalization, and filtering to obtain a clean VCF for downstream population genetics or GWAS.
---

# variant-call-filter

## When to use

- Upstream read QC and reference alignment are already finished
- The goal is a clean, traceable cohort VCF
- You want to retain raw calls, normalized VCFs, and filtered VCFs

## Quick start

- Prepare `bam_list.tsv` as described in `references/input_contract.md`
- Copy and edit `assets/variant_call_filter.env.template`
- Run `bash scripts/run_variant_call_filter.sh --config assets/variant_call_filter.env.template`

## Default workflow

1. Validate the BAM manifest and reference genome
2. Choose a caller such as `bcftools` or `GATK` based on project policy
3. Normalize the raw VCF, split multiallelic sites, and collect basic statistics
4. Apply quality, depth, missingness, and allele-frequency filters
5. Export filter reports, site statistics, and the final cohort VCF

## Primary outputs

- `bam_manifest.tsv`
- `bam_paths.list`
- `workflow_plan.txt`
- `raw_calls/`
- `normalized/`
- `filtered/`

## Read only when needed

- BAM manifest rules: `references/input_contract.md`
- Caller selection and filter threshold guidance: `references/workflow.md`
