---
name: read-qc-alignment
description: Use when starting from raw FASTQ reads and you need resequencing QC, trimming, reference alignment, sorted BAM generation, and mapping QC for breeding projects.
---

# read-qc-alignment

## When to use

- The input is still raw sequencing data, usually paired-end `FASTQ.gz`
- The goal is to produce high-quality BAM files for downstream variant calling
- You want to preserve both pre-trimming and post-trimming QC outputs

## Quick start

- Prepare the sample sheet and reference genome described in `references/input_contract.md`
- Copy and edit `assets/read_qc_alignment.env.template`
- Run `bash scripts/run_read_qc_alignment.sh --config assets/read_qc_alignment.env.template`

## Default workflow

1. Validate the sample sheet, FASTQ paths, and reference genome indexes
2. Run raw-read QC and collect FastQC-style metrics
3. Trim adapters and low-quality bases with `fastp` or an equivalent tool
4. Align cleaned reads with `bwa mem` or another short-read aligner
5. Sort, index, and summarize BAM alignment metrics
6. Export a BAM manifest ready for downstream variant calling

## Primary outputs

- `sample_manifest.tsv`
- `workflow_plan.txt`
- `qc_raw/` and `qc_trimmed/`
- `bam/` and `stats/`

## Read only when needed

- Input columns and directory requirements: `references/input_contract.md`
- Tool choices and key checkpoints: `references/workflow.md`
