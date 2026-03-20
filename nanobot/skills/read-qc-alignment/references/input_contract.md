# Input Contract

## Required files

- `sample_sheet.tsv`
- Reference genome `genome.fa`
- A directory containing paired-end raw FASTQ files

## `sample_sheet.tsv` columns

- The header must be `sample_id	read1	read2`
- `read1` and `read2` may be file names relative to `FASTQ_DIR` or absolute paths
- `sample_id` values must be unique

## Minimum validation checks

- Sample names stay consistent across all intermediate files
- The number of FASTQ pairs matches the number of sample sheet rows
- The reference genome already has `bwa` and `samtools faidx` indexes

## Recommended output layout

- `outputs/read_qc_alignment/qc_raw`
- `outputs/read_qc_alignment/qc_trimmed`
- `outputs/read_qc_alignment/bam`
- `outputs/read_qc_alignment/stats`
