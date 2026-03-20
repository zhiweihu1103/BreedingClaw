# Workflow Notes

## Default tools

- Raw-read QC: `fastqc`
- Trimming: `fastp`
- Alignment: `bwa mem`
- Sorting and indexing: `samtools`

## Key decisions

- Prefer `bwa mem` for short-read population resequencing unless the project specifies another aligner
- If adapter sets or lane naming differ, normalize names before execution
- If joint calling is planned later, keep read group naming consistent at this stage

## Delivery standard

- Every sample produces at least one sorted and indexed BAM
- Every sample has a `flagstat` file or an equivalent mapping summary
- Both pre-trimming and post-trimming QC outputs are retained for traceability
