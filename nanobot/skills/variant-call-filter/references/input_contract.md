# Input Contract

## Required files

- `bam_list.tsv`
- Reference genome `genome.fa`

## `bam_list.tsv` columns

- The header must be `sample_id	bam`
- `bam` may be an absolute path or a path relative to the working directory
- All BAM files must already be sorted and indexed

## Recommended minimal metadata

- Each sample appears only once
- Do not mix BAM files generated against different reference versions
- Retain batch information for downstream troubleshooting

## Final VCF recommendations

- Keep raw calls, normalized results, and filtered results together
- Record all filter thresholds in the workflow configuration
