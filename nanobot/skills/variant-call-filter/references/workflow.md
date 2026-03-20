# Workflow Notes

## Caller selection

- Use `bcftools` first when the cohort size is moderate and fast iteration matters
- Switch to `GATK HaplotypeCaller` plus joint genotyping when the project already follows GATK Best Practices

## Filtering guidance

- Normalize sites before applying quality filters
- Include at least `QUAL`, depth, and missingness in the baseline filters
- Population datasets often also need MAF or MAC constraints

## Delivery standard

- The filtered VCF can be read cleanly by `plink` or `bcftools stats`
- Site-count changes remain explainable across each stage
- Sample order remains traceable back to the original BAM manifest
