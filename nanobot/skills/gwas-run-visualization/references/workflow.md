# Workflow Notes

## Batch execution guidance

- Generate a unified trait list before writing per-trait phenotype vectors
- Give every trait its own output prefix to avoid overwriting
- Separate logs from association outputs so failed traits are easier to trace

## Minimum visualization set

- Manhattan plot
- QQ plot
- Lead SNP table
- Per-trait significant-hit summary
- Use `scripts/plot_manhattan_qq_reference.R` as the reference implementation
- Use `references/plot_style.md` as the style contract for colors, thresholds, and layout

## Delivery standard

- Every trait maps to one unique result prefix
- Result tables contain chromosome, position, p-value, and SNP identifiers
- Visualization input columns remain consistent across all traits
