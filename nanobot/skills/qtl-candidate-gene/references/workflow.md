# Workflow Notes

## Interval definition strategy

- Prefer LD-supported boundaries derived from the lead SNP and nearby genotype data
- If a precomputed LD table is available, reuse it instead of recalculating
- If LD is not available yet, use a fixed window as the fallback interval
- Nearby peaks on the same trait and chromosome should be merged or resolved to one lead signal

## Candidate gene guidance

- Keep all genes inside each QTL interval
- Reserve separate columns for representative genes, nearest genes, and functional notes
- If annotation databases are available, append functional keywords in dedicated columns

## Plotting guidance

- Use `scripts/plot_ld_qtl_summary_reference.R` as the reference renderer for layered LD/QTL summary figures
- Use `references/plot_input_contract.md` to prepare plotting tables
- Use `references/plot_style.md` as the style contract for colors, panel layout, and highlight logic

## Delivery standard

- Every QTL can be traced back to its lead SNP and significance threshold
- The candidate gene table records the applied window or LD rule
- All outputs are ready for downstream plotting or report assembly
