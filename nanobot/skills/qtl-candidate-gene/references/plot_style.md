# Plot Style Reference

This style contract is extracted from the LD/QTL summary figure used in the `GWAS_codex` workflow that produces files such as `AcCa_13_0__H.ld_qtl_summary.pdf`.

## Overall Figure Structure

- One genome-wide Manhattan overview panel on top
- Two stacked local panels per selected QTL:
  - local Manhattan panel
  - gene track panel
- Figure size follows `14 x (4.6 + 2.7 * n_loci)` inches
- Panel height ratios follow `2.9` for the top panel and `2.1 : 0.8` for each local Manhattan : gene pair

## Genome-wide Manhattan Style

- Scatter points use alternating chromosome colors:
  - odd chromosomes: `#4C78A8`
  - even chromosomes: `#F58518`
- Point size is visually compact (`cex ≈ 0.35`)
- Point transparency is about `0.75`
- Background points are sampled for plotting speed, while QTL intervals and stronger points are retained
- Each QTL interval is highlighted with a translucent rectangle at alpha `0.12`
- QTL colors cycle in this order:
  - `#D73027`
  - `#4575B4`
  - `#1A9850`
  - `#984EA3`
  - `#FF7F00`
- Lead SNPs are shown as filled circles with black borders
- Locus labels are placed just above the lead point
- The genome-wide significance line is a dark red dashed line:
  - color: `#b2182b`
  - line width: `0.8`

## Local Manhattan Style

- The QTL interval is shaded with the same locus color at alpha `0.12`
- Local points are colored by LD (`r2`) bins computed from the lead SNP whenever genotype data are available:
  - `r2 >= 0.8`: `#b2182b`
  - `0.6 <= r2 < 0.8`: `#ef8a62`
  - `0.4 <= r2 < 0.6`: `#fddbc7`
  - `0.2 <= r2 < 0.4`: `#67a9cf`
  - `r2 < 0.2` or missing: `#d1d1d1`
- Local point size is medium (`cex ≈ 0.75`)
- Point transparency is about `0.9`
- Background points are sampled for plotting speed, while the lead SNP, higher-LD points, and stronger local signals are retained
- The lead SNP is emphasized as a black diamond with a white border
- Panel titles are left-aligned and include:
  - locus ID
  - chromosome and QTL interval
  - lead SNP ID
  - QTL type

## Gene Track Style

- Genes are drawn as horizontal segments
- Overlapping genes are stacked into lanes with `15000 bp` padding
- Non-highlight genes use:
  - color `#555555`
  - line width `1.4`
- Highlight genes use the locus color and line width `2.5`
- Highlight labels are placed above the start coordinate on `+` strand genes and above the end coordinate on `-` strand genes
- If no genes are available, the panel shows a neutral message in gray

## Output Character

- The final figure is a layered, publication-style summary panel
- It emphasizes three linked scales simultaneously:
  - whole-genome significance context
  - local LD-supported signal structure
  - candidate gene context
