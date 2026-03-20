# Plot Style Reference

This reference captures the visual style extracted from the original `Manhattan_QQ_plot.R` and keeps it as the default plotting template for `gwas-run-visualization`.

## Overall Style

- White background using `theme_bw()`
- No panel grid lines
- Bold, centered plot titles
- Clean, publication-style layout with minimal decoration

## Manhattan Plot Style

- One point layer colored by chromosome
- Point size: `1.5`
- Chromosome labels shown as `Chr1`, `Chr2`, and so on
- The x-axis is built from cumulative chromosome offsets
- The legend is hidden
- Background points are rendered with sampled-background thinning, while stronger signal points are retained
- Two dashed threshold lines are shown:
  - suggestive line at `-log10(1 / n_snps)`
  - Bonferroni line at `-log10(0.05 / n_snps)` in red

## Manhattan Color Palette

The chromosome palette follows the original qualitative sequence:

- `#e41a1c`
- `#377eb8`
- `#4daf4a`
- `#984ea3`
- `#ff7f00`
- `#addd33`
- `#a65628`
- `#f781bf`
- `#999999`

If the chromosome count exceeds nine, the palette is repeated in order.

## QQ Plot Style

- Point color: `steelblue`
- Point size: `1.2`
- The strongest tail is retained while background points are sampled for speed
- A red dashed `y = x` reference line
- White background with grid removed
- Square standalone panel in the default export set

## Layout and Export Style

- Combined figure layout uses `cowplot::plot_grid()`
- Relative panel widths: Manhattan `4`, QQ `1.5`
- Combined figure size: `12 x 6` inches
- Standalone Manhattan size: `10 x 6` inches
- Standalone QQ size: `5 x 5` inches
- Default exports include both `PDF` and `PNG`

## Data Handling Rules Preserved from the Original Script

- Replace `P = 0` with the minimum non-zero p-value from the QQ input
- Use the QQ input row count to define both suggestive and Bonferroni thresholds
- For very large QQ inputs, retain all strongest `1%` of points and downsample the remainder to a maximum of `200000` points
- Export plotting tables alongside the figures for traceability
