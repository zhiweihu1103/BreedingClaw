# Plot Input Contract

## Minimal plotting tables

### `chrom_sizes`

- `chrom`
- `length_bp`

Controls chromosome offsets in the genome-wide Manhattan overview panel.

### `global_manhattan`

- `snp_id`
- `chrom`
- `bp`
- `p_value` or `mlog10_p`

Provides the genome-wide Manhattan scatter points.

### `qtl_regions`

- `locus_id`
- `chrom`
- `lead_snp_id`
- `lead_bp`
- `lead_p` or `lead_mlog10_p`
- `qtl_start`
- `qtl_end`
- `panel_start` (optional)
- `panel_end` (optional)
- `qtl_type` (optional)

Defines the highlighted QTL intervals and local plotting windows.

### `local_manhattan`

- `locus_id`
- `snp_id`
- `chrom`
- `bp`
- `mlog10_p` or `p_value`
- `r2` (optional)
- `is_lead` (optional)

Provides local association points for each QTL zoom panel. In the upgraded LD-aware workflow, `r2` should come from real lead-SNP LD calculations whenever genotype inputs are available.

### `gene_models`

- `gene_id`
- `chrom`
- `start_bp`
- `end_bp`
- `strand`

Provides the lower gene track panel for each locus.

### `highlight_genes` (optional)

- `locus_id`
- `gene_id`
- `rank` (optional)
- `label` (optional)

Controls which genes are highlighted and optionally relabeled.

## Rendering behavior

- If `chrom_sizes` is missing, chromosome lengths are inferred from `global_manhattan` and `qtl_regions`.
- If `panel_start` and `panel_end` are missing, they are expanded from `qtl_start` and `qtl_end` using the configured flank distance.
- If `r2` is missing, local points are drawn with the neutral low-LD color.
- If genotype prefixes are provided to the builder, `r2` is expected to reflect real PLINK-based LD rather than fixed-window placeholders.
- If `gene_models` is missing, the lower panel is rendered as an empty track with a message.

## Reference renderer

- `scripts/plot_ld_qtl_summary_reference.R`

This renderer is an R translation of the layered LD/QTL summary figure style extracted from the `GWAS_codex` project outputs.
