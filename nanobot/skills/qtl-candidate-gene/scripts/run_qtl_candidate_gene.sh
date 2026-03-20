#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
  cat <<'EOF'
Usage: bash scripts/run_qtl_candidate_gene.sh --config path/to/qtl_candidate_gene.env
EOF
}

config=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      config="${2:-}"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

[[ -n "$config" ]] || { usage; exit 1; }
[[ -f "$config" ]] || { echo "Config not found: $config" >&2; exit 1; }

set -a
source "$config"
set +a

required_vars=(GWAS_TABLE GENE_ANNOTATION OUTDIR P_THRESHOLD WINDOW_BP)
for var_name in "${required_vars[@]}"; do
  [[ -n "${!var_name:-}" ]] || { echo "Missing variable: $var_name" >&2; exit 1; }
done

TRAIT_LABEL="${TRAIT_LABEL:-$(basename "$GWAS_TABLE" | sed 's/\.gwas_summary\.tsv$//' | sed 's/\.all\.QQ\.pmap\.txt$//')}"
REFERENCE_FAI="${REFERENCE_FAI:-}"
MAX_LOCI="${MAX_LOCI:-3}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
R_BIN="${R_BIN:-Rscript}"
PLOT_SCRIPT="${PLOT_SCRIPT:-${script_dir}/plot_ld_qtl_summary_reference.R}"
QTL_BUILDER="${QTL_BUILDER:-${script_dir}/build_qtl_fixed_window_bundle.py}"
LD_TABLE="${LD_TABLE:-}"
GENO_TFILE_PREFIX="${GENO_TFILE_PREFIX:-}"
GENO_BFILE_PREFIX="${GENO_BFILE_PREFIX:-}"
BED_CACHE_PREFIX="${BED_CACHE_PREFIX:-}"
PLINK_BIN="${PLINK_BIN:-plink}"
PLINK_THREADS="${PLINK_THREADS:-1}"
LD_R2="${LD_R2:-0.2}"
LOCAL_FLANK_BP="${LOCAL_FLANK_BP:-200000}"
MIN_QTL_PAD_BP="${MIN_QTL_PAD_BP:-50000}"

mkdir -p "$OUTDIR/qtl_regions" "$OUTDIR/candidates" "$OUTDIR/plots" "$OUTDIR/logs"

builder_cmd=(
  "$PYTHON_BIN" "$QTL_BUILDER"
  --gwas-table "$GWAS_TABLE"
  --gene-annotation "$GENE_ANNOTATION"
  --outdir "$OUTDIR/qtl_regions"
  --p-threshold "$P_THRESHOLD"
  --window-bp "$WINDOW_BP"
  --trait-label "$TRAIT_LABEL"
  --reference-fai "$REFERENCE_FAI"
  --max-loci "$MAX_LOCI"
  --ld-r2 "$LD_R2"
  --local-flank-bp "$LOCAL_FLANK_BP"
  --min-qtl-pad-bp "$MIN_QTL_PAD_BP"
  --plink-bin "$PLINK_BIN"
  --plink-threads "$PLINK_THREADS"
)

if [[ -n "$LD_TABLE" && -f "$LD_TABLE" ]]; then
  builder_cmd+=(--ld-table "$LD_TABLE")
fi
if [[ -n "$GENO_TFILE_PREFIX" ]]; then
  builder_cmd+=(--geno-tfile-prefix "$GENO_TFILE_PREFIX")
fi
if [[ -n "$GENO_BFILE_PREFIX" ]]; then
  builder_cmd+=(--geno-bfile-prefix "$GENO_BFILE_PREFIX")
fi
if [[ -n "$BED_CACHE_PREFIX" ]]; then
  builder_cmd+=(--bed-cache-prefix "$BED_CACHE_PREFIX")
fi

"${builder_cmd[@]}"

cp "$OUTDIR/qtl_regions/${TRAIT_LABEL}.top10_candidate_genes.tsv" "$OUTDIR/candidates/" 2>/dev/null || true

"$R_BIN" "$PLOT_SCRIPT" \
  --trait-label "$TRAIT_LABEL" \
  --chrom-sizes "$OUTDIR/qtl_regions/plot_tables/chrom_sizes.tsv" \
  --global-manhattan "$OUTDIR/qtl_regions/plot_tables/global_manhattan.tsv" \
  --qtl-regions "$OUTDIR/qtl_regions/plot_tables/qtl_regions.tsv" \
  --local-manhattan "$OUTDIR/qtl_regions/plot_tables/local_manhattan.tsv" \
  --gene-models "$OUTDIR/qtl_regions/plot_tables/gene_models.tsv" \
  --highlight-genes "$OUTDIR/qtl_regions/plot_tables/highlight_genes.tsv" \
  --out "$OUTDIR/plots/${TRAIT_LABEL}.ld_qtl_summary.pdf" \
  --max-loci "$MAX_LOCI"

cat > "$OUTDIR/workflow_plan.txt" <<EOF
1. Build significant loci from ${GWAS_TABLE} using p <= ${P_THRESHOLD}.
2. Refine each locus with real LD from ${LD_TABLE:-${GENO_BFILE_PREFIX:-${GENO_TFILE_PREFIX:-<not_provided>}}} when available, otherwise fall back to fixed windows.
3. Intersect the refined QTL intervals with ${GENE_ANNOTATION} to rank candidate genes.
4. Export QTL summary tables and plotting tables to ${OUTDIR}/qtl_regions.
5. Export candidate gene tables to ${OUTDIR}/candidates.
6. Render the LD/QTL summary figure to ${OUTDIR}/plots/${TRAIT_LABEL}.ld_qtl_summary.pdf using ${PLOT_SCRIPT}.
EOF

printf 'Prepared %s, %s, and %s\n' "$OUTDIR/qtl_regions/${TRAIT_LABEL}.qtl_summary.tsv" "$OUTDIR/candidates/${TRAIT_LABEL}.top10_candidate_genes.tsv" "$OUTDIR/plots/${TRAIT_LABEL}.ld_qtl_summary.pdf"
