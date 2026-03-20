#!/usr/bin/env bash
set -euo pipefail

resolve_bin() {
  local current_value="$1"
  shift

  if [[ -n "$current_value" ]] && command -v "$current_value" >/dev/null 2>&1; then
    printf '%s\n' "$current_value"
    return 0
  fi

  local candidate
  for candidate in "$@"; do
    if [[ -n "$candidate" ]] && command -v "$candidate" >/dev/null 2>&1; then
      printf '%s\n' "$candidate"
      return 0
    fi
  done

  return 1
}

usage() {
  cat <<'EOF'
Usage: bash scripts/run_genotype_kinship_prep.sh --config path/to/genotype_kinship_prep.env
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

required_vars=(VCF_FILE PHENO_TSV OUTDIR PREFIX THREADS PCA_COMPONENTS)
for var_name in "${required_vars[@]}"; do
  [[ -n "${!var_name:-}" ]] || { echo "Missing variable: $var_name" >&2; exit 1; }
done

EMMAX_BIN="${EMMAX_BIN:-}"
KINSHIP_BIN="${KINSHIP_BIN:-}"
if ! EMMAX_BIN="$(resolve_bin "$EMMAX_BIN" emmax-intel64 emmax)"; then
  echo "Unable to locate the EMMAX association binary. Set EMMAX_BIN or add emmax-intel64 to PATH." >&2
  exit 1
fi
if ! KINSHIP_BIN="$(resolve_bin "$KINSHIP_BIN" emmax-kin-intel64 emmax_kin-intel64 emmax_kin emmax-kin)"; then
  echo "Unable to locate the EMMAX kinship binary. Set KINSHIP_BIN or add emmax-kin-intel64 to PATH." >&2
  exit 1
fi

mkdir -p "$OUTDIR/intermediate" "$OUTDIR/covariates" "$OUTDIR/emmax" "$OUTDIR/logs"

awk -F'\t' '
NR==1 {
  if ($1 != "sample_id") {
    print "phenotype header must start with sample_id" > "/dev/stderr"
    exit 2
  }
  for (i = 2; i <= NF; i++) {
    print $i
  }
}
' "$PHENO_TSV" > "$OUTDIR/trait_list.txt"

cat > "$OUTDIR/workflow_plan.txt" <<EOF
plink --vcf ${VCF_FILE} --make-bed --out ${OUTDIR}/intermediate/${PREFIX}
plink --vcf ${VCF_FILE} --recode12 --transpose --output-missing-genotype 0 --out ${OUTDIR}/emmax/${PREFIX}
plink --bfile ${OUTDIR}/intermediate/${PREFIX} --pca ${PCA_COMPONENTS} header tabs --out ${OUTDIR}/covariates/${PREFIX}
${KINSHIP_BIN} -v -s -d 10 ${OUTDIR}/emmax/${PREFIX}
${EMMAX_BIN} -t ${OUTDIR}/emmax/${PREFIX} -p <trait.pheno.tsv> -k <kinship.kf> -c <covariates.tsv> -o <trait_prefix>
EOF

printf 'Prepared %s and %s\n' "$OUTDIR/trait_list.txt" "$OUTDIR/workflow_plan.txt"
