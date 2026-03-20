#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

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
Usage: bash scripts/run_gwas_visualization.sh --config path/to/gwas_run_visualization.env
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

required_vars=(TPED_PREFIX KINSHIP_FILE PHENO_TSV OUTDIR)
for var_name in "${required_vars[@]}"; do
  [[ -n "${!var_name:-}" ]] || { echo "Missing variable: $var_name" >&2; exit 1; }
done

EMMAX_BIN="${EMMAX_BIN:-}"
if ! EMMAX_BIN="$(resolve_bin "$EMMAX_BIN" emmax-intel64 emmax)"; then
  echo "Unable to locate the EMMAX association binary. Set EMMAX_BIN or add emmax-intel64 to PATH." >&2
  exit 1
fi

TFAM_FILE="${TFAM_FILE:-${TPED_PREFIX}.tfam}"
RUN_EMMAX="${RUN_EMMAX:-false}"
JOBS="${JOBS:-1}"
PLOT_RESULTS="${PLOT_RESULTS:-true}"
P_THRESHOLD="${P_THRESHOLD:-0.001}"
PLOT_SCRIPT="${PLOT_SCRIPT:-${script_dir}/plot_manhattan_qq_reference.R}"
PHENO_MISSING_VALUE="${PHENO_MISSING_VALUE:-NA}"

[[ -f "$TFAM_FILE" ]] || { echo "TFAM file not found: $TFAM_FILE" >&2; exit 1; }
[[ -f "$KINSHIP_FILE" ]] || { echo "Kinship file not found: $KINSHIP_FILE" >&2; exit 1; }
[[ -f "$PHENO_TSV" ]] || { echo "Phenotype matrix not found: $PHENO_TSV" >&2; exit 1; }
if [[ -n "${COVARIATE_TSV:-}" ]]; then
  [[ -f "$COVARIATE_TSV" ]] || { echo "Covariate file not found: $COVARIATE_TSV" >&2; exit 1; }
fi

mkdir -p "$OUTDIR/traits" "$OUTDIR/assoc" "$OUTDIR/plots" "$OUTDIR/logs"

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

command_file="$OUTDIR/run_emmax_commands.sh"
printf '%s\n' '#!/usr/bin/env bash' 'set -euo pipefail' > "$command_file"

manifest_file="$OUTDIR/trait_manifest.tsv"
printf 'trait_id\tphenotype_file\n' > "$manifest_file"

while IFS= read -r trait; do
  trait="${trait//$'\r'/}"
  [[ -n "$trait" ]] || continue
  python3 - "$PHENO_TSV" "$TFAM_FILE" "$OUTDIR/traits/${trait}.pheno.tsv" "$trait" "$PHENO_MISSING_VALUE" <<'PY'
import csv
import sys

pheno_tsv, tfam_file, out_file, trait, missing_value = sys.argv[1:6]
values = {}

with open(pheno_tsv, "r", encoding="utf-8", errors="replace") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    if "sample_id" not in reader.fieldnames:
        raise SystemExit("Phenotype matrix must contain a sample_id column")
    if trait not in reader.fieldnames:
        raise SystemExit(f"Trait not found in phenotype matrix: {trait}")
    for row in reader:
        sample = row["sample_id"].strip()
        value = row[trait].strip()
        values[sample] = value if value != "" else "NA"

with open(tfam_file, "r", encoding="utf-8", errors="replace") as tfam, open(out_file, "w", encoding="utf-8") as out:
    for line in tfam:
        line = line.strip()
        if not line:
            continue
        fields = line.split()
        fid, iid = fields[0], fields[1]
        out.write(f"{fid}\t{iid}\t{values.get(iid, missing_value)}\n")
PY

  printf '%s\t%s\n' "$trait" "$OUTDIR/traits/${trait}.pheno.tsv" >> "$manifest_file"

  if [[ -n "${COVARIATE_TSV:-}" ]]; then
    printf '%q ' "$EMMAX_BIN" -t "$TPED_PREFIX" -p "$OUTDIR/traits/${trait}.pheno.tsv" -k "$KINSHIP_FILE" -c "$COVARIATE_TSV" -o "$OUTDIR/assoc/${trait}" >> "$command_file"
  else
    printf '%q ' "$EMMAX_BIN" -t "$TPED_PREFIX" -p "$OUTDIR/traits/${trait}.pheno.tsv" -k "$KINSHIP_FILE" -o "$OUTDIR/assoc/${trait}" >> "$command_file"
  fi
  printf '\n' >> "$command_file"
done < "$OUTDIR/trait_list.txt"

chmod +x "$command_file"

cat > "$OUTDIR/plot_plan.txt" <<EOF
Expected per-trait summary columns:
trait	chrom	pos	snp_id	pvalue	beta	stderr

Expected visualization outputs:
plots/<trait>.manhattan.png
plots/<trait>.qq.png
assoc/<trait>.lead_snps.tsv

Reference plotting implementation:
skills/gwas-run-visualization/scripts/plot_manhattan_qq_reference.R

Reference style notes:
skills/gwas-run-visualization/references/plot_style.md

Resolved EMMAX binary:
${EMMAX_BIN}
EOF

if [[ "$RUN_EMMAX" == "true" ]]; then
  export http_proxy="${http_proxy:-http://127.0.0.1:7890}"
  export https_proxy="${https_proxy:-http://127.0.0.1:7890}"
  grep -v '^[[:space:]]*$' "$command_file" | tail -n +3 | xargs -I{} -P "$JOBS" bash -lc "{}"

  while IFS= read -r trait; do
    trait="${trait//$'\r'/}"
    [[ -n "$trait" ]] || continue
    ps_file="$OUTDIR/assoc/${trait}.ps"
    [[ -f "$ps_file" ]] || { echo "Missing EMMAX output: $ps_file" >&2; exit 1; }

    awk -F '[:\t]' 'BEGIN{OFS="\t"; print "SNP\tCHR\tBP\tP"} $5 ~ /^[0-9eE.+-]+$/ {print $1 ":" $2, $1, $2, $5}' "$ps_file" > "$OUTDIR/assoc/${trait}.all.QQ.pmap.txt"
    awk -F '[:\t]' -v p_threshold="$P_THRESHOLD" 'BEGIN{OFS="\t"; print "SNP\tCHR\tBP\tP"} $5 ~ /^[0-9eE.+-]+$/ && ($5 + 0) < p_threshold {print $1 ":" $2, $1, $2, $5}' "$ps_file" > "$OUTDIR/assoc/${trait}.sig.Manhattan.pmap.txt"
    awk -F '[:\t]' -v trait="$trait" 'BEGIN{OFS="\t"; print "trait\tchrom\tpos\tsnp_id\tpvalue\tbeta\tstderr"} $5 ~ /^[0-9eE.+-]+$/ {print trait, $1, $2, $1 ":" $2, $5, $3, $4}' "$ps_file" > "$OUTDIR/assoc/${trait}.gwas_summary.tsv"

    if [[ "$PLOT_RESULTS" == "true" ]]; then
      Rscript "$PLOT_SCRIPT" \
        --prefix "$trait" \
        --sig-file "$OUTDIR/assoc/${trait}.sig.Manhattan.pmap.txt" \
        --qq-file "$OUTDIR/assoc/${trait}.all.QQ.pmap.txt" \
        --outdir "$OUTDIR/plots"
    fi
  done < "$OUTDIR/trait_list.txt"
fi

printf 'Prepared %s, %s, and %s\n' "$OUTDIR/trait_list.txt" "$manifest_file" "$command_file"
