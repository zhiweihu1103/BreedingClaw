#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: bash scripts/run_variant_call_filter.sh --config path/to/variant_call_filter.env
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

required_vars=(BAM_LIST REFERENCE_FASTA OUTDIR CALLER MIN_QUAL MIN_DP MAX_MISSING)
for var_name in "${required_vars[@]}"; do
  [[ -n "${!var_name:-}" ]] || { echo "Missing variable: $var_name" >&2; exit 1; }
done

mkdir -p "$OUTDIR/raw_calls" "$OUTDIR/normalized" "$OUTDIR/filtered" "$OUTDIR/stats"

awk -F'\t' '
NR==1 {
  if ($1 != "sample_id" || $2 != "bam") {
    print "bam_list header must be: sample_id bam" > "/dev/stderr"
    exit 2
  }
}
NR>1 {
  print $1 "\t" $2
}
' "$BAM_LIST" > "$OUTDIR/bam_manifest.tsv"

awk -F'\t' 'NR>1 { print $2 }' "$BAM_LIST" > "$OUTDIR/bam_paths.list"

case "$CALLER" in
  bcftools)
    cat > "$OUTDIR/workflow_plan.txt" <<EOF
bcftools mpileup -f ${REFERENCE_FASTA} -b ${OUTDIR}/bam_paths.list -Ou | bcftools call -mv -Oz -o ${OUTDIR}/raw_calls/cohort.raw.vcf.gz
bcftools norm -f ${REFERENCE_FASTA} -m -any ${OUTDIR}/raw_calls/cohort.raw.vcf.gz -Oz -o ${OUTDIR}/normalized/cohort.norm.vcf.gz
bcftools view -i 'QUAL>=${MIN_QUAL} && INFO/DP>=${MIN_DP}' ${OUTDIR}/normalized/cohort.norm.vcf.gz -Oz -o ${OUTDIR}/filtered/cohort.filtered.vcf.gz
bcftools stats ${OUTDIR}/filtered/cohort.filtered.vcf.gz > ${OUTDIR}/stats/cohort.filtered.stats.txt
EOF
    ;;
  gatk)
    cat > "$OUTDIR/workflow_plan.txt" <<EOF
gatk HaplotypeCaller -R ${REFERENCE_FASTA} -I <sample.bam> -O ${OUTDIR}/raw_calls/<sample>.g.vcf.gz -ERC GVCF
gatk CombineGVCFs -R ${REFERENCE_FASTA} --variant <sample1.g.vcf.gz> --variant <sample2.g.vcf.gz> -O ${OUTDIR}/raw_calls/cohort.g.vcf.gz
gatk GenotypeGVCFs -R ${REFERENCE_FASTA} -V ${OUTDIR}/raw_calls/cohort.g.vcf.gz -O ${OUTDIR}/raw_calls/cohort.raw.vcf.gz
bcftools norm -f ${REFERENCE_FASTA} -m -any ${OUTDIR}/raw_calls/cohort.raw.vcf.gz -Oz -o ${OUTDIR}/normalized/cohort.norm.vcf.gz
bcftools view -i 'QUAL>=${MIN_QUAL} && INFO/DP>=${MIN_DP}' ${OUTDIR}/normalized/cohort.norm.vcf.gz -Oz -o ${OUTDIR}/filtered/cohort.filtered.vcf.gz
EOF
    ;;
  *)
    echo "Unsupported CALLER: $CALLER" >&2
    exit 1
    ;;
esac

printf 'Prepared %s, %s, and %s\n' "$OUTDIR/bam_manifest.tsv" "$OUTDIR/bam_paths.list" "$OUTDIR/workflow_plan.txt"
