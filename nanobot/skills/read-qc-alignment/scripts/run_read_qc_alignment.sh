#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: bash scripts/run_read_qc_alignment.sh --config path/to/read_qc_alignment.env
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

required_vars=(SAMPLE_SHEET FASTQ_DIR REFERENCE_FASTA OUTDIR THREADS ALIGNER)
for var_name in "${required_vars[@]}"; do
  [[ -n "${!var_name:-}" ]] || { echo "Missing variable: $var_name" >&2; exit 1; }
done

mkdir -p "$OUTDIR/qc_raw" "$OUTDIR/qc_trimmed" "$OUTDIR/trimmed" "$OUTDIR/bam" "$OUTDIR/stats"

awk -F'\t' '
NR==1 {
  if ($1 != "sample_id" || $2 != "read1" || $3 != "read2") {
    print "sample_sheet header must be: sample_id read1 read2" > "/dev/stderr"
    exit 2
  }
}
NR>1 {
  print $1 "\t" $2 "\t" $3
}
' "$SAMPLE_SHEET" > "$OUTDIR/sample_manifest.tsv"

cat > "$OUTDIR/workflow_plan.txt" <<EOF
fastqc -t ${THREADS} ${FASTQ_DIR}/<read1> ${FASTQ_DIR}/<read2> -o ${OUTDIR}/qc_raw
fastp -i ${FASTQ_DIR}/<read1> -I ${FASTQ_DIR}/<read2> -o ${OUTDIR}/trimmed/<sample>_R1.fq.gz -O ${OUTDIR}/trimmed/<sample>_R2.fq.gz
${ALIGNER} mem -t ${THREADS} ${REFERENCE_FASTA} ${OUTDIR}/trimmed/<sample>_R1.fq.gz ${OUTDIR}/trimmed/<sample>_R2.fq.gz | samtools sort -@ ${THREADS} -o ${OUTDIR}/bam/<sample>.sorted.bam
samtools index ${OUTDIR}/bam/<sample>.sorted.bam
samtools flagstat ${OUTDIR}/bam/<sample>.sorted.bam > ${OUTDIR}/stats/<sample>.flagstat.txt
EOF

printf 'Prepared %s and %s\n' "$OUTDIR/sample_manifest.tsv" "$OUTDIR/workflow_plan.txt"
