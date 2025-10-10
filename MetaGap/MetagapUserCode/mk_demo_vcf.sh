#!/usr/bin/env bash
# mk_demo_vcfs.sh
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: mk_demo_vcf.sh [OPTIONS] [START_POS [STEP]]

Generate a set of single-sample VCFs for exercising the merge pipeline.

Options:
  --mode MODE      Generation mode. Supported values:
                     shared  - Every file contains the same variants (default).
                     partial - Files share some variants but also contain
                                sample-specific or subset-only records.
  --start-pos POS  Override the starting position for the first variant.
  --step STEP      Override the distance between consecutive variants.
  -h, --help       Display this help message and exit.

The legacy positional arguments START_POS and STEP are still accepted for
compatibility. Environment variables START_POS, STEP, MODE, and VARS_PER_FILE
may also be set prior to invocation.
USAGE
}

START_POS=${START_POS:-14370}
STEP=${STEP:-50}

if [[ -z "${VARS_PER_FILE+x}" ]]; then
  VARS_PER_FILE=3
  VARS_PER_FILE_WAS_DEFAULT=1
else
  VARS_PER_FILE_WAS_DEFAULT=0
fi

MODE=${MODE:-shared}
CHROM=20

declare -a POSITIONAL_ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode)
      [[ $# -ge 2 ]] || { echo "--mode requires an argument" >&2; exit 1; }
      MODE=$2
      shift 2
      ;;
    --start-pos)
      [[ $# -ge 2 ]] || { echo "--start-pos requires an argument" >&2; exit 1; }
      START_POS=$2
      shift 2
      ;;
    --step)
      [[ $# -ge 2 ]] || { echo "--step requires an argument" >&2; exit 1; }
      STEP=$2
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --)
      shift
      POSITIONAL_ARGS+=("$@")
      break
      ;;
    -*)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 1
      ;;
    *)
      POSITIONAL_ARGS+=("$1")
      shift
      ;;
  esac
done

if [[ ${#POSITIONAL_ARGS[@]} -gt 0 ]]; then
  START_POS=${POSITIONAL_ARGS[0]}
fi

if [[ ${#POSITIONAL_ARGS[@]} -gt 1 ]]; then
  STEP=${POSITIONAL_ARGS[1]}
fi

if [[ ${#POSITIONAL_ARGS[@]} -gt 2 ]]; then
  echo "Too many positional arguments" >&2
  usage >&2
  exit 1
fi

case "$MODE" in
  shared|partial)
    ;;
  *)
    echo "Unsupported mode: $MODE" >&2
    usage >&2
    exit 1
    ;;
esac

mkdir -p demo_vcfs

if ! [[ $SAMPLE_COUNT =~ ^[0-9]+$ ]] || (( SAMPLE_COUNT <= 0 )); then
  echo "Sample count must be a positive integer (received: $SAMPLE_COUNT)" >&2
  exit 1
fi

REFS=(G T A C)
ALTS=(A C G T)
GENOTYPES=(0/0 0/1 1/1 ./.)

mk() {
  local out=$1 sample=$2
  printf '%s' $'##fileformat=VCFv4.2\n##source=MetaGapDemo\n##reference=GRCh37\n##contig=<ID=20>\n##FILTER=<ID=q10,Description="Quality below 10">\n##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">\n##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\n' > "$out"
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' "$sample" >> "$out"
}

write_variants() {
  local out=$1 sample=$2 recs_name=$3
  local -n recs=$recs_name
  if [[ ${#recs[@]} -gt VARS_PER_FILE ]]; then
    echo "Expected at most ${VARS_PER_FILE} records for $sample but found ${#recs[@]}" >&2
    exit 1
generate_variant_fields() {
  local sample_idx=$1 variant_idx=$2

  local id="rs$((6054257 + sample_idx * VARS_PER_FILE + variant_idx))"
  local ref=${REFS[variant_idx % ${#REFS[@]}]}
  local alt=${ALTS[((sample_idx + variant_idx) % ${#ALTS[@]})]}
  if [[ $alt == $ref ]]; then
    alt=${ALTS[((sample_idx + variant_idx + 1) % ${#ALTS[@]})]}
  fi

  local qual=$((45 + ((sample_idx * 11 + variant_idx * 7) % 30)))
  local filter=PASS
  if (( (sample_idx + variant_idx) % 5 == 0 )); then
    filter=q10
  fi

  local dp=$((8 + ((sample_idx * 5 + variant_idx * 3) % 12)))
  local info="NS=1;DP=$dp"
  local format="GT:DP:GQ"

  local gt=${GENOTYPES[((sample_idx + variant_idx) % ${#GENOTYPES[@]})]}
  local gq=$((35 + ((sample_idx * 3 + variant_idx * 5) % 40)))
  local sample_data="$gt:$dp:$gq"

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
    "$id" "$ref" "$alt" "$qual" "$filter" "$info" "$format" "$sample_data"
}

write_variants() {
  local out=$1 sample=$2 sample_idx=$3

  mk "$out" "$sample"
  for ((i = 0; i < VARS_PER_FILE; i++)); do
    local rec="${recs[$i]-}"
    if [[ -z "$rec" ]]; then
      continue
    fi
    local pos=$((START_POS + i * STEP))
    IFS=$'\t' read -r id ref alt qual filter info format sample_data <<<"$(generate_variant_fields "$sample_idx" "$i")"
    printf '%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$CHROM" "$pos" "$id" "$ref" "$alt" "$qual" "$filter" "$info" "$format" "$sample_data" >> "$out"
  done
}

for ((sample_idx = 1; sample_idx <= SAMPLE_COUNT; sample_idx++)); do
  sample_name=$(printf 'SAMPLE_%03d' "$sample_idx")
  out_file="demo_vcfs/${sample_idx}.vcf"
  write_variants "$out_file" "$sample_name" "$sample_idx"
done

if command -v bgzip >/dev/null 2>&1 && command -v tabix >/dev/null 2>&1; then
  for f in demo_vcfs/*.vcf; do bgzip -f "$f"; tabix -f -p vcf "$f.gz"; done
  echo "Wrote $SAMPLE_COUNT compressed sample VCF files in demo_vcfs/"
else
  echo "bgzip/tabix not found; skipping compression" >&2
  echo "Wrote $SAMPLE_COUNT sample VCF files in demo_vcfs/"
fi
