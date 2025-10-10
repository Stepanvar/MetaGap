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
  fi

  mk "$out" "$sample"
  for ((i = 0; i < VARS_PER_FILE; i++)); do
    local rec="${recs[$i]-}"
    if [[ -z "$rec" ]]; then
      continue
    fi
    local pos=$((START_POS + i * STEP))
    IFS='|' read -r id ref alt qual filter info format sample_data <<<"$rec"
    printf '%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$CHROM" "$pos" "$id" "$ref" "$alt" "$qual" "$filter" "$info" "$format" "$sample_data" >> "$out"
  done
}

if [[ $MODE == shared ]]; then
  records_a=(
    'rs6054257|G|A|60|PASS|NS=1;DP=10|GT:DP:GQ|0/1:10:50'
    '.|T|A|40|PASS|NS=1;DP=12|GT:DP:GQ|0/0:12:60'
    'rs6040355|A|G|50|PASS|NS=1;DP=8|GT:DP:GQ|1/1:8:40'
  )

  records_b=(
    'rs6054257|G|A|44|PASS|NS=1;DP=9|GT:DP:GQ|0/0:9:50'
    '.|T|A|35|PASS|NS=1;DP=11|GT:DP:GQ|0/1:11:45'
    'rs6040355|A|G|55|PASS|NS=1;DP=7|GT:DP:GQ|0/1:7:48'
  )

  records_c=(
    'rs6054257|G|A|70|PASS|NS=1;DP=12|GT:DP:GQ|1/1:12:60'
    '.|T|A|42|PASS|NS=1;DP=10|GT:DP:GQ|0/0:10:55'
    'rs6040355|A|G|52|PASS|NS=1;DP=9|GT:DP:GQ|0/1:9:52'
  )

  records_d=(
    'rs6054257|G|A|50|PASS|NS=1;DP=0|GT:DP:GQ|./.:.:.'
    '.|T|A|60|PASS|NS=1;DP=13|GT:DP:GQ|1/1:13:62'
    'rs6040355|A|G|45|PASS|NS=1;DP=8|GT:DP:GQ|0/0:8:45'
  )
else
  if [[ $VARS_PER_FILE_WAS_DEFAULT -eq 1 ]]; then
    VARS_PER_FILE=4
  fi

  records_a=(
    'rs6054257|G|A|60|PASS|NS=1;DP=10|GT:DP:GQ|0/1:10:50'
    'rsSubsetAC|C|T|58|PASS|NS=1;DP=14|GT:DP:GQ|1/1:14:60'
    'rs6040355|A|G|50|PASS|NS=1;DP=8|GT:DP:GQ|1/1:8:40'
    ''
  )

  records_b=(
    'rs6054257|G|A|44|PASS|NS=1;DP=9|GT:DP:GQ|0/0:9:50'
    ''
    'rs6040355|A|G|55|PASS|NS=1;DP=7|GT:DP:GQ|0/1:7:48'
    ''
  )

  records_c=(
    'rs6054257|G|A|70|PASS|NS=1;DP=12|GT:DP:GQ|1/1:12:60'
    'rsSubsetAC|C|T|62|PASS|NS=1;DP=13|GT:DP:GQ|0/1:13:55'
    ''
    ''
  )

  records_d=(
    'rs6054257|G|A|50|PASS|NS=1;DP=0|GT:DP:GQ|./.:.:.'
    ''
    'rs6040355|A|G|45|PASS|NS=1;DP=8|GT:DP:GQ|0/0:8:45'
    'rsUniqueD|T|C|60|PASS|NS=1;DP=11|GT:DP:GQ|0/1:11:62'
  )
fi

write_variants demo_vcfs/1.vcf SAMPLE_A records_a
write_variants demo_vcfs/2.vcf SAMPLE_B records_b
write_variants demo_vcfs/3.vcf SAMPLE_C records_c
write_variants demo_vcfs/4.vcf SAMPLE_D records_d

if command -v bgzip >/dev/null 2>&1 && command -v tabix >/dev/null 2>&1; then
  for f in demo_vcfs/*.vcf; do bgzip -f "$f"; tabix -f -p vcf "$f.gz"; done
  echo "Wrote: demo_vcfs/{1..4}.vcf.gz (mode=$MODE)"
else
  echo "bgzip/tabix not found; skipping compression" >&2
  echo "Wrote: demo_vcfs/{1..4}.vcf (mode=$MODE)"
fi
