#!/usr/bin/env bash
# mk_demo_vcfs.sh
set -euo pipefail

START_POS=${1:-14370}
STEP=${2:-50}
VARS_PER_FILE=${VARS_PER_FILE:-3}
CHROM=20

mkdir -p demo_vcfs

mk() {
  local out=$1 sample=$2
  printf '%s' $'##fileformat=VCFv4.2\n##source=MetaGapDemo\n##reference=GRCh37\n##contig=<ID=20>\n##FILTER=<ID=q10,Description="Quality below 10">\n##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">\n##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\n' > "$out"
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' "$sample" >> "$out"
}

write_variants() {
  local out=$1 sample=$2 recs_name=$3
  local -n recs=$recs_name
  if [[ ${#recs[@]} -ne VARS_PER_FILE ]]; then
    echo "Expected ${VARS_PER_FILE} records for $sample but found ${#recs[@]}" >&2
    exit 1
  fi

  mk "$out" "$sample"
  for ((i = 0; i < VARS_PER_FILE; i++)); do
    local pos=$((START_POS + i * STEP))
    IFS='|' read -r id ref alt qual filter info format sample_data <<<"${recs[$i]}"
    printf '%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$CHROM" "$pos" "$id" "$ref" "$alt" "$qual" "$filter" "$info" "$format" "$sample_data" >> "$out"
  done
}

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

write_variants demo_vcfs/1.vcf SAMPLE_A records_a
write_variants demo_vcfs/2.vcf SAMPLE_B records_b
write_variants demo_vcfs/3.vcf SAMPLE_C records_c
write_variants demo_vcfs/4.vcf SAMPLE_D records_d

if command -v bgzip >/dev/null 2>&1 && command -v tabix >/dev/null 2>&1; then
  for f in demo_vcfs/*.vcf; do bgzip -f "$f"; tabix -f -p vcf "$f.gz"; done
  echo "Wrote: demo_vcfs/{1..4}.vcf.gz"
else
  echo "bgzip/tabix not found; skipping compression" >&2
  echo "Wrote: demo_vcfs/{1..4}.vcf"
fi
