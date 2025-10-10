#!/usr/bin/env bash
set -euo pipefail

# Tools
for t in bcftools bgzip tabix; do
  command -v "$t" >/dev/null || { echo "Missing tool: $t"; exit 1; }
done

# Usage
IN_DIR="${1:-}"; META_FILE="${2:-}"; OUT_DIR="${3:-}"
if [[ -z "$IN_DIR" || -z "$META_FILE" || -z "$OUT_DIR" ]]; then
  echo "Usage: $0 <input_vcf_directory> <metadata_template_file> <output_directory>"
  exit 1
fi
mkdir -p "$OUT_DIR"

# Tunables (override via env)
QUAL_THRESHOLD="${QUAL_THRESHOLD:-30}"
AN_THRESHOLD="${AN_THRESHOLD:-50}"

# Paths
MERGED_VCF="$OUT_DIR/merged_raw.vcf.gz"
JOINT_VCF="$OUT_DIR/cohort_joint.vcf.gz"
FILTERED_VCF="$OUT_DIR/cohort_filtered.vcf.gz"
ANON_VCF="$OUT_DIR/cohort_anon.vcf.gz"
FINAL_VCF_GZ="$OUT_DIR/cohort_final.vcf.gz"
META_CLEAN="$OUT_DIR/meta_clean.hdr"

echo "Normalizing metadata header..."
# Ensure one '##' directive per line, strip CRLF, and drop duplicate fileformat if present
tr -d '\r' < "$META_FILE" \
  | sed 's/##/\n##/g' \
  | sed '/^\s*$/d' \
  | grep -v '^##fileformat' > "$META_CLEAN"

echo "Collecting and preparing VCFs from: $IN_DIR"
mapfile -d '' -t VCF_LIST < <(find "$IN_DIR" -type f \( -name "*.vcf" -o -name "*.vcf.gz" \) -print0)
(( ${#VCF_LIST[@]} )) || { echo "ERROR: No VCF files found in $IN_DIR"; exit 1; }

PREPARED=()
for f in "${VCF_LIST[@]}"; do
  g="$f"
  if [[ "$f" == *.vcf ]]; then bgzip -f "$f"; g="$f.gz"; fi
  tabix -f -p vcf "$g"
  bcftools view -h "$g" >/dev/null
  PREPARED+=("$g")
done

echo "Merging..."
bcftools merge --no-version -m all -Oz -o "$MERGED_VCF" "${PREPARED[@]}"

echo "Recomputing AC/AN/AF..."
bcftools +fill-tags --no-version "$MERGED_VCF" -Oz -o "$JOINT_VCF" -- -t AC,AN,AF

echo "Filtering (PASS, QUAL>${QUAL_THRESHOLD}, AN>=${AN_THRESHOLD})..."
bcftools view --no-version -f PASS -i "QUAL>${QUAL_THRESHOLD} & INFO/AN>=${AN_THRESHOLD}" -Oz -o "$FILTERED_VCF" "$JOINT_VCF"

echo "Removing genotypes and sample columns..."
bcftools view --no-version -G -Oz -o "$ANON_VCF" "$FILTERED_VCF"

echo "Appending custom metadata..."
bcftools annotate --no-version -h "$META_CLEAN" -Oz -o "$FINAL_VCF_GZ" "$ANON_VCF"
tabix -f -p vcf "$FINAL_VCF_GZ"

# Cleanup
rm -f "$MERGED_VCF" "$JOINT_VCF" "$FILTERED_VCF" "$ANON_VCF" "$META_CLEAN"

echo "Done. Final joint-called VCF: $FINAL_VCF_GZ"