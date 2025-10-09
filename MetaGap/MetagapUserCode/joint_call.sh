#!/usr/bin/env bash
set -euo pipefail

# activate project env
if ! command -v conda >/dev/null 2>&1; then
  echo "conda not found in PATH"; exit 1
fi
# shell integration for non-interactive shells
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$(pwd)/.conda/env"

# hard fail if tools missing
for t in bcftools bgzip tabix; do
  command -v "$t" >/dev/null || { echo "Missing tool: $t"; exit 1; }
done


# Minimal Joint Calling Pipeline Script
# Usage: joint_call_mvp.sh <input_vcf_dir> <metadata_template.txt> <output_dir>

# Input arguments
IN_DIR="${1:-}"
META_FILE="${2:-}"
OUT_DIR="${3:-}"

if [[ -z "$IN_DIR" || -z "$META_FILE" || -z "$OUT_DIR" ]]; then
    echo "Usage: $0 <input_vcf_directory> <metadata_template_file> <output_directory>"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# Temporary file paths
MERGED_VCF="$OUT_DIR/merged_raw.vcf.gz"
JOINT_VCF="$OUT_DIR/cohort_joint.vcf.gz"
FILTERED_VCF="$OUT_DIR/cohort_filtered.vcf.gz"
FINAL_VCF="$OUT_DIR/cohort_final.vcf"
FINAL_VCF_GZ="$FINAL_VCF.gz"

echo "Collecting VCF files from input directory: $IN_DIR"
# Find all .vcf or .vcf.gz files in the input directory
VCF_LIST=()
while IFS= read -r -d $'\0' file; do
    VCF_LIST+=("$file")
done < <(find "$IN_DIR" -type f \( -name "*.vcf" -o -name "*.vcf.gz" \) -print0)

if [[ "${#VCF_LIST[@]}" -eq 0 ]]; then
    echo "ERROR: No VCF files found in $IN_DIR"
    exit 1
fi

echo "Ensuring all VCFs are bgzipped and indexed..."
# Build a fresh list of gzipped, indexed inputs
PREPARED_LIST=()
for f in "${VCF_LIST[@]}"; do
  g="$f"
  if [[ "$f" == *.vcf ]]; then
    bgzip -f "$f"
    g="$f.gz"
  fi
  tabix -f -p vcf "$g"
  PREPARED_LIST+=("$g")
done
VCF_LIST=("${PREPARED_LIST[@]}")

echo "Check *.vcf files"
for g in "${VCF_LIST[@]}"; do
  bcftools view -h "$g" >/dev/null || exit 1
done

echo "Found ${#VCF_LIST[@]} VCF files. Merging them into a combined VCF..."
# Merge all VCFs, preserving all alleles/genotypes (-m all)
# Output as compressed VCF to MERGED_VCF
bcftools merge -m all -Oz -o "$MERGED_VCF" "${VCF_LIST[@]}"

echo "Recalculating AC, AN, AF tags across the merged cohort..."
# Use bcftools +fill-tags to recalc AC, AN, AF in INFO for the merged file
bcftools +fill-tags "$MERGED_VCF" -Oz -o "$JOINT_VCF" -- -t AC,AN,AF

echo "Applying quality filters to remove low-confidence variants..."
# Filter variants: QUAL > 30, total alleles > 50, and PASS only
bcftools view -i 'QUAL>30 && INFO/AN>50 && FILTER="PASS"' -Oz -o "$FILTERED_VCF" "$JOINT_VCF"

echo "Removing individual sample genotype columns (anonymizing data)..."
# Extract header and data without sample columns.
#  -H: print data lines only (no header), then cut first 8 columns (up to INFO).
bcftools view -h "$FILTERED_VCF" > "$OUT_DIR/header.temp"
bcftools view -H "$FILTERED_VCF" | cut -f1-8 > "$OUT_DIR/body.temp"

# Prepare final VCF header with custom metadata
echo "Combining custom metadata header with VCF header..."
{
    # Print original fileformat line first
    grep -m1 '^##fileformat' "$OUT_DIR/header.temp"
    # Print custom metadata lines
    cat "$META_FILE"
    # Print the rest of the original header lines (skip fileformat and column header)
    grep -v '^##fileformat' "$OUT_DIR/header.temp" | grep -v '^#CHROM'
    # Print the final column header line (CHROM...INFO only, no FORMAT or samples)
    echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
} > "$OUT_DIR/header_final.temp"

# Combine header and body to create final VCF
cat "$OUT_DIR/header_final.temp" "$OUT_DIR/body.temp" > "$FINAL_VCF"

echo "Compressing and indexing the final VCF..."
bgzip -f "$FINAL_VCF"   # compress to .vcf.gz
tabix -p vcf -f "$FINAL_VCF_GZ"  # index the compressed VCF

# Clean up temporary files
rm -f "$OUT_DIR/header.temp" "$OUT_DIR/body.temp" "$OUT_DIR/header_final.temp" "$MERGED_VCF" "$JOINT_VCF" "$FILTERED_VCF"

echo "Done. Final joint-called VCF: $FINAL_VCF_GZ"