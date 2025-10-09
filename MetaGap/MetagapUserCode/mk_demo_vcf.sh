# mk_demo_vcfs.sh
set -euo pipefail
mkdir -p demo_vcfs

mk() {
  out=$1; sample=$2
  printf '%s' $'##fileformat=VCFv4.2\n##source=MetaGapDemo\n##reference=GRCh37\n##contig=<ID=20>\n##FILTER=<ID=q10,Description="Quality below 10">\n##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">\n##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\n' > "$out"
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' "$sample" >> "$out"
}

# A
mk demo_vcfs/1.vcf SAMPLE_A
printf $'20\t14370\trs6054257\tG\tA\t60\tPASS\tNS=1;DP=10\tGT:DP:GQ\t0/1:10:50\n20\t17330\t.\tT\tA\t40\tPASS\tNS=1;DP=12\tGT:DP:GQ\t0/0:12:60\n20\t1110696\trs6040355\tA\tG\t50\tPASS\tNS=1;DP=8\tGT:DP:GQ\t1/1:8:40\n' >> demo_vcfs/1.vcf

# B
mk demo_vcfs/2.vcf SAMPLE_B
printf $'20\t14370\trs6054257\tG\tA\t44\tPASS\tNS=1;DP=9\tGT:DP:GQ\t0/0:9:50\n20\t17330\t.\tT\tA\t35\tPASS\tNS=1;DP=11\tGT:DP:GQ\t0/1:11:45\n20\t1110696\trs6040355\tA\tG\t55\tPASS\tNS=1;DP=7\tGT:DP:GQ\t0/1:7:48\n' >> demo_vcfs/2.vcf

# C
mk demo_vcfs/3.vcf SAMPLE_C
printf $'20\t14370\trs6054257\tG\tA\t70\tPASS\tNS=1;DP=12\tGT:DP:GQ\t1/1:12:60\n20\t17330\t.\tT\tA\t42\tPASS\tNS=1;DP=10\tGT:DP:GQ\t0/0:10:55\n20\t1110696\trs6040355\tA\tG\t52\tPASS\tNS=1;DP=9\tGT:DP:GQ\t0/1:9:52\n' >> demo_vcfs/3.vcf

# D
mk demo_vcfs/4.vcf SAMPLE_D
printf $'20\t14370\trs6054257\tG\tA\t50\tPASS\tNS=1;DP=0\tGT:DP:GQ\t./.:.:.\n20\t17330\t.\tT\tA\t60\tPASS\tNS=1;DP=13\tGT:DP:GQ\t1/1:13:62\n20\t1110696\trs6040355\tA\tG\t45\tPASS\tNS=1;DP=8\tGT:DP:GQ\t0/0:8:45\n' >> demo_vcfs/4.vcf

# compress + index
for f in demo_vcfs/*.vcf; do bgzip -f "$f"; tabix -f -p vcf "$f.gz"; done
echo "Wrote: demo_vcfs/{1..4}.vcf.gz"
