#!/usr/bin/env bash
set -euo pipefail

############################################
# Configuration
############################################
REF="hg38.fa"
CHR="chrM"
START=1000
END=15500
COVERAGE=2000          # increase to 2000+ if you want extreme
PREFIX="chrM_heavy_stress"

############################################
# Dependency check
############################################
for tool in samtools bwa bcftools art_illumina python3; do
    command -v $tool >/dev/null 2>&1 || { echo "$tool not found"; exit 1; }
done

############################################
# 1. Extract chrM
############################################
echo "Extracting chrM..."
samtools faidx $REF $CHR > ${CHR}.hg38.fa
samtools faidx ${CHR}.hg38.fa
bwa index ${CHR}.hg38.fa

############################################
# 2. Deterministic heavy mutation
############################################
echo "Generating dense deterministic variants..."

python3 <<EOF
start = $START
end = $END

def mutate(base):
    # deterministic base flip (no randomness)
    mapping = {"A":"C", "C":"G", "G":"T", "T":"A"}
    return mapping.get(base.upper(), "A")

with open("${CHR}.hg38.fa") as f:
    header = f.readline()
    seq = list("".join(line.strip() for line in f))

i = start
while i < end and i < len(seq):

    # SNP every 20 bp
    if i % 20 == 0:
        seq[i] = mutate(seq[i])

    i += 1

with open("${CHR}.mut.fa", "w") as out:
    out.write(header)
    for j in range(0, len(seq), 60):
        out.write("".join(seq[j:j+60]) + "\n")
EOF

############################################
# 3. Simulate ultra-high depth reads
############################################
echo "Simulating ${COVERAGE}x coverage..."

art_illumina \
    -ss HS25 \
    -i ${CHR}.mut.fa \
    -l 150 \
    -f ${COVERAGE} \
    -m 200 \
    -s 10 \
    -o ${PREFIX}

############################################
# 4. Align to original reference
############################################
echo "Aligning reads..."

bwa mem ${CHR}.hg38.fa \
    ${PREFIX}1.fq ${PREFIX}2.fq | \
    samtools sort -o ${PREFIX}.sorted.bam

samtools index ${PREFIX}.sorted.bam

############################################
# 5. Call variants
############################################
echo "Calling variants..."

bcftools mpileup --max-depth 100000 --no-BAQ -q 0 -Q 0 -f ${CHR}.hg38.fa ${PREFIX}.sorted.bam | \
bcftools call -mv --ploidy 1 -Oz -o ${PREFIX}.vcf.gz

tabix -p vcf ${PREFIX}.vcf.gz

############################################
echo ""
echo "✅ DONE"
echo "Variants densely distributed in:"
echo "  ${CHR}:${START}-${END}"
echo ""
echo "High-depth BAM:"
echo "  ${PREFIX}.sorted.bam"
echo ""
echo "Variant count:"
bcftools view -H ${PREFIX}.vcf.gz | wc -l

bgzip -d -f ${PREFIX}.vcf.gz

### Command for generating vcf for chr11 HG00096
# bcftools mpileup --max-depth 100000 --no-BAQ -q 0 -Q 0 -f human_g1k_v37.fasta HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam | \
# bcftools call -mv -Oz -o HG00096.chr11.vcf.gz
#
# tabix -p vcf HG00096.chr11.vcf.gz
#
# echo "Variant count:"
# bcftools view -H HG00096.chr11.vcf.gz | wc -l
#