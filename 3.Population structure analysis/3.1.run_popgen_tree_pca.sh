#!/bin/bash
set -euo pipefail

# software
vcftools="/project/data/tool/vcftools"
bcftools="/project/data/tool/bcftools"
tassel="/project/data/tool/TASSEL5/run_pipeline.pl"
fasttree="/project/data/tool/FastTree"
plink="/project/data/tool/plink"

# input/output
INVCF="/project/data/gatk_call/vcf/filtered/all.pass.snp.vcf.gz"
OUTDIR="/project/data/popgen"
PREFIX="rice.maf05.miss95"

mkdir -p "${OUTDIR}"

# 1. SNP filtering
"${vcftools}" --gzvcf "${INVCF}" \
    --maf 0.05 \
    --max-missing 0.95 \
    --recode --recode-INFO-all \
    --out "${OUTDIR}/${PREFIX}"

bgzip -f "${OUTDIR}/${PREFIX}.recode.vcf"
tabix -f -p vcf "${OUTDIR}/${PREFIX}.recode.vcf.gz"

# 2. sort VCF
"${bcftools}" sort "${OUTDIR}/${PREFIX}.recode.vcf.gz" \
    -Oz -o "${OUTDIR}/${PREFIX}.sorted.vcf.gz"
tabix -f -p vcf "${OUTDIR}/${PREFIX}.sorted.vcf.gz"

# 3. VCF -> PHYLIP for phylogeny
gunzip -c "${OUTDIR}/${PREFIX}.sorted.vcf.gz" > "${OUTDIR}/${PREFIX}.sorted.vcf"

"${tassel}" -Xms4g -Xmx32g \
    -importGuess "${OUTDIR}/${PREFIX}.sorted.vcf" \
    -ExportPlugin \
    -saveAs "${OUTDIR}/${PREFIX}.phy" \
    -format Phylip_Inter

# 4. maximum-likelihood tree with FastTree
"${fasttree}" -gtr -nt "${OUTDIR}/${PREFIX}.phy" > "${OUTDIR}/${PREFIX}.fasttree.nwk"

# 5. PCA in PLINK
"${plink}" \
    --vcf "${OUTDIR}/${PREFIX}.sorted.vcf.gz" \
    --double-id \
    --allow-extra-chr \
    --set-missing-var-ids @:# \
    --pca 10 \
    --out "${OUTDIR}/${PREFIX}.pca"

echo "Done."