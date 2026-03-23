#!/bin/bash
set -euo pipefail

# 软件路径
bwa="/project/data/tool/bwa-mem2"
samtools="/project/data/tool/samtools"
gatk="/project/data/tool/gatk"
picard="/project/data/tool/picard.jar"
snpeff="/project/data/tool/snpEff/snpEff.jar"
snpeff_cfg="/project/data/tool/snpEff/snpEff.config"

# 输入输出路径
REF="/project/data/ref/MSU7.fa"
FQ_DIR="/project/data/fastq"
BAM_DIR="/project/data/bam"
OUT_DIR="/project/data/gatk_call"
SAMPLE_LIST="/project/data/sample.list"
SNPEFF_DB="rice"

TMP="${OUT_DIR}/tmp"
GVCF="${OUT_DIR}/gvcf"
COMBINE="${OUT_DIR}/combine"
GENO="${OUT_DIR}/geno"
VCF="${OUT_DIR}/vcf"

CHRS=(Chr1 Chr2 Chr3 Chr4 Chr5 Chr6 Chr7 Chr8 Chr9 Chr10 Chr11 Chr12)

mkdir -p "$BAM_DIR" "$TMP" "$GVCF" "$COMBINE" "$GENO" "$VCF/split" "$VCF/filtered"
for chr in "${CHRS[@]}"; do
    mkdir -p "$GVCF/$chr"
done

# 1) 参考基因组建索引
[[ -f "${REF}.0123" || -f "${REF}.bwt.2bit.64" || -f "${REF}.bwt" ]] || "$bwa" index "$REF"
[[ -f "${REF}.fai" ]] || "$samtools" faidx "$REF"
[[ -f "${REF%.*}.dict" ]] || "$gatk" CreateSequenceDictionary -R "$REF"

# 2) 比对、排序、去重复、建 BAM 索引
while read -r s; do
    [[ -z "$s" ]] && continue

    sorted="$BAM_DIR/$s.sorted.bam"
    markdup="$BAM_DIR/$s.sorted.markdup.bam"

    if [[ ! -f "$markdup" ]]; then
        "$bwa" mem -t 24 \
            -R "@RG\tID:${s}\tSM:${s}\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
            "$REF" "$FQ_DIR/${s}_R1.fq.gz" "$FQ_DIR/${s}_R2.fq.gz" | \
        "$samtools" sort -@ 16 -o "$sorted" -

        java -jar "$picard" MarkDuplicates \
            I="$sorted" \
            O="$markdup" \
            M="$BAM_DIR/$s.markdup.metrics.txt" \
            CREATE_INDEX=true \
            REMOVE_DUPLICATES=false \
            TMP_DIR="$TMP"

        rm -f "$sorted"
    fi

    # 3) 每个样品、每条染色体生成 gVCF
    for chr in "${CHRS[@]}"; do
        out="$GVCF/$chr/$s.$chr.g.vcf.gz"
        [[ -f "$out" ]] && continue

        "$gatk" --java-options "-Xmx12g -Djava.io.tmpdir=$TMP" HaplotypeCaller \
            -R "$REF" \
            -I "$markdup" \
            -L "$chr" \
            -ERC GVCF \
            -O "$out"
    done
done < "$SAMPLE_LIST"

# 4) 每条染色体联合分型
for chr in "${CHRS[@]}"; do
    args=()
    while read -r s; do
        [[ -z "$s" ]] && continue
        args+=("-V" "$GVCF/$chr/$s.$chr.g.vcf.gz")
    done < "$SAMPLE_LIST"

    "$gatk" --java-options "-Xmx20g -Djava.io.tmpdir=$TMP" CombineGVCFs \
        -R "$REF" \
        "${args[@]}" \
        -O "$COMBINE/$chr.g.vcf.gz"

    "$gatk" --java-options "-Xmx20g -Djava.io.tmpdir=$TMP" GenotypeGVCFs \
        -R "$REF" \
        -V "$COMBINE/$chr.g.vcf.gz" \
        -O "$GENO/$chr.raw.vcf.gz"
done

# 5) 合并全基因组 VCF
merge_args=()
for chr in "${CHRS[@]}"; do
    merge_args+=("-I" "$GENO/$chr.raw.vcf.gz")
done

"$gatk" --java-options "-Xmx20g -Djava.io.tmpdir=$TMP" MergeVcfs \
    "${merge_args[@]}" \
    -O "$VCF/all.raw.vcf.gz"

# 6) 拆分 SNP 和 INDEL
"$gatk" --java-options "-Xmx12g -Djava.io.tmpdir=$TMP" SelectVariants \
    -R "$REF" \
    -V "$VCF/all.raw.vcf.gz" \
    --select-type-to-include SNP \
    -O "$VCF/split/all.raw.snp.vcf.gz"

"$gatk" --java-options "-Xmx12g -Djava.io.tmpdir=$TMP" SelectVariants \
    -R "$REF" \
    -V "$VCF/all.raw.vcf.gz" \
    --select-type-to-include INDEL \
    -O "$VCF/split/all.raw.indel.vcf.gz"

# 7) 硬过滤
"$gatk" --java-options "-Xmx12g -Djava.io.tmpdir=$TMP" VariantFiltration \
    -R "$REF" \
    -V "$VCF/split/all.raw.snp.vcf.gz" \
    -O "$VCF/filtered/all.raw.snp.filtered.vcf.gz" \
    --filter-name "SNP_HardFilter" \
    --filter-expression "QD < 2.0 || QUAL < 30.0 || SOR > 3.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"

"$gatk" --java-options "-Xmx12g -Djava.io.tmpdir=$TMP" VariantFiltration \
    -R "$REF" \
    -V "$VCF/split/all.raw.indel.vcf.gz" \
    -O "$VCF/filtered/all.raw.indel.filtered.vcf.gz" \
    --filter-name "INDEL_HardFilter" \
    --filter-expression "QD < 2.0 || QUAL < 30.0 || FS > 200.0 || ReadPosRankSum < -20.0"

# 8) 提取 PASS 位点
"$gatk" --java-options "-Xmx12g -Djava.io.tmpdir=$TMP" SelectVariants \
    -R "$REF" \
    -V "$VCF/filtered/all.raw.snp.filtered.vcf.gz" \
    --exclude-filtered true \
    -O "$VCF/filtered/all.pass.snp.vcf.gz"

"$gatk" --java-options "-Xmx12g -Djava.io.tmpdir=$TMP" SelectVariants \
    -R "$REF" \
    -V "$VCF/filtered/all.raw.indel.filtered.vcf.gz" \
    --exclude-filtered true \
    -O "$VCF/filtered/all.pass.indel.vcf.gz"

# 9) 合并 PASS SNP 和 INDEL
"$gatk" --java-options "-Xmx12g -Djava.io.tmpdir=$TMP" MergeVcfs \
    -I "$VCF/filtered/all.pass.snp.vcf.gz" \
    -I "$VCF/filtered/all.pass.indel.vcf.gz" \
    -O "$VCF/filtered/all.pass.vcf.gz"

# 10) snpEff 注释
java -Xmx16g -jar "$snpeff" \
    -c "$snpeff_cfg" \
    -ud 5000 \
    -csvStats "$VCF/filtered/all.pass.csv" \
    -htmlStats "$VCF/filtered/all.pass.html" \
    -o vcf \
    "$SNPEFF_DB" \
    "$VCF/filtered/all.pass.vcf.gz" \
    > "$VCF/filtered/all.pass.ann.vcf"

echo "Done."