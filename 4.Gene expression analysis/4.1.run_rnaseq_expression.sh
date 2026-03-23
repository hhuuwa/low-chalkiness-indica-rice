#!/bin/bash
set -euo pipefail

# software
fastp="/project/data/tool/fastp"
star="/project/data/tool/STAR"
samtools="/project/data/tool/samtools"
featureCounts="/project/data/tool/featureCounts"

# reference
REF_FA="/project/data/ref/MSU7.fa"
GTF="/project/data/ref/MSU7.gtf"
STAR_INDEX="/project/data/ref/MSU7_STAR_index"

# input/output
SAMPLE_LIST="/project/data/rnaseq/sample_list.tsv"
OUTDIR="/project/data/rnaseq"
THREADS=16

CLEAN_DIR="${OUTDIR}/clean"
ALIGN_DIR="${OUTDIR}/align"
COUNT_DIR="${OUTDIR}/count"
LOG_DIR="${OUTDIR}/log"

mkdir -p "$CLEAN_DIR" "$ALIGN_DIR" "$COUNT_DIR" "$LOG_DIR"

# 1. build STAR index (run once)
if [ ! -d "$STAR_INDEX" ]; then
    mkdir -p "$STAR_INDEX"
    "$star" \
        --runThreadN "$THREADS" \
        --runMode genomeGenerate \
        --genomeDir "$STAR_INDEX" \
        --genomeFastaFiles "$REF_FA" \
        --sjdbGTFfile "$GTF" \
        --sjdbOverhang 149
fi

# 2. fastp + STAR alignment
tail -n +2 "$SAMPLE_LIST" | while IFS=$'\t' read -r sample r1 r2; do
    echo "Processing $sample ..."

    clean_r1="${CLEAN_DIR}/${sample}_clean_R1.fq.gz"
    clean_r2="${CLEAN_DIR}/${sample}_clean_R2.fq.gz"

    # fastp
    "$fastp" \
        -i "$r1" \
        -I "$r2" \
        -o "$clean_r1" \
        -O "$clean_r2" \
        -h "${LOG_DIR}/${sample}.fastp.html" \
        -j "${LOG_DIR}/${sample}.fastp.json" \
        -w "$THREADS"

    # STAR unique mapping only
    "$star" \
        --runThreadN "$THREADS" \
        --genomeDir "$STAR_INDEX" \
        --readFilesIn "$clean_r1" "$clean_r2" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${ALIGN_DIR}/${sample}." \
        --outSAMtype BAM SortedByCoordinate \
        --outFilterMultimapNmax 1 \
        --outSAMattributes NH HI AS nM

    mv "${ALIGN_DIR}/${sample}.Aligned.sortedByCoord.out.bam" \
       "${ALIGN_DIR}/${sample}.uniq.sorted.bam"

    "$samtools" index "${ALIGN_DIR}/${sample}.uniq.sorted.bam"
done

# 3. featureCounts on uniquely mapped paired reads
"$featureCounts" \
    -T "$THREADS" \
    -p \
    -B \
    -C \
    -a "$GTF" \
    -o "${COUNT_DIR}/gene_counts.txt" \
    "${ALIGN_DIR}"/*.uniq.sorted.bam