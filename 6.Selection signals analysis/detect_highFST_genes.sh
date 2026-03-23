#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 5 ]]; then
    echo "Usage: $0 <candidate_gene_list.txt> <all.vcf.gz> <guangdong_indica.list> <international_indica.list> <outdir> [gff=/project/data/ref/MSU7.gff]"
    exit 1
fi

CANDIDATE_GENES="$1"
VCF="$2"
GD_LIST="$3"
INT_LIST="$4"
OUTDIR="$5"
GFF="${6:-/project/data/ref/MSU7.gff}"

WINDOW_SIZE=100000
STEP_SIZE=10000
FST_COL=5   # 5 = WEIGHTED_FST in standard vcftools windowed FST output; change to 6 for MEAN_FST

mkdir -p "${OUTDIR}"/{fst,bed,genes,tmp}

# Check dependencies
command -v vcftools >/dev/null 2>&1 || { echo "ERROR: vcftools not found"; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo "ERROR: bedtools not found"; exit 1; }
command -v awk >/dev/null 2>&1 || { echo "ERROR: awk not found"; exit 1; }
command -v sort >/dev/null 2>&1 || { echo "ERROR: sort not found"; exit 1; }

echo "[INFO] Step 1: calculate windowed FST ..."
vcftools \
    --gzvcf "${VCF}" \
    --weir-fst-pop "${GD_LIST}" \
    --weir-fst-pop "${INT_LIST}" \
    --fst-window-size "${WINDOW_SIZE}" \
    --fst-window-step "${STEP_SIZE}" \
    --out "${OUTDIR}/fst/GD_vs_INT"

FST_FILE="${OUTDIR}/fst/GD_vs_INT.windowed.weir.fst"

if [[ ! -s "${FST_FILE}" ]]; then
    echo "ERROR: FST output file not found: ${FST_FILE}"
    exit 1
fi

echo "[INFO] Step 2: get top 10% FST threshold ..."

# Extract valid FST values, remove nan/NA, sort ascending
awk -v col="${FST_COL}" '
BEGIN{FS=OFS="\t"}
NR==1{next}
$col != "nan" && $col != "NA" && $col != "" {
    print $col
}
' "${FST_FILE}" | sort -g > "${OUTDIR}/tmp/fst_values.sorted.txt"

N=$(wc -l < "${OUTDIR}/tmp/fst_values.sorted.txt" | awk '{print $1}')
if [[ "${N}" -eq 0 ]]; then
    echo "ERROR: no valid FST values found."
    exit 1
fi

# 90th percentile cutoff = threshold for top 10%
RANK=$(awk -v n="${N}" 'BEGIN{
    r = int(0.9 * n)
    if (r < 1) r = 1
    print r
}')
THRESHOLD=$(awk -v rank="${RANK}" 'NR==rank{print; exit}' "${OUTDIR}/tmp/fst_values.sorted.txt")

echo "[INFO] Number of valid windows: ${N}"
echo "[INFO] 90th percentile FST threshold: ${THRESHOLD}"

echo "${THRESHOLD}" > "${OUTDIR}/fst/top10_threshold.txt"

echo "[INFO] Step 3: convert top 10% FST windows to BED ..."

awk -v col="${FST_COL}" -v thr="${THRESHOLD}" '
BEGIN{FS=OFS="\t"}
NR==1{next}
$col != "nan" && $col != "NA" && $col != "" && $col >= thr {
    start0 = $2 - 1
    if (start0 < 0) start0 = 0
    print $1, start0, $3, $col
}
' "${FST_FILE}" \
| sort -k1,1 -k2,2n -k3,3n \
> "${OUTDIR}/bed/top10_fst_windows.bed"

# Merge overlapping top 10% windows
bedtools merge \
    -i "${OUTDIR}/bed/top10_fst_windows.bed" \
> "${OUTDIR}/bed/top10_fst_windows.merged.bed"

echo "[INFO] Step 4: extract candidate gene coordinates from GFF ..."

# Build candidate gene BED from GFF
# Output BED columns: chr, start0, end1, gene_id, strand
awk -v OFS="\t" '
BEGIN{
    while((getline < "'"${CANDIDATE_GENES}"'") > 0){
        gsub(/\r/, "", $1)
        if($1!="") ids[$1]=1
    }
}
$3=="gene"{
    id=""
    name=""

    n=split($9, a, ";")
    for(i=1;i<=n;i++){
        if(a[i] ~ /^ID=/){
            id = substr(a[i], 4)
        } else if(a[i] ~ /^Name=/){
            name = substr(a[i], 6)
        } else if(a[i] ~ /^gene_id=/){
            if(id=="") id = substr(a[i], 9)
        }
    }

    sub(/^gene:/, "", id)
    sub(/^gene:/, "", name)

    hit=""
    if(id in ids) hit=id
    else if(name in ids) hit=name

    if(hit!=""){
        start0 = $4 - 1
        if(start0 < 0) start0 = 0
        print $1, start0, $5, hit, $7
    }
}
' "${GFF}" \
| sort -k1,1 -k2,2n -k3,3n -u \
> "${OUTDIR}/genes/candidate_genes.bed"

if [[ ! -s "${OUTDIR}/genes/candidate_genes.bed" ]]; then
    echo "ERROR: no candidate genes were matched in GFF."
    exit 1
fi

echo "[INFO] Step 5: intersect candidate genes with top 10% FST windows ..."

bedtools intersect \
    -a "${OUTDIR}/genes/candidate_genes.bed" \
    -b "${OUTDIR}/bed/top10_fst_windows.merged.bed" \
    -wa -u \
> "${OUTDIR}/genes/candidate_genes_in_top10FST.bed"

cut -f4 "${OUTDIR}/genes/candidate_genes_in_top10FST.bed" \
| sort -u \
> "${OUTDIR}/genes/candidate_genes_in_top10FST.txt"

echo "[INFO] Step 6: summary ..."

TOTAL_CAND=$(wc -l < "${OUTDIR}/genes/candidate_genes.bed" | awk '{print $1}')
PASS_CAND=$(wc -l < "${OUTDIR}/genes/candidate_genes_in_top10FST.txt" | awk '{print $1}')
TOP_WIN=$(wc -l < "${OUTDIR}/bed/top10_fst_windows.bed" | awk '{print $1}')
MERGED_TOP_WIN=$(wc -l < "${OUTDIR}/bed/top10_fst_windows.merged.bed" | awk '{print $1}')

{
    echo -e "Metric\tValue"
    echo -e "Total_candidate_genes\t${TOTAL_CAND}"
    echo -e "Candidate_genes_in_top10FST\t${PASS_CAND}"
    echo -e "Top10_windows_raw\t${TOP_WIN}"
    echo -e "Top10_windows_merged\t${MERGED_TOP_WIN}"
    echo -e "Top10_FST_threshold\t${THRESHOLD}"
} > "${OUTDIR}/summary.tsv"

echo "[INFO] Done."
echo "[INFO] Main outputs:"
echo "  ${OUTDIR}/fst/GD_vs_INT.windowed.weir.fst"
echo "  ${OUTDIR}/fst/top10_threshold.txt"
echo "  ${OUTDIR}/bed/top10_fst_windows.merged.bed"
echo "  ${OUTDIR}/genes/candidate_genes_in_top10FST.txt"
echo "  ${OUTDIR}/summary.tsv"