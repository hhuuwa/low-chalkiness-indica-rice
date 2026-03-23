#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 4 ]]; then
    echo "Usage: $0 <gene_id.txt> <gff> <all.vcf.gz> <outdir> [promoter_bp]"
    exit 1
fi

GENE_LIST="$1"
GFF="$2"
VCF="$3"
OUTDIR="$4"
PROMOTER_BP="${5:-2000}"

BED_DIR="${OUTDIR}/bed"
VCF_DIR="${OUTDIR}/vcf"
TMP_DIR="${OUTDIR}/tmp"

mkdir -p "${BED_DIR}" "${VCF_DIR}" "${TMP_DIR}"

# Check required programs
command -v awk >/dev/null 2>&1 || { echo "ERROR: awk not found"; exit 1; }
command -v bcftools >/dev/null 2>&1 || { echo "ERROR: bcftools not found"; exit 1; }

# Index the input VCF if no index is found
if [[ ! -f "${VCF}.tbi" && ! -f "${VCF}.csi" ]]; then
    echo "[INFO] Indexing VCF..."
    bcftools index -t "${VCF}"
fi

echo "[INFO] Step 1: build promoter + CDS BED files ..."

awk -v OFS="\t" -v outdir="${BED_DIR}" -v promoter_bp="${PROMOTER_BP}" '
function get_attr(attr, key,   n, i, a, kv) {
    n = split(attr, a, ";")
    for (i = 1; i <= n; i++) {
        split(a[i], kv, "=")
        if (kv[1] == key) return kv[2]
    }
    return ""
}
BEGIN {
    # Load target gene IDs
    while ((getline < "'"${GENE_LIST}"'") > 0) {
        gsub(/\r/, "", $0)
        if ($1 != "") target[$1] = 1
    }
}
# First pass: store coordinates/strand for target genes
# and build transcript-to-gene mapping
FNR == NR {
    if ($0 ~ /^#/ || NF < 9) next

    feature = $3
    attr = $9

    if (feature == "gene") {
        id   = get_attr(attr, "ID")
        name = get_attr(attr, "Name")

        hit = ""
        if (id in target) hit = id
        else if (name in target) hit = name

        if (hit != "") {
            gene_chr[hit]    = $1
            gene_start[hit]  = $4
            gene_end[hit]    = $5
            gene_strand[hit] = $7
            found[hit] = 1
        }
    }
    else if (feature == "mRNA" || feature == "transcript") {
        tx_id  = get_attr(attr, "ID")
        parent = get_attr(attr, "Parent")

        if (parent in found) {
            tx2gene[tx_id] = parent
        }
    }
    next
}
# Second pass: write promoter and CDS intervals into BED files
{
    if ($0 ~ /^#/ || NF < 9) next

    feature = $3
    attr = $9

    # Output promoter region once for each target gene
    if (feature == "gene") {
        id   = get_attr(attr, "ID")
        name = get_attr(attr, "Name")

        gene = ""
        if (id in found) gene = id
        else if (name in found) gene = name

        if (gene != "" && !(gene in promoter_done)) {
            chr    = gene_chr[gene]
            start  = gene_start[gene]
            end    = gene_end[gene]
            strand = gene_strand[gene]

            # Strand-specific upstream promoter region
            # BED format is 0-based, half-open
            if (strand == "+") {
                p_start0 = start - promoter_bp - 1
                if (p_start0 < 0) p_start0 = 0
                p_end = start - 1
            } else {
                p_start0 = end
                p_end = end + promoter_bp
            }

            bedfile = outdir "/" gene ".bed"
            print chr, p_start0, p_end, gene, "PROMOTER", strand >> bedfile
            promoter_done[gene] = 1
        }
    }

    # Output CDS intervals
    if (feature == "CDS") {
        parent = get_attr(attr, "Parent")
        n = split(parent, pp, ",")

        for (i = 1; i <= n; i++) {
            gene = ""
            if (pp[i] in tx2gene) gene = tx2gene[pp[i]]
            else if (pp[i] in found) gene = pp[i]   # compatible with CDS directly linked to gene

            if (gene != "") {
                chr    = $1
                start0 = $4 - 1
                end1   = $5
                strand = $7
                bedfile = outdir "/" gene ".bed"
                print chr, start0, end1, gene, "CDS", strand >> bedfile
            }
        }
    }
}
' "${GFF}" "${GFF}"

# Sort and deduplicate each BED file
for bed in "${BED_DIR}"/*.bed; do
    [[ -e "$bed" ]] || continue
    sort -k1,1 -k2,2n -k3,3n -u "$bed" > "${TMP_DIR}/$(basename "$bed")"
    mv "${TMP_DIR}/$(basename "$bed")" "$bed"
done

# Generate a summary table for gene and promoter coordinates
echo -e "GeneID\tChr\tGeneStart\tGeneEnd\tStrand\tPromoterRegion" > "${OUTDIR}/gene_region_summary.tsv"

awk -v OFS="\t" -v promoter_bp="${PROMOTER_BP}" '
function get_attr(attr, key,   n, i, a, kv) {
    n = split(attr, a, ";")
    for (i = 1; i <= n; i++) {
        split(a[i], kv, "=")
        if (kv[1] == key) return kv[2]
    }
    return ""
}
BEGIN {
    while ((getline < "'"${GENE_LIST}"'") > 0) {
        gsub(/\r/, "", $0)
        if ($1 != "") target[$1] = 1
    }
}
$0 !~ /^#/ && NF >= 9 && $3=="gene" {
    id   = get_attr($9, "ID")
    name = get_attr($9, "Name")

    hit = ""
    if (id in target) hit = id
    else if (name in target) hit = name

    if (hit != "") {
        if ($7 == "+") {
            ps = $4 - promoter_bp
            if (ps < 1) ps = 1
            pe = $4 - 1
        } else {
            ps = $5 + 1
            pe = $5 + promoter_bp
        }
        print hit, $1, $4, $5, $7, $1 ":" ps "-" pe
    }
}
' "${GFF}" >> "${OUTDIR}/gene_region_summary.tsv"

echo "[INFO] Step 2: extract VCF for each gene ..."

while read -r gene; do
    [[ -z "${gene}" ]] && continue

    bed="${BED_DIR}/${gene}.bed"
    outvcf="${VCF_DIR}/${gene}.vcf.gz"

    if [[ ! -s "${bed}" ]]; then
        echo "[WARN] ${gene}: no BED generated, skip."
        continue
    fi

    echo "[INFO] Extracting ${gene} ..."
    bcftools view \
        -T "${bed}" \
        -Oz \
        -o "${outvcf}" \
        "${VCF}"

    nvar=$(bcftools view -H "${outvcf}" | wc -l | awk "{print \$1}")

    if [[ "${nvar}" -eq 0 ]]; then
        echo "[WARN] ${gene}: no variants found, remove empty VCF."
        rm -f "${outvcf}" "${outvcf}.tbi" "${outvcf}.csi"
        continue
    fi

    bcftools index -t "${outvcf}"
done < "${GENE_LIST}"

# Merge all gene BED files into one combined BED file
cat "${BED_DIR}"/*.bed 2>/dev/null | sort -k1,1 -k2,2n -k3,3n -u > "${OUTDIR}/all.promoter_cds.bed" || true

echo "[INFO] Done."
echo "[INFO] BED dir : ${BED_DIR}"
echo "[INFO] VCF dir : ${VCF_DIR}"
echo "[INFO] Summary : ${OUTDIR}/gene_region_summary.tsv"