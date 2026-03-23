library(GenomicFeatures)
library(GenomicRanges)

# input
gtf_file <- "/project/data/ref/MSU7.gtf"
count_file <- "/project/data/rnaseq/count/gene_counts.txt"

# read featureCounts matrix
dat <- read.table(count_file, header = TRUE, sep = "\t", comment.char = "#", check.names = FALSE)
anno <- dat[, 1:6]
cnt <- dat[, 7:ncol(dat), drop = FALSE]

# simplify sample names
colnames(cnt) <- basename(colnames(cnt))
colnames(cnt) <- sub("\\.uniq\\.sorted\\.bam$", "", colnames(cnt))

# gene lengths from GTF (non-overlapping exonic length per gene)
txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")
exons_by_gene <- exonsBy(txdb, by = "gene")
gene_len <- sapply(exons_by_gene, function(x) sum(width(reduce(x))))

# keep common genes
common_genes <- intersect(anno$Geneid, names(gene_len))
cnt2 <- cnt[match(common_genes, anno$Geneid), , drop = FALSE]
rownames(cnt2) <- common_genes
gene_len2 <- gene_len[common_genes]

# TPM calculation
rpk <- sweep(cnt2, 1, gene_len2 / 1000, "/")
tpm <- sweep(rpk, 2, colSums(rpk), "/") * 1e6

# output
write.table(
  data.frame(Geneid = rownames(tpm), tpm, check.names = FALSE),
  file = "/project/data/rnaseq/count/gene_tpm.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)