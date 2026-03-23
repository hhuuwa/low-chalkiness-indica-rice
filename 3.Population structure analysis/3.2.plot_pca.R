library(ggplot2)

prefix <- "/project/data/popgen/rice.maf05.miss95.pca"

eigvec <- read.table(paste0(prefix, ".eigenvec"), header = FALSE, stringsAsFactors = FALSE)
eigval <- scan(paste0(prefix, ".eigenval"))

colnames(eigvec) <- c("FID", "IID", paste0("PC", 1:(ncol(eigvec) - 2)))
var_exp <- eigval / sum(eigval) * 100

p <- ggplot(eigvec, aes(x = PC1, y = PC2)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(
    x = paste0("PC1 (", round(var_exp[1], 2), "%)"),
    y = paste0("PC2 (", round(var_exp[2], 2), "%)")
  )

ggsave("/project/data/popgen/rice.maf05.miss95.pca.PC1_PC2.pdf", p, width = 6, height = 5)

write.table(
  data.frame(PC = paste0("PC", 1:length(var_exp)), VarianceExplained = var_exp),
  "/project/data/popgen/rice.maf05.miss95.pca.variance.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)