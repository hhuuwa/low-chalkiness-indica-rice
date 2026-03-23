#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(geneHapR)
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop(
    paste(
      "Usage: Rscript step2b_genehapr_network_temporal_pedigree.R",
      "<expanded_vcf_dir> <sample_meta.tsv> <outdir> <genes>",
      "[group_col=Group] [ped_vcf_dir=NA]"
    )
  )
}

expanded_vcf_dir <- args[1]
meta_file        <- args[2]
outdir           <- args[3]
genes_arg        <- args[4]
group_col        <- ifelse(length(args) >= 5, args[5], "Group")
ped_vcf_dir      <- ifelse(length(args) >= 6, args[6], NA)

genes <- strsplit(genes_arg, ",")[[1]]
genes <- trimws(genes)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "per_gene"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "temporal"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "pedigree"), recursive = TRUE, showWarnings = FALSE)

meta <- fread(meta_file, data.table = FALSE)
meta$ID <- as.character(meta$ID)

#-----------------------------------
# helper functions
#-----------------------------------
is_true <- function(x) {
  tolower(as.character(x)) %in% c("1", "true", "t", "yes", "y")
}

split_accessions <- function(x) {
  if (length(x) == 0 || is.na(x) || x == "") return(character(0))
  y <- unlist(strsplit(as.character(x), "[,;[:space:]]+"))
  y[y != ""]
}

hap_to_df <- function(hap_obj) {
  as.data.frame(as.matrix(hap_obj), check.names = FALSE, stringsAsFactors = FALSE)
}

get_site_cols <- function(df) {
  lead <- df[[1]]
  chrom_row <- which(lead == "CHROM")
  if (length(chrom_row) != 1) stop("Cannot find CHROM row in hap object.")
  idx <- which(!is.na(df[chrom_row, ]) & df[chrom_row, ] != "")
  setdiff(idx, 1)
}

get_aux_cols <- function(df) {
  setdiff(seq_len(ncol(df)), c(1, get_site_cols(df)))
}

get_acc_col <- function(df) {
  aux <- get_aux_cols(df)
  hit <- aux[tolower(colnames(df)[aux]) == "accession"]
  if (length(hit) == 1) return(hit)
  if (length(aux) >= 1) return(aux[1])
  stop("Cannot find accession column in hapSummary.")
}

get_hap_rows <- function(df) {
  lead <- df[[1]]
  which(!lead %in% c("CHROM", "POS", "INFO", "ALLELE"))
}

hap_strings_from_summary <- function(hapSummary) {
  df <- hap_to_df(hapSummary)
  hap_rows <- get_hap_rows(df)
  site_cols <- get_site_cols(df)

  data.frame(
    Hap = df[hap_rows, 1],
    HapString = apply(df[hap_rows, site_cols, drop = FALSE], 1, paste, collapse = "|"),
    stringsAsFactors = FALSE
  )
}

sample_hap_from_summary <- function(hapSummary) {
  df <- hap_to_df(hapSummary)
  hap_rows <- get_hap_rows(df)
  acc_col <- get_acc_col(df)

  res <- vector("list", length(hap_rows))
  for (i in seq_along(hap_rows)) {
    r <- hap_rows[i]
    hap_name <- df[r, 1]
    accs <- split_accessions(df[r, acc_col])

    if (length(accs) == 0) {
      res[[i]] <- NULL
    } else {
      res[[i]] <- data.frame(
        ID = accs,
        Hap = hap_name,
        stringsAsFactors = FALSE
      )
    }
  }
  bind_rows(res)
}

make_hapNet_safe <- function(hapSummary, AccINFO, group_col = "Group") {
  exports <- getNamespaceExports("geneHapR")
  if ("get_hapNet" %in% exports) {
    get_hapNet(hapSummary, AccINFO = AccINFO, groupName = group_col)
  } else if ("network" %in% exports) {
    network(hapSummary, AccINFO = AccINFO, groupName = group_col)
  } else {
    stop("Neither get_hapNet nor network found in geneHapR.")
  }
}

map_query_haps_to_reference <- function(ref_hapSummary, query_hapSummary) {
  ref_def <- hap_strings_from_summary(ref_hapSummary)
  qry_def <- hap_strings_from_summary(query_hapSummary)
  qry_sample_hap <- sample_hap_from_summary(query_hapSummary)

  qry_map <- qry_def %>%
    left_join(ref_def, by = "HapString", suffix = c(".query", ".ref")) %>%
    mutate(RefHap = ifelse(is.na(Hap.ref), "Novel", Hap.ref)) %>%
    select(QueryHap = Hap.query, HapString, RefHap)

  qry_sample_hap %>%
    left_join(qry_map, by = c("Hap" = "QueryHap"))
}

#-----------------------------------
# main analysis for expanded panel
#-----------------------------------
AccINFO <- import_AccINFO(meta_file)
ref_hapSummary_list <- list()
all_sample_hap <- list()

for (gene in genes) {
  message("[INFO] Processing expanded panel: ", gene)

  vcf_file <- file.path(expanded_vcf_dir, paste0(gene, ".vcf.gz"))
  if (!file.exists(vcf_file)) {
    stop(paste("Missing VCF:", vcf_file))
  }

  gene_dir <- file.path(outdir, "per_gene", gene)
  dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE)

  vcf <- import_vcf(vcf_file)

  hapResult <- vcf2hap(
    vcf,
    hapPrefix = "H",
    hetero_remove = TRUE,
    na_drop = TRUE
  )

  hapSummary <- hap_summary(hapResult)

  write.table(
    as.matrix(hapResult),
    file = file.path(gene_dir, paste0(gene, ".expanded.hapResult.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
  )

  write.table(
    as.matrix(hapSummary),
    file = file.path(gene_dir, paste0(gene, ".expanded.hapSummary.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
  )

  sample_hap <- sample_hap_from_summary(hapSummary) %>%
    mutate(Gene = gene) %>%
    left_join(meta, by = "ID")

  write.table(
    sample_hap,
    file = file.path(gene_dir, paste0(gene, ".expanded.sample_hap.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  hapNet <- make_hapNet_safe(hapSummary, AccINFO, group_col = group_col)

  pdf(file.path(gene_dir, paste0(gene, ".expanded.hapNet.pdf")), width = 8, height = 6)
  plotHapNet(
    hapNet,
    size = "freq",
    scale = "log2",
    cex = 0.9,
    col.link = 2,
    link.width = 2,
    show.mutation = 2
  )
  title(main = gene)
  dev.off()

  pdf(file.path(gene_dir, paste0(gene, ".expanded.hapTable.pdf")), width = 10, height = 4.5)
  plotHapTable(hapSummary)
  dev.off()

  ref_hapSummary_list[[gene]] <- hapSummary
  all_sample_hap[[gene]] <- sample_hap
}

all_sample_hap_df <- bind_rows(all_sample_hap)

write.table(
  all_sample_hap_df,
  file = file.path(outdir, "all_genes.expanded.sample_hap.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

#-----------------------------------
# temporal haplotype frequency
#-----------------------------------
if (all(c("IsRepresentativeChinese", "ReleaseYear") %in% colnames(meta))) {

  temporal_df <- all_sample_hap_df %>%
    filter(is_true(IsRepresentativeChinese), !is.na(ReleaseYear), ReleaseYear != "NA") %>%
    mutate(
      ReleaseYear = as.numeric(ReleaseYear),
      Decade = paste0(floor(ReleaseYear / 10) * 10, "s")
    ) %>%
    count(Gene, Decade, Hap, name = "Count") %>%
    group_by(Gene, Decade) %>%
    mutate(Frequency = Count / sum(Count)) %>%
    ungroup()

  write.table(
    temporal_df,
    file = file.path(outdir, "temporal", "representative_chinese.haplotype_frequency.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  for (gene in genes) {
    pdat <- temporal_df %>% filter(Gene == gene)
    if (nrow(pdat) == 0) next

    p <- ggplot(pdat, aes(x = Decade, y = Frequency, fill = Hap)) +
      geom_col() +
      theme_bw(base_size = 12) +
      labs(
        title = paste0(gene, " haplotype frequency across decades"),
        x = "Release decade",
        y = "Frequency"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()
      )

    ggsave(
      filename = file.path(outdir, "temporal", paste0(gene, ".representative_chinese.temporal_frequency.pdf")),
      plot = p, width = 7, height = 4.5
    )
  }
}

#-----------------------------------
# Huanghuazhan pedigree
#-----------------------------------
# Case 1: pedigree samples already included in expanded panel
if ("IsHHZPedigree" %in% colnames(meta)) {
  hhz_from_main <- all_sample_hap_df %>%
    filter(is_true(IsHHZPedigree)) %>%
    select(ID, Gene, Hap)

  if (nrow(hhz_from_main) > 0) {
    write.table(
      hhz_from_main,
      file = file.path(outdir, "pedigree", "HHZ_pedigree.hap_assignment.from_expanded_panel.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )

    hhz_matrix <- hhz_from_main %>%
      tidyr::pivot_wider(names_from = Gene, values_from = Hap)

    write.table(
      hhz_matrix,
      file = file.path(outdir, "pedigree", "HHZ_pedigree.hap_matrix.from_expanded_panel.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
  }
}

# Case 2: pedigree is provided as separate per-gene VCFs
if (!is.na(ped_vcf_dir)) {
  ped_all <- list()

  for (gene in genes) {
    ped_vcf <- file.path(ped_vcf_dir, paste0(gene, ".vcf.gz"))
    if (!file.exists(ped_vcf)) next

    message("[INFO] Processing HHZ pedigree: ", gene)

    ped_vcf_obj <- import_vcf(ped_vcf)
    ped_hapResult <- vcf2hap(
      ped_vcf_obj,
      hapPrefix = "H",
      hetero_remove = TRUE,
      na_drop = TRUE
    )
    ped_hapSummary <- hap_summary(ped_hapResult)

    mapped <- map_query_haps_to_reference(ref_hapSummary_list[[gene]], ped_hapSummary) %>%
      mutate(Gene = gene)

    write.table(
      mapped,
      file = file.path(outdir, "pedigree", paste0(gene, ".HHZ_pedigree.refHap.tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )

    ped_all[[gene]] <- mapped
  }

  if (length(ped_all) > 0) {
    ped_all_df <- bind_rows(ped_all)

    write.table(
      ped_all_df,
      file = file.path(outdir, "pedigree", "HHZ_pedigree.hap_assignment.from_separate_vcf.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
  }
}

message("[INFO] All analyses completed.")