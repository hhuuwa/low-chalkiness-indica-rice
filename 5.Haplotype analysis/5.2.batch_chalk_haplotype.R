#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(geneHapR)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggsignif)
})

theme_set(
  theme_bw(base_size = 13) +
    theme(
      text = element_text(face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop(
    paste(
      "Usage: Rscript batch_chalk_haplotype.R <vcf_dir> <accinfo.tsv> <pheno.tsv> <outdir>",
      "[trait_col=Chalkiness] [acc_id_col=ID] [subpop_col=Subpopulation]",
      "[gd_group_name=Guangdong_indica] [hap_prefix=H]",
      "[hap_freq_ratio_cutoff=0.05] [gd_freq_cutoff=0.5]",
      "[fdr_cutoff=0.05] [p_cutoff=0.05]"
    )
  )
}

vcf_dir <- args[1]
accinfo_file <- args[2]
pheno_file <- args[3]
outdir <- args[4]

trait_col <- ifelse(length(args) >= 5, args[5], "Chalkiness")
acc_id_col <- ifelse(length(args) >= 6, args[6], "ID")
subpop_col <- ifelse(length(args) >= 7, args[7], "Subpopulation")
gd_group_name <- ifelse(length(args) >= 8, args[8], "Guangdong_indica")
hap_prefix <- ifelse(length(args) >= 9, args[9], "H")

hap_freq_ratio_cutoff <- ifelse(length(args) >= 10, as.numeric(args[10]), 0.05)
gd_freq_cutoff <- ifelse(length(args) >= 11, as.numeric(args[11]), 0.50)
fdr_cutoff <- ifelse(length(args) >= 12, as.numeric(args[12]), 0.05)
p_cutoff <- ifelse(length(args) >= 13, as.numeric(args[13]), 0.05)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "per_gene"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "summary"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "logs"), recursive = TRUE, showWarnings = FALSE)

# =========================
# Helper functions
# =========================

strip_meta_rows <- function(x) {
  df <- as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(df) <= 4) return(NULL)
  df[-(1:4), , drop = FALSE]
}

safe_gene_id <- function(vcf_file) {
  x <- basename(vcf_file)
  x <- sub("\\.vcf\\.gz$", "", x, ignore.case = TRUE)
  x <- sub("\\.vcf$", "", x, ignore.case = TRUE)
  x
}

make_pairwise_tests <- function(df, trait_name, p_cutoff = 0.05) {
  # df must contain: Hap, Value
  hap_stat <- df %>%
    group_by(Hap) %>%
    summarise(
      n = n(),
      mean_trait = mean(Value, na.rm = TRUE),
      sd_trait = sd(Value, na.rm = TRUE),
      .groups = "drop"
    )

  haps <- hap_stat$Hap
  pair_res <- list()
  idx <- 1

  if (length(haps) >= 2) {
    for (i in 1:(length(haps) - 1)) {
      for (j in (i + 1):length(haps)) {
        h1 <- haps[i]
        h2 <- haps[j]

        x1 <- df$Value[df$Hap == h1]
        x2 <- df$Value[df$Hap == h2]

        if (length(x1) < 2 || length(x2) < 2) next

        tt <- try(t.test(x1, x2), silent = TRUE)
        if (inherits(tt, "try-error")) next

        pair_res[[idx]] <- data.frame(
          Hap1 = h1,
          Hap2 = h2,
          n1 = length(x1),
          n2 = length(x2),
          mean1 = mean(x1, na.rm = TRUE),
          mean2 = mean(x2, na.rm = TRUE),
          p_value = tt$p.value,
          stringsAsFactors = FALSE
        )
        idx <- idx + 1
      }
    }
  }

  if (length(pair_res) == 0) {
    pair_df <- data.frame(
      Hap1 = character(),
      Hap2 = character(),
      n1 = integer(),
      n2 = integer(),
      mean1 = numeric(),
      mean2 = numeric(),
      p_value = numeric(),
      stringsAsFactors = FALSE
    )
  } else {
    pair_df <- bind_rows(pair_res)
  }

  is_superior_hap <- function(h, pair_df, n_expected, p_cutoff = 0.05) {
    cmp <- pair_df %>% filter(Hap1 == h | Hap2 == h)
    if (nrow(cmp) == 0) return(FALSE)
    if (nrow(cmp) != n_expected) return(FALSE)

    all(apply(cmp, 1, function(z) {
      z <- as.list(z)
      if (z$Hap1 == h) {
        as.numeric(z$mean1) < as.numeric(z$mean2) &&
          as.numeric(z$p_value) < p_cutoff
      } else {
        as.numeric(z$mean2) < as.numeric(z$mean1) &&
          as.numeric(z$p_value) < p_cutoff
      }
    }))
  }

  n_expected <- max(nrow(hap_stat) - 1, 0)

  superior_df <- hap_stat %>%
    rowwise() %>%
    mutate(
      superior_haplotype = is_superior_hap(Hap, pair_df, n_expected, p_cutoff)
    ) %>%
    ungroup()

  list(pair_df = pair_df, superior_df = superior_df)
}

plot_hap_composition_within_hap <- function(df, subpop_col, geneID, outfile) {
  tmp <- df %>%
    filter(!is.na(Hap), Hap != "", !is.na(.data[[subpop_col]])) %>%
    group_by(Hap, Group = .data[[subpop_col]]) %>%
    summarise(Freq = n(), .groups = "drop") %>%
    group_by(Hap) %>%
    mutate(Sum_Freq = sum(Freq),
           Percent = Freq / Sum_Freq,
           Label = paste0(round(100 * Percent, 1), "% (", Freq, ")")) %>%
    ungroup()

  if (nrow(tmp) == 0) return(NULL)

  hap_order <- df %>%
    filter(!is.na(Hap), Hap != "") %>%
    count(Hap, sort = TRUE) %>%
    pull(Hap)

  tmp$Hap <- factor(tmp$Hap, levels = hap_order)

  xlab_df <- tmp %>%
    distinct(Hap, Sum_Freq) %>%
    arrange(Hap) %>%
    mutate(labelX = paste0(Hap, " (", Sum_Freq, ")"))

  p <- ggplot(tmp, aes(x = Hap, y = Percent, fill = Group)) +
    geom_col(width = 0.8) +
    geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 3) +
    scale_x_discrete(labels = xlab_df$labelX) +
    labs(title = geneID, x = "", y = "Percent within haplotype") +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  ggsave(outfile, p, width = 7, height = 5)
}

plot_hap_composition_within_group <- function(df, subpop_col, geneID, outfile) {
  tmp <- df %>%
    filter(!is.na(Hap), Hap != "", !is.na(.data[[subpop_col]])) %>%
    group_by(Group = .data[[subpop_col]], Hap) %>%
    summarise(Freq = n(), .groups = "drop") %>%
    group_by(Group) %>%
    mutate(Sum_Freq = sum(Freq),
           Percent = Freq / Sum_Freq,
           Label = paste0(round(100 * Percent, 1), "% (", Freq, ")")) %>%
    ungroup()

  if (nrow(tmp) == 0) return(NULL)

  group_order <- tmp %>%
    group_by(Group) %>%
    summarise(Total = sum(Freq), .groups = "drop") %>%
    arrange(desc(Total)) %>%
    pull(Group)

  tmp$Group <- factor(tmp$Group, levels = group_order)

  xlab_df <- tmp %>%
    distinct(Group, Sum_Freq) %>%
    arrange(Group) %>%
    mutate(labelX = paste0(Group, " (", Sum_Freq, ")"))

  p <- ggplot(tmp, aes(x = Group, y = Percent, fill = Hap)) +
    geom_col(width = 0.8) +
    geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 3) +
    scale_x_discrete(labels = xlab_df$labelX) +
    labs(title = geneID, x = "", y = "Percent within subpopulation") +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  ggsave(outfile, p, width = 7, height = 5)
}

plot_trait_boxplot <- function(df, trait_col, subpop_col, geneID, outfile) {
  df_trait <- df %>%
    filter(!is.na(Hap), Hap != "", !is.na(.data[[trait_col]])) %>%
    transmute(
      Hap = as.character(Hap),
      Trait = trait_col,
      Value = as.numeric(.data[[trait_col]]),
      Group = as.character(.data[[subpop_col]])
    )

  if (nrow(df_trait) == 0) return(NULL)

  df_all <- df_trait %>%
    mutate(Group = "All")

  df_plot <- bind_rows(df_trait, df_all)

  hap_levels <- df_plot %>%
    group_by(Hap) %>%
    summarise(mean_value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
    arrange(mean_value) %>%
    pull(Hap)

  df_plot$Hap <- factor(df_plot$Hap, levels = hap_levels)

  pairlist <- NULL
  if (length(hap_levels) >= 2 && length(hap_levels) <= 8) {
    pairlist <- combn(hap_levels, 2, simplify = FALSE)
  }

  p <- ggplot(df_plot, aes(x = Hap, y = Value, fill = Hap, color = Hap)) +
    geom_boxplot(alpha = 0.7, width = 0.45, outlier.shape = NA) +
    facet_grid(Trait ~ Group, scales = "free_y") +
    labs(title = geneID, x = "", y = trait_col) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.18)))

  if (!is.null(pairlist)) {
    p <- p + geom_signif(
      comparisons = pairlist,
      step_increase = 0.12,
      map_signif_level = FALSE,
      test = "t.test",
      textsize = 3,
      size = 0.25,
      color = "black"
    )
  }

  ggsave(outfile, p, width = 8, height = 4.5)
}

analyze_one_gene <- function(vcf_file, accinfo, pheno) {
  geneID <- safe_gene_id(vcf_file)
  gene_dir <- file.path(outdir, "per_gene", geneID)
  dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE)

  message("[INFO] Processing ", geneID)

  vcf <- import_vcf(vcf_file)
  hapResult <- vcf2hap(
    vcf,
    hapPrefix = hap_prefix,
    hetero_remove = TRUE,
    na_drop = FALSE
  )

  hap_df_raw <- strip_meta_rows(hapResult)
  if (is.null(hap_df_raw) || !all(c("Accession", "Hap") %in% colnames(hap_df_raw))) {
    stop("No usable hapResult rows found.")
  }

  hap_df_raw <- hap_df_raw %>%
    transmute(
      Accession = as.character(Accession),
      Hap = as.character(Hap)
    ) %>%
    filter(!is.na(Accession), Accession != "", !is.na(Hap), Hap != "") %>%
    distinct()

  if (nrow(hap_df_raw) == 0) {
    stop("No valid accessions remained after haplotype inference.")
  }

  cutfreq <- max(1, round(nrow(hap_df_raw) * hap_freq_ratio_cutoff))

  phapResult <- filter_hap(
    hapResult,
    rm.mode = c("freq"),
    freq.min = cutfreq
  )

  phap_df <- strip_meta_rows(phapResult)
  if (is.null(phap_df) || !all(c("Accession", "Hap") %in% colnames(phap_df))) {
    stop("No usable rows remained after haplotype frequency filtering.")
  }

  phap_df <- phap_df %>%
    transmute(
      Accession = as.character(Accession),
      Hap = as.character(Hap)
    ) %>%
    filter(!is.na(Accession), Accession != "", !is.na(Hap), Hap != "") %>%
    distinct()

  if (nrow(phap_df) == 0) {
    stop("No valid accessions remained after haplotype filtering.")
  }

  phapSummary <- hap_summary(phapResult, hapPrefix = hap_prefix)

  write.hap(phapResult, file = file.path(gene_dir, paste0(geneID, ".hapResult.tsv")), sep = "\t")
  write.hap(phapSummary, file = file.path(gene_dir, paste0(geneID, ".hapSummary.tsv")), sep = "\t")

  pdf(file.path(gene_dir, paste0(geneID, ".hapTable.pdf")), width = 8, height = 6)
  plotHapTable(
    phapSummary,
    hapPrefix = hap_prefix,
    angle = 45,
    displayIndelSize = 0,
    title = geneID
  )
  dev.off()

  merged <- phap_df %>%
    left_join(accinfo, by = c("Accession" = acc_id_col)) %>%
    left_join(pheno, by = c("Accession" = acc_id_col))

  fwrite(
    merged,
    file = file.path(gene_dir, paste0(geneID, ".allele.summary.tsv")),
    sep = "\t"
  )

  # Raw enrichment p-values; BH correction will be applied globally after all genes finish
  N_total <- merged %>%
    filter(!is.na(Hap), Hap != "") %>%
    nrow()

  K_gd_total <- merged %>%
    filter(!is.na(Hap), Hap != "", .data[[subpop_col]] == gd_group_name) %>%
    nrow()

  enrich_df <- merged %>%
    filter(!is.na(Hap), Hap != "") %>%
    mutate(is_gd = .data[[subpop_col]] == gd_group_name) %>%
    group_by(Hap) %>%
    summarise(
      n_hap = n(),
      k_gd = sum(is_gd, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Gene = geneID,
      N_total = N_total,
      K_gd_total = K_gd_total,
      freq_total = n_hap / N_total,
      freq_gd = ifelse(K_gd_total > 0, k_gd / K_gd_total, NA_real_),
      p_hyper = ifelse(
        K_gd_total > 0,
        phyper(k_gd - 1, K_gd_total, N_total - K_gd_total, n_hap, lower.tail = FALSE),
        NA_real_
      )
    ) %>%
    select(Gene, Hap, n_hap, k_gd, N_total, K_gd_total, freq_total, freq_gd, p_hyper)

  fwrite(
    enrich_df,
    file = file.path(gene_dir, paste0(geneID, ".hap.enrichment.raw.tsv")),
    sep = "\t"
  )

  # Trait comparison and superior haplotype
  trait_df <- merged %>%
    filter(!is.na(Hap), Hap != "", !is.na(.data[[trait_col]])) %>%
    transmute(
      Hap = as.character(Hap),
      Accession = as.character(Accession),
      Value = as.numeric(.data[[trait_col]])
    )

  if (nrow(trait_df) > 0) {
    tst <- make_pairwise_tests(trait_df, trait_col, p_cutoff = p_cutoff)
    pair_df <- tst$pair_df %>% mutate(Gene = geneID, .before = 1)
    superior_df <- tst$superior_df %>% mutate(Gene = geneID, .before = 1)
  } else {
    pair_df <- data.frame(
      Gene = character(),
      Hap1 = character(),
      Hap2 = character(),
      n1 = integer(),
      n2 = integer(),
      mean1 = numeric(),
      mean2 = numeric(),
      p_value = numeric(),
      stringsAsFactors = FALSE
    )
    superior_df <- data.frame(
      Gene = geneID,
      Hap = character(),
      n = integer(),
      mean_trait = numeric(),
      sd_trait = numeric(),
      superior_haplotype = logical(),
      stringsAsFactors = FALSE
    )
  }

  fwrite(
    pair_df,
    file = file.path(gene_dir, paste0(geneID, ".pairwise_ttest.tsv")),
    sep = "\t"
  )

  fwrite(
    superior_df,
    file = file.path(gene_dir, paste0(geneID, ".superior_haplotype.tsv")),
    sep = "\t"
  )

  plot_hap_composition_within_hap(
    merged,
    subpop_col = subpop_col,
    geneID = geneID,
    outfile = file.path(gene_dir, paste0(geneID, ".hap.class1.pdf"))
  )

  plot_hap_composition_within_group(
    merged,
    subpop_col = subpop_col,
    geneID = geneID,
    outfile = file.path(gene_dir, paste0(geneID, ".hap.class2.pdf"))
  )

  plot_trait_boxplot(
    merged,
    trait_col = trait_col,
    subpop_col = subpop_col,
    geneID = geneID,
    outfile = file.path(gene_dir, paste0(geneID, ".hap.trait.boxplot.pdf"))
  )

  list(
    gene = geneID,
    merged = merged,
    enrich = enrich_df,
    pair = pair_df,
    superior = superior_df,
    status = data.frame(
      Gene = geneID,
      Status = "OK",
      Accessions_after_filter = nrow(phap_df),
      Haplotypes = dplyr::n_distinct(phap_df$Hap),
      stringsAsFactors = FALSE
    )
  )
}

# =========================
# Input tables
# =========================
accinfo <- fread(accinfo_file, data.table = FALSE)
pheno <- fread(pheno_file, data.table = FALSE)

if (!(acc_id_col %in% colnames(accinfo))) {
  stop(paste("Column", acc_id_col, "not found in accinfo file."))
}
if (!(subpop_col %in% colnames(accinfo))) {
  stop(paste("Column", subpop_col, "not found in accinfo file."))
}
if (!(acc_id_col %in% colnames(pheno))) {
  names(pheno)[1] <- acc_id_col
}
if (!(trait_col %in% colnames(pheno))) {
  stop(paste("Column", trait_col, "not found in phenotype file."))
}

accinfo[[acc_id_col]] <- as.character(accinfo[[acc_id_col]])
pheno[[acc_id_col]] <- as.character(pheno[[acc_id_col]])
pheno[[trait_col]] <- suppressWarnings(as.numeric(pheno[[trait_col]]))

# =========================
# List VCF files
# =========================
vcf_files <- list.files(
  vcf_dir,
  pattern = "\\.vcf(\\.gz)?$",
  full.names = TRUE
)

if (length(vcf_files) == 0) {
  stop("No VCF files found in vcf_dir.")
}

# =========================
# Batch run
# =========================
all_enrich <- list()
all_pair <- list()
all_superior <- list()
all_status <- list()
failed <- list()

ie <- 1
ip <- 1
isup <- 1
istat <- 1
ifail <- 1

for (vf in vcf_files) {
  res <- tryCatch(
    analyze_one_gene(vf, accinfo, pheno),
    error = function(e) {
      geneID <- safe_gene_id(vf)
      failed[[ifail]] <<- data.frame(
        Gene = geneID,
        Status = "FAILED",
        Message = as.character(e$message),
        stringsAsFactors = FALSE
      )
      ifail <<- ifail + 1

      all_status[[istat]] <<- data.frame(
        Gene = geneID,
        Status = "FAILED",
        Accessions_after_filter = NA_integer_,
        Haplotypes = NA_integer_,
        stringsAsFactors = FALSE
      )
      istat <<- istat + 1
      NULL
    }
  )

  if (!is.null(res)) {
    all_enrich[[ie]] <- res$enrich
    ie <- ie + 1

    all_pair[[ip]] <- res$pair
    ip <- ip + 1

    all_superior[[isup]] <- res$superior
    isup <- isup + 1

    all_status[[istat]] <- res$status
    istat <- istat + 1
  }
}

enrich_all <- if (length(all_enrich) > 0) bind_rows(all_enrich) else data.frame()
pair_all <- if (length(all_pair) > 0) bind_rows(all_pair) else data.frame()
superior_all <- if (length(all_superior) > 0) bind_rows(all_superior) else data.frame()
status_all <- if (length(all_status) > 0) bind_rows(all_status) else data.frame()
failed_all <- if (length(failed) > 0) bind_rows(failed) else data.frame()

# =========================
# Global BH correction across all gene-haplotype tests
# =========================
if (nrow(enrich_all) > 0) {
  enrich_all$FDR <- NA_real_
  ok <- !is.na(enrich_all$p_hyper)
  enrich_all$FDR[ok] <- p.adjust(enrich_all$p_hyper[ok], method = "BH")

  enrich_all <- enrich_all %>%
    mutate(
      Guangdong_characteristic = (!is.na(freq_gd) & freq_gd > gd_freq_cutoff) &
        (!is.na(FDR) & FDR < fdr_cutoff)
    )
}

if (nrow(enrich_all) > 0 && nrow(superior_all) > 0) {
  final_all <- enrich_all %>%
    left_join(superior_all, by = c("Gene", "Hap")) %>%
    mutate(
      superior_haplotype = ifelse(is.na(superior_haplotype), FALSE, superior_haplotype),
      Candidate_type = case_when(
        Guangdong_characteristic & superior_haplotype ~ "Guangdong_characteristic_and_superior",
        Guangdong_characteristic & !superior_haplotype ~ "Guangdong_characteristic_only",
        !Guangdong_characteristic & superior_haplotype ~ "Superior_only",
        TRUE ~ "Others"
      )
    )
} else if (nrow(enrich_all) > 0) {
  final_all <- enrich_all %>%
    mutate(
      n = NA_integer_,
      mean_trait = NA_real_,
      sd_trait = NA_real_,
      superior_haplotype = FALSE,
      Candidate_type = ifelse(Guangdong_characteristic, "Guangdong_characteristic_only", "Others")
    )
} else {
  final_all <- data.frame()
}

# =========================
# Write global summaries
# =========================
if (nrow(enrich_all) > 0) {
  fwrite(enrich_all, file = file.path(outdir, "summary", "01_all_haplotype_enrichment.tsv"), sep = "\t")
}
if (nrow(pair_all) > 0) {
  fwrite(pair_all, file = file.path(outdir, "summary", "02_all_pairwise_ttest.tsv"), sep = "\t")
}
if (nrow(superior_all) > 0) {
  fwrite(superior_all, file = file.path(outdir, "summary", "03_all_superior_haplotypes.tsv"), sep = "\t")
}
if (nrow(final_all) > 0) {
  fwrite(final_all, file = file.path(outdir, "summary", "04_all_final_haplotype_summary.tsv"), sep = "\t")

  fwrite(
    final_all %>% filter(Guangdong_characteristic),
    file = file.path(outdir, "summary", "05_guangdong_characteristic_haplotypes.tsv"),
    sep = "\t"
  )

  fwrite(
    final_all %>% filter(superior_haplotype),
    file = file.path(outdir, "summary", "06_superior_haplotypes.tsv"),
    sep = "\t"
  )

  fwrite(
    final_all %>% filter(Guangdong_characteristic & superior_haplotype),
    file = file.path(outdir, "summary", "07_guangdong_characteristic_and_superior.tsv"),
    sep = "\t"
  )
}

if (nrow(status_all) > 0) {
  fwrite(status_all, file = file.path(outdir, "logs", "gene_status.tsv"), sep = "\t")
}
if (nrow(failed_all) > 0) {
  fwrite(failed_all, file = file.path(outdir, "logs", "failed_genes.tsv"), sep = "\t")
}

# =========================
# Write per-gene final summaries using global FDR
# =========================
if (nrow(final_all) > 0) {
  split_final <- split(final_all, final_all$Gene)
  for (g in names(split_final)) {
    gdir <- file.path(outdir, "per_gene", g)
    dir.create(gdir, recursive = TRUE, showWarnings = FALSE)
    fwrite(
      split_final[[g]],
      file = file.path(gdir, paste0(g, ".final_haplotype_summary.tsv")),
      sep = "\t"
    )
  }
}

message("[INFO] Batch analysis finished.")