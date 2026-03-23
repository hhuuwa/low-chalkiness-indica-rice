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
      "Usage: Rscript step2_major_haplotype.R <functional_vcf_dir> <accinfo.tsv> <pheno.tsv> <outdir>",
      "[trait_col=Chalkiness] [acc_id_col=ID] [subpop_col=Subpopulation]",
      "[hap_prefix=H] [major_freq_cutoff=0.05] [major_count_cutoff=16] [ttest_p_cutoff=0.05]"
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
hap_prefix <- ifelse(length(args) >= 8, args[8], "H")
major_freq_cutoff <- ifelse(length(args) >= 9, as.numeric(args[9]), 0.05)
major_count_cutoff <- ifelse(length(args) >= 10, as.numeric(args[10]), 16)
ttest_p_cutoff <- ifelse(length(args) >= 11, as.numeric(args[11]), 0.05)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "per_gene"), recursive = TRUE, showWarnings = FALSE)

# =========================
# Helper functions
# =========================
safe_gene_id <- function(vcf_file) {
  x <- basename(vcf_file)
  x <- sub("\\.functional\\.vcf\\.gz$", "", x, ignore.case = TRUE)
  x <- sub("\\.vcf\\.gz$", "", x, ignore.case = TRUE)
  x <- sub("\\.vcf$", "", x, ignore.case = TRUE)
  x
}

strip_meta_rows <- function(x) {
  df <- as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(df) <= 4) return(NULL)
  df[-(1:4), , drop = FALSE]
}

make_pairwise_tests <- function(df, p_cutoff = 0.05) {
  # df must contain Hap and Value
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
          Significant = tt$p.value < p_cutoff,
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
      Significant = logical(),
      stringsAsFactors = FALSE
    )
  } else {
    pair_df <- bind_rows(pair_res)
  }

  list(hap_stat = hap_stat, pair_df = pair_df)
}

plot_major_hap_distribution_within_hap <- function(df, subpop_col, geneID, outfile) {
  tmp <- df %>%
    filter(!is.na(Hap), Hap != "", !is.na(.data[[subpop_col]])) %>%
    group_by(Hap, Group = .data[[subpop_col]]) %>%
    summarise(Freq = n(), .groups = "drop") %>%
    group_by(Hap) %>%
    mutate(
      Sum_Freq = sum(Freq),
      Percent = Freq / Sum_Freq,
      Label = paste0(round(100 * Percent, 1), "% (", Freq, ")")
    ) %>%
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
    labs(
      title = geneID,
      x = "",
      y = "Percent within major haplotype"
    ) +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  ggsave(outfile, p, width = 7, height = 5)
}

plot_major_hap_distribution_within_group <- function(df, subpop_col, geneID, outfile) {
  tmp <- df %>%
    filter(!is.na(Hap), Hap != "", !is.na(.data[[subpop_col]])) %>%
    group_by(Group = .data[[subpop_col]], Hap) %>%
    summarise(Freq = n(), .groups = "drop") %>%
    group_by(Group) %>%
    mutate(
      Sum_Freq = sum(Freq),
      Percent = Freq / Sum_Freq,
      Label = paste0(round(100 * Percent, 1), "% (", Freq, ")")
    ) %>%
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
    labs(
      title = geneID,
      x = "",
      y = "Percent within subpopulation"
    ) +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  ggsave(outfile, p, width = 7, height = 5)
}

plot_major_hap_trait_boxplot <- function(df, trait_col, geneID, outfile) {
  df_plot <- df %>%
    filter(!is.na(Hap), Hap != "", !is.na(.data[[trait_col]])) %>%
    transmute(
      Hap = as.character(Hap),
      Value = as.numeric(.data[[trait_col]])
    )

  if (nrow(df_plot) == 0) return(NULL)

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
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.75) +
    labs(
      title = geneID,
      x = "Major haplotype",
      y = trait_col
    ) +
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

  ggsave(outfile, p, width = 7, height = 5)
}

# =========================
# Read input tables
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
  pattern = "\\.functional\\.vcf\\.gz$",
  full.names = TRUE
)

if (length(vcf_files) == 0) {
  stop("No functional VCF files found.")
}

# =========================
# Batch analysis
# =========================
all_freq <- list()
all_major <- list()
all_status <- list()
all_pairwise <- list()
all_pheno_summary <- list()

i1 <- 1
i2 <- 1
i3 <- 1
i4 <- 1
i5 <- 1

for (vf in vcf_files) {
  geneID <- safe_gene_id(vf)
  message("[INFO] Processing ", geneID)

  gdir <- file.path(outdir, "per_gene", geneID)
  dir.create(gdir, recursive = TRUE, showWarnings = FALSE)

  res <- tryCatch({
    vcf <- import_vcf(vf)

    hapResult <- vcf2hap(
      vcf,
      hapPrefix = hap_prefix,
      hetero_remove = TRUE,
      na_drop = FALSE
    )

    hap_df <- strip_meta_rows(hapResult)
    if (is.null(hap_df) || !all(c("Accession", "Hap") %in% colnames(hap_df))) {
      stop("No usable haplotype rows produced.")
    }

    hap_assign <- hap_df[, c("Accession", "Hap"), drop = FALSE] %>%
      mutate(
        Accession = as.character(Accession),
        Hap = as.character(Hap)
      ) %>%
      filter(Accession != "", Hap != "")

    if (nrow(hap_assign) == 0) stop("No valid haplotype assignments.")

    n_total <- nrow(hap_assign)

    hap_freq <- hap_assign %>%
      count(Hap, name = "Count") %>%
      mutate(
        Gene = geneID,
        Frequency = Count / n_total,
        Major_haplotype = Frequency > major_freq_cutoff & Count >= major_count_cutoff
      ) %>%
      select(Gene, Hap, Count, Frequency, Major_haplotype)

    major_haps <- hap_freq %>%
      filter(Major_haplotype)

    hap_assign_major <- hap_assign %>%
      semi_join(major_haps, by = "Hap")

    # geneHapR summaries
    write.hap(hapResult, file = file.path(gdir, paste0(geneID, ".hapResult.tsv")), sep = "\t")
    hapSummary <- hap_summary(hapResult, hapPrefix = hap_prefix)
    write.hap(hapSummary, file = file.path(gdir, paste0(geneID, ".hapSummary.tsv")), sep = "\t")

    # Basic exports
    fwrite(hap_assign, file.path(gdir, paste0(geneID, ".hap_assignment.tsv")), sep = "\t")
    fwrite(hap_freq, file.path(gdir, paste0(geneID, ".hap_frequency.tsv")), sep = "\t")
    fwrite(major_haps, file.path(gdir, paste0(geneID, ".major_haplotypes.tsv")), sep = "\t")
    fwrite(hap_assign_major, file.path(gdir, paste0(geneID, ".major_haplotype_assignment.tsv")), sep = "\t")

    # Frequency plot
    p_freq <- ggplot(hap_freq, aes(x = reorder(Hap, -Count), y = Frequency)) +
      geom_col() +
      geom_hline(yintercept = major_freq_cutoff, linetype = 2) +
      geom_text(
        aes(label = paste0(Count, " (", round(Frequency * 100, 1), "%)")),
        vjust = -0.3,
        size = 3
      ) +
      labs(title = geneID, x = "Haplotype", y = "Frequency") +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )

    ggsave(file.path(gdir, paste0(geneID, ".hap_frequency.pdf")), p_freq, width = 6, height = 4)

    # Merge with subgroup and phenotype
    merged_major <- hap_assign_major %>%
      left_join(accinfo, by = c("Accession" = acc_id_col)) %>%
      left_join(pheno, by = c("Accession" = acc_id_col))

    fwrite(
      merged_major,
      file.path(gdir, paste0(geneID, ".major_haplotype_merged.tsv")),
      sep = "\t"
    )

    # Group distribution
    plot_major_hap_distribution_within_hap(
      merged_major,
      subpop_col = subpop_col,
      geneID = geneID,
      outfile = file.path(gdir, paste0(geneID, ".major_hap.class1.pdf"))
    )

    plot_major_hap_distribution_within_group(
      merged_major,
      subpop_col = subpop_col,
      geneID = geneID,
      outfile = file.path(gdir, paste0(geneID, ".major_hap.class2.pdf"))
    )

    # Phenotype summary + pairwise t-test
    trait_df <- merged_major %>%
      filter(!is.na(Hap), Hap != "", !is.na(.data[[trait_col]])) %>%
      transmute(
        Hap = as.character(Hap),
        Accession = as.character(Accession),
        Value = as.numeric(.data[[trait_col]])
      )

    if (nrow(trait_df) > 0) {
      tt <- make_pairwise_tests(trait_df, p_cutoff = ttest_p_cutoff)
      pheno_summary <- tt$hap_stat %>%
        mutate(Gene = geneID, .before = 1)
      pairwise_df <- tt$pair_df %>%
        mutate(Gene = geneID, .before = 1)
    } else {
      pheno_summary <- data.frame(
        Gene = geneID,
        Hap = character(),
        n = integer(),
        mean_trait = numeric(),
        sd_trait = numeric(),
        stringsAsFactors = FALSE
      )
      pairwise_df <- data.frame(
        Gene = geneID,
        Hap1 = character(),
        Hap2 = character(),
        n1 = integer(),
        n2 = integer(),
        mean1 = numeric(),
        mean2 = numeric(),
        p_value = numeric(),
        Significant = logical(),
        stringsAsFactors = FALSE
      )
    }

    fwrite(
      pheno_summary,
      file.path(gdir, paste0(geneID, ".major_haplotype_pheno_summary.tsv")),
      sep = "\t"
    )

    fwrite(
      pairwise_df,
      file.path(gdir, paste0(geneID, ".major_haplotype_pairwise_ttest.tsv")),
      sep = "\t"
    )

    plot_major_hap_trait_boxplot(
      merged_major,
      trait_col = trait_col,
      geneID = geneID,
      outfile = file.path(gdir, paste0(geneID, ".major_hap.", trait_col, ".boxplot.pdf"))
    )

    list(
      hap_freq = hap_freq,
      major_haps = major_haps,
      pheno_summary = pheno_summary,
      pairwise_df = pairwise_df,
      status = data.frame(
        Gene = geneID,
        N_accessions = n_total,
        N_haplotypes = dplyr::n_distinct(hap_assign$Hap),
        N_major_haplotypes = nrow(major_haps),
        stringsAsFactors = FALSE
      )
    )
  }, error = function(e) {
    list(
      hap_freq = NULL,
      major_haps = NULL,
      pheno_summary = NULL,
      pairwise_df = NULL,
      status = data.frame(
        Gene = geneID,
        N_accessions = NA,
        N_haplotypes = NA,
        N_major_haplotypes = NA,
        Message = as.character(e$message),
        stringsAsFactors = FALSE
      )
    )
  })

  if (!is.null(res$hap_freq)) {
    all_freq[[i1]] <- res$hap_freq
    i1 <- i1 + 1
  }
  if (!is.null(res$major_haps)) {
    all_major[[i2]] <- res$major_haps
    i2 <- i2 + 1
  }
  if (!is.null(res$pheno_summary)) {
    all_pheno_summary[[i3]] <- res$pheno_summary
    i3 <- i3 + 1
  }
  if (!is.null(res$pairwise_df)) {
    all_pairwise[[i4]] <- res$pairwise_df
    i4 <- i4 + 1
  }
  all_status[[i5]] <- res$status
  i5 <- i5 + 1
}

freq_all <- if (length(all_freq) > 0) bind_rows(all_freq) else data.frame()
major_all <- if (length(all_major) > 0) bind_rows(all_major) else data.frame()
pheno_all <- if (length(all_pheno_summary) > 0) bind_rows(all_pheno_summary) else data.frame()
pairwise_all <- if (length(all_pairwise) > 0) bind_rows(all_pairwise) else data.frame()
status_all <- if (length(all_status) > 0) bind_rows(all_status) else data.frame()

if (nrow(freq_all) > 0) {
  fwrite(freq_all, file.path(outdir, "all_haplotype_frequencies.tsv"), sep = "\t")
}
if (nrow(major_all) > 0) {
  fwrite(major_all, file.path(outdir, "all_major_haplotypes.tsv"), sep = "\t")
}
if (nrow(pheno_all) > 0) {
  fwrite(pheno_all, file.path(outdir, paste0("all_major_haplotype_", trait_col, "_summary.tsv")), sep = "\t")
}
if (nrow(pairwise_all) > 0) {
  fwrite(pairwise_all, file.path(outdir, paste0("all_major_haplotype_", trait_col, "_pairwise_ttest.tsv")), sep = "\t")
}
if (nrow(status_all) > 0) {
  fwrite(status_all, file.path(outdir, "run_status.tsv"), sep = "\t")
}

message("[INFO] Finished.")