#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(agricolae)
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
      "Usage: Rscript step3_haplotype_combination.R <step2_outdir> <accinfo.tsv> <pheno.tsv> <outdir>",
      "[trait_col=Chalkiness] [acc_id_col=ID] [subpop_col=Subpopulation]",
      "[gd_group=Guangdong_indica] [int_group=International_indica]",
      "[genes=Chalk5,Wx,OsDER1,OsATG8b] [major_combo_count_cutoff=3] [alpha=0.05]"
    )
  )
}

step2_outdir <- args[1]
accinfo_file <- args[2]
pheno_file <- args[3]
outdir <- args[4]

trait_col <- ifelse(length(args) >= 5, args[5], "Chalkiness")
acc_id_col <- ifelse(length(args) >= 6, args[6], "ID")
subpop_col <- ifelse(length(args) >= 7, args[7], "Subpopulation")
gd_group <- ifelse(length(args) >= 8, args[8], "Guangdong_indica")
int_group <- ifelse(length(args) >= 9, args[9], "International_indica")
genes_arg <- ifelse(length(args) >= 10, args[10], "Chalk5,Wx,OsDER1,OsATG8b")
major_combo_count_cutoff <- ifelse(length(args) >= 11, as.numeric(args[11]), 3)
alpha_level <- ifelse(length(args) >= 12, as.numeric(args[12]), 0.05)

genes <- strsplit(genes_arg, ",")[[1]]
genes <- trimws(genes)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE)

# =========================
# Helper functions
# =========================
read_one_gene_assignment <- function(step2_outdir, gene) {
  f <- file.path(step2_outdir, "per_gene", gene, paste0(gene, ".major_haplotype_assignment.tsv"))
  if (!file.exists(f)) {
    stop(paste("Missing file:", f))
  }

  df <- fread(f, data.table = FALSE)
  if (!all(c("Accession", "Hap") %in% colnames(df))) {
    stop(paste("Columns Accession/Hap not found in:", f))
  }

  df %>%
    transmute(
      Accession = as.character(Accession),
      !!gene := as.character(Hap)
    ) %>%
    distinct()
}

run_anova_duncan <- function(df_pop, combo_col = "ComboID", trait_col = "Trait", alpha = 0.05) {
  # Need at least 2 groups and enough observations
  n_groups <- dplyr::n_distinct(df_pop[[combo_col]])
  if (n_groups < 2) {
    return(list(
      anova = data.frame(),
      groups = data.frame(),
      means = data.frame()
    ))
  }

  fit <- aov(stats::as.formula(paste(trait_col, "~", combo_col)), data = df_pop)
  anova_df <- as.data.frame(anova(fit))
  anova_df$Term <- rownames(anova_df)
  rownames(anova_df) <- NULL

  duncan_out <- agricolae::duncan.test(
    fit,
    combo_col,
    alpha = alpha,
    group = TRUE,
    console = FALSE
  )

  groups_df <- duncan_out$groups
  groups_df$ComboID <- rownames(groups_df)
  rownames(groups_df) <- NULL

  means_df <- duncan_out$means
  means_df$ComboID <- rownames(means_df)
  rownames(means_df) <- NULL

  list(
    anova = anova_df,
    groups = groups_df,
    means = means_df
  )
}

make_population_boxplot <- function(df_pop, group_letters, pop_name, trait_col, outfile) {
  if (nrow(df_pop) == 0) return(NULL)

  stat_df <- df_pop %>%
    group_by(ComboID) %>%
    summarise(
      n = n(),
      mean_trait = mean(Trait, na.rm = TRUE),
      ymax = max(Trait, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(group_letters, by = "ComboID")

  combo_order <- stat_df %>%
    arrange(mean_trait) %>%
    pull(ComboID)

  df_pop$ComboID <- factor(df_pop$ComboID, levels = combo_order)
  stat_df$ComboID <- factor(stat_df$ComboID, levels = combo_order)

  y_range <- range(df_pop$Trait, na.rm = TRUE)
  y_span <- y_range[2] - y_range[1]
  if (!is.finite(y_span) || y_span == 0) y_span <- 1

  stat_df <- stat_df %>%
    mutate(
      label_y = ymax + 0.08 * y_span,
      groups = ifelse(is.na(groups), "", as.character(groups))
    )

  p <- ggplot(df_pop, aes(x = ComboID, y = Trait, fill = ComboID, color = ComboID)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 1.6, alpha = 0.75) +
    geom_text(
      data = stat_df,
      aes(x = ComboID, y = label_y, label = groups),
      inherit.aes = FALSE,
      size = 4
    ) +
    labs(
      title = pop_name,
      x = "Major haplotype combination",
      y = trait_col
    ) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.18)))

  ggsave(outfile, p, width = 8, height = 5)
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
# 1. Read major haplotype assignments for the four genes
# =========================
assign_list <- lapply(genes, function(g) read_one_gene_assignment(step2_outdir, g))
combo_df <- Reduce(function(x, y) full_join(x, y, by = "Accession"), assign_list)

# Keep only accessions with non-missing major haplotypes at all four genes
combo_df_complete <- combo_df %>%
  filter(if_all(all_of(genes), ~ !is.na(.) & . != ""))

# Build combination strings
combo_df_complete <- combo_df_complete %>%
  mutate(
    Combination = apply(
      combo_df_complete[, genes, drop = FALSE],
      1,
      function(z) paste(paste0(genes, "=", z), collapse = "|")
    )
  )

# Assign stable combination IDs by overall frequency
combo_id_map <- combo_df_complete %>%
  count(Combination, name = "Global_count", sort = TRUE) %>%
  mutate(ComboID = sprintf("HC%03d", seq_len(n())))

combo_df_complete <- combo_df_complete %>%
  left_join(combo_id_map, by = "Combination")

# Merge population and phenotype
combo_merged <- combo_df_complete %>%
  left_join(accinfo, by = c("Accession" = acc_id_col)) %>%
  left_join(pheno, by = c("Accession" = acc_id_col))

# Restrict to the two target populations
combo_merged <- combo_merged %>%
  filter(.data[[subpop_col]] %in% c(gd_group, int_group)) %>%
  mutate(
    Population = dplyr::case_when(
      .data[[subpop_col]] == gd_group ~ "Guangdong_indica",
      .data[[subpop_col]] == int_group ~ "International_indica",
      TRUE ~ as.character(.data[[subpop_col]])
    ),
    Trait = as.numeric(.data[[trait_col]])
  )

# =========================
# 2. Frequency tables
# =========================
freq_by_pop <- combo_merged %>%
  group_by(Population, ComboID, Combination) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Population) %>%
  mutate(Frequency = Count / sum(Count)) %>%
  ungroup()

global_combo_count <- combo_merged %>%
  count(ComboID, Combination, name = "Global_count") %>%
  mutate(Major_combination = Global_count > major_combo_count_cutoff)

major_combo_defs <- global_combo_count %>%
  filter(Major_combination)

freq_by_pop <- freq_by_pop %>%
  left_join(global_combo_count, by = c("ComboID", "Combination"))

# Keep major combinations for downstream analysis
combo_major <- combo_merged %>%
  semi_join(major_combo_defs, by = c("ComboID", "Combination"))

# =========================
# 3. Phenotype analysis within each population
# =========================
anova_list <- list()
duncan_list <- list()
pheno_summary_list <- list()

for (pop_name in c("Guangdong_indica", "International_indica")) {
  df_pop <- combo_major %>%
    filter(Population == pop_name, !is.na(Trait), ComboID != "")

  if (nrow(df_pop) == 0) next

  # Per-combination phenotype summary
  pheno_summary <- df_pop %>%
    group_by(Population, ComboID, Combination) %>%
    summarise(
      n = n(),
      mean_trait = mean(Trait, na.rm = TRUE),
      sd_trait = sd(Trait, na.rm = TRUE),
      .groups = "drop"
    )

  # ANOVA + Duncan
  stats_out <- run_anova_duncan(df_pop, combo_col = "ComboID", trait_col = "Trait", alpha = alpha_level)

  anova_df <- stats_out$anova
  if (nrow(anova_df) > 0) {
    anova_df$Population <- pop_name
    anova_list[[pop_name]] <- anova_df
  }

  duncan_groups <- stats_out$groups
  if (nrow(duncan_groups) > 0) {
    duncan_groups$Population <- pop_name
    duncan_groups <- duncan_groups %>%
      rename(Duncan_group = groups)
    duncan_list[[pop_name]] <- duncan_groups
  } else {
    duncan_groups <- data.frame(
      ComboID = character(),
      Duncan_group = character(),
      Population = character(),
      stringsAsFactors = FALSE
    )
  }

  pheno_summary <- pheno_summary %>%
    left_join(duncan_groups %>% select(ComboID, Duncan_group), by = "ComboID")

  pheno_summary_list[[pop_name]] <- pheno_summary

  # Boxplot with grouping letters
  plot_groups <- duncan_groups %>%
    transmute(ComboID, groups = Duncan_group)

  make_population_boxplot(
    df_pop = df_pop,
    group_letters = plot_groups,
    pop_name = pop_name,
    trait_col = trait_col,
    outfile = file.path(outdir, "plots", paste0(pop_name, ".major_haplotype_combination.", trait_col, ".boxplot.pdf"))
  )
}

anova_all <- if (length(anova_list) > 0) bind_rows(anova_list) else data.frame()
duncan_all <- if (length(duncan_list) > 0) bind_rows(duncan_list) else data.frame()
pheno_summary_all <- if (length(pheno_summary_list) > 0) bind_rows(pheno_summary_list) else data.frame()

# =========================
# 4. Export tables
# =========================
fwrite(combo_df_complete, file.path(outdir, "tables", "01_accession_haplotype_combinations.tsv"), sep = "\t")
fwrite(combo_id_map, file.path(outdir, "tables", "02_haplotype_combination_definitions.tsv"), sep = "\t")
fwrite(freq_by_pop, file.path(outdir, "tables", "03_haplotype_combination_frequency_by_population.tsv"), sep = "\t")
fwrite(major_combo_defs, file.path(outdir, "tables", "04_major_haplotype_combinations.tsv"), sep = "\t")
fwrite(combo_major, file.path(outdir, "tables", "05_major_haplotype_combination_merged.tsv"), sep = "\t")

if (nrow(pheno_summary_all) > 0) {
  fwrite(pheno_summary_all, file.path(outdir, "tables", paste0("06_major_haplotype_combination_", trait_col, "_summary.tsv")), sep = "\t")
}
if (nrow(anova_all) > 0) {
  fwrite(anova_all, file.path(outdir, "tables", paste0("07_major_haplotype_combination_", trait_col, "_anova.tsv")), sep = "\t")
}
if (nrow(duncan_all) > 0) {
  fwrite(duncan_all, file.path(outdir, "tables", paste0("08_major_haplotype_combination_", trait_col, "_duncan_groups.tsv")), sep = "\t")
}

message("[INFO] Finished combination analysis.")