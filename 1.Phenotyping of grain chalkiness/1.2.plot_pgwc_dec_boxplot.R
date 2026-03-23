#!/usr/bin/env Rscript

library(ggplot2)
library(ggpubr)
library(cowplot)

#---------------------------
# input files
#---------------------------
pgwc_file <- "311.PGWC.phe"
dec_file  <- "311.DEC.phe"
out_file  <- "PGWC_DEC.diff.pdf"

#---------------------------
# read data
#---------------------------
pgwc <- read.table(pgwc_file, sep = "\t", header = TRUE, check.names = FALSE)
dec  <- read.table(dec_file,  sep = "\t", header = TRUE, check.names = FALSE)

#---------------------------
# settings
#---------------------------
mycol <- c("#37BEDC", "#FE6500")

my_theme <- theme(
  panel.grid = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  axis.text.x = element_text(size = 12, face = "bold", hjust = 0.5),
  axis.text.y = element_text(size = 12, face = "bold", hjust = 0.5),
  axis.title.y = element_text(size = 13, face = "bold", vjust = 0.2),
  plot.title = element_text(hjust = 0.5),
  axis.title.x = element_blank()
)

#---------------------------
# function
#---------------------------
make_plot <- function(dat, yvar) {
  if (!"subpop" %in% colnames(dat)) {
    stop("Column 'subpop' not found.")
  }
  if (!yvar %in% colnames(dat)) {
    stop(paste("Column", yvar, "not found."))
  }

  subdat <- dat[, c("subpop", yvar)]
  colnames(subdat) <- c("subpop", "trait")
  subdat <- subdat[!is.na(subdat$subpop) & !is.na(subdat$trait), ]

  grp <- unique(subdat$subpop)
  if (length(grp) != 2) {
    stop(paste("Trait", yvar, ": subpop must contain exactly 2 groups. Current group number =", length(grp)))
  }

  tt <- t.test(trait ~ subpop, data = subdat)
  ptxt <- format(tt$p.value, scientific = TRUE, digits = 3)

  p <- ggboxplot(
    subdat,
    x = "subpop",
    y = "trait",
    color = "subpop",
    palette = mycol,
    legend = "none",
    add = "jitter",
    add.params = list(size = 3, jitter = 0.1, alpha = 0.6),
    outlier.shape = NA,
    main = paste0("t-test: p-value = ", ptxt)
  ) +
    ylab(yvar) +
    my_theme

  return(p)
}

#---------------------------
# make plots
#---------------------------
p1 <- make_plot(pgwc, "PGWC1")
p2 <- make_plot(pgwc, "PGWC2")
p3 <- make_plot(pgwc, "PGWC_blue")

p4 <- make_plot(dec, "DEC1")
p5 <- make_plot(dec, "DEC2")
p6 <- make_plot(dec, "DEC_BLUE")

#---------------------------
# combine
#---------------------------
p_all <- plot_grid(
  p1, p2, p3,
  p4, p5, p6,
  labels = c("A", "B", "C", "D", "E", "F"),
  ncol = 3,
  nrow = 2
)

ggsave(out_file, plot = p_all, width = 14, height = 10)