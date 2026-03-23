# =========================================================
# Mixed-effects model analysis and BLUE distribution plotting
# =========================================================
#
# This script:
#   1. Reads a phenotype file containing multiple replicate values per sample
#   2. Reshapes the data into long format
#   3. Fits a mixed-effects model with SampleID as a fixed effect
#      and replicate index as a random effect
#   4. Estimates LSMeans (used here as BLUE values) for each sample
#   5. Exports the BLUE results and a histogram of their distribution
#
# Required arguments:
#   input_file     : Input phenotype file in TXT format
#   output_prefix  : Prefix for output files
#   trait_name     : Trait name used in plot labels and titles
# =========================================================

# Load required libraries
library(argparser, quietly = TRUE)
library(magrittr)
library(lsmeans)
library(lme4)
library(ggplot2)

# ---------------------------------------------------------
# Parse command-line arguments
# ---------------------------------------------------------
p <- arg_parser("Script for mixed-effects model analysis and histogram generation")

p <- add_argument(
  p, "input_file",
  help = "Path to the input file (TXT format)",
  type = "character"
)

p <- add_argument(
  p, "output_prefix",
  help = "Prefix for the output files",
  type = "character"
)

p <- add_argument(
  p, "trait_name",
  help = "Name of the trait to analyze",
  type = "character"
)

args <- parse_args(p)

# Assign input arguments
input_file_path <- args$input_file
output_file_prefix <- args$output_prefix
trait_name <- args$trait_name

# Define output file names
lsmeans_output_file_path <- paste0(output_file_prefix, "_Trait_Blue.txt")
histogram_output_file_path <- paste0(output_file_prefix, "_Trait_Blue_histogram.pdf")

# ---------------------------------------------------------
# Read input data
# ---------------------------------------------------------
# The input file is expected to contain:
#   - a header row
#   - the first column as SampleID
#   - the remaining columns as replicate trait values
input_lines <- readLines(input_file_path)

# Identify the maximum number of replicate values across all samples
# This is used to standardize row lengths when reshaping the data
max_length <- max(sapply(input_lines[-1], function(line) {
  length(strsplit(line, "\\s+")[[1]]) - 1
}))

# ---------------------------------------------------------
# Reshape data into long format
# ---------------------------------------------------------
# Create an empty data frame to store:
#   SampleID    : sample identifier
#   rep         : replicate index
#   Trait_value : observed trait value
output_df <- data.frame(
  SampleID = character(),
  rep = integer(),
  Trait_value = character(),
  stringsAsFactors = FALSE
)

# Process each sample row
for (line in input_lines[-1]) {
  data <- strsplit(line, "\\s+")[[1]]
  sample_id <- data[1]

  # Fill missing replicate positions with 0 so that all samples
  # have the same number of replicate entries
  values <- c(data[-1], rep("0", max_length - (length(data) - 1)))

  # Append sample data in long format
  output_df <- rbind(
    output_df,
    data.frame(
      SampleID = sample_id,
      rep = seq_along(values),
      Trait_value = values,
      stringsAsFactors = FALSE
    )
  )
}

# Convert columns to appropriate data types
output_df$Trait_value <- as.numeric(output_df$Trait_value)
output_df$rep <- as.integer(output_df$rep)

# ---------------------------------------------------------
# Fit mixed-effects model
# ---------------------------------------------------------
# Model structure:
#   Trait_value ~ SampleID + (1 | rep)
#
# where:
#   - SampleID is treated as a fixed effect
#   - rep is treated as a random effect
m1 <- lmer(Trait_value ~ SampleID + (1 | rep), data = output_df)

# Print model summary to screen
summary(m1)

# ---------------------------------------------------------
# Estimate LSMeans / BLUE values
# ---------------------------------------------------------
# Calculate estimated marginal means for each SampleID
re <- lsmeans(m1, "SampleID")

# Export BLUE results
write.table(
  as.data.frame(re),
  lsmeans_output_file_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Convert results to a data frame for downstream plotting
re_df <- as.data.frame(re)

# ---------------------------------------------------------
# Summarize BLUE distribution
# ---------------------------------------------------------
# Extract BLUE values and remove missing entries
dat <- na.omit(re_df$lsmean)

# Calculate the mean of BLUE values
m <- round(mean(dat), 4)

# Define a function for population standard deviation
stdev_p <- function(x) {
  sqrt(mean((x - mean(x))^2))
}

# Calculate population standard deviation of BLUE values
v <- round(stdev_p(dat), 4)

# ---------------------------------------------------------
# Plot histogram of BLUE values
# ---------------------------------------------------------
p1 <- ggplot(re_df, aes(x = lsmean, fill = "orange")) +
  geom_histogram(linewidth = 0.5, color = "black") +
  ggtitle(paste(trait_name, "mean:", m, "sd:", v, sep = " ")) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  labs(x = trait_name)

# Save histogram to PDF
ggsave(
  histogram_output_file_path,
  plot = p1,
  width = 5,
  height = 5
)