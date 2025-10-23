#!/usr/bin/env Rscript
#
# Script: Identify iHS selection candidates across all chromosomes
# Purpose: Extract significant iHS signals and compare X vs autosomes
#
# Author: MIT IHH PIB Project
# Date: 2025-01-23

library(tidyverse)
library(scales)

# Setup paths
base_dir <- "/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib"
ihs_dir <- file.path(base_dir, "results/ihs")
output_dir <- file.path(base_dir, "results/analysis")

# Create output directories
dir.create(file.path(output_dir, "candidates"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)

cat("===========================================\n")
cat("iHS Selection Candidate Identification\n")
cat("===========================================\n\n")

# Define chromosomes
chroms <- c(1:22, "X")

cat("Step 1: Processing each chromosome separately...\n")

# Define thresholds based on literature recommendations
# Voight et al. 2006: standardized iHS ~ N(0,1) under neutrality
# Common thresholds:
#   |iHS| > 2.0  ≈ top ~2.3% (2 std dev from mean)
#   |iHS| > 2.5  ≈ top ~0.6% (used in Paul et al. 2024)
#   |iHS| > 3.0  ≈ top ~0.13% (very stringent)
thresholds <- c(
  "liberal" = 2.0,     # Voight et al. 2006 suggestion
  "moderate" = 2.5,    # Paul et al. 2024, more stringent
  "stringent" = 3.0    # Very stringent, strongest signals
)

# Will also compute p-values from standard normal distribution
# and empirical percentiles for comparison

# Create directory for per-chromosome candidates
dir.create(file.path(output_dir, "candidates/per_chromosome"),
           recursive = TRUE, showWarnings = FALSE)

# Process each chromosome and extract candidates
chr_stats_list <- list()
all_candidates <- list()

for (chr in chroms) {
  file <- file.path(ihs_dir, paste0("ALL.chr", chr, ".ihs.tsv"))

  if (!file.exists(file)) {
    cat(sprintf("  WARNING: Missing file for chr%s\n", chr))
    next
  }

  cat(sprintf("  Processing chr%s...", chr))

  # Read chromosome data
  chr_data <- read_tsv(file, col_types = cols(), show_col_types = FALSE) %>%
    separate(Location, into = c("chr", "pos", "ref", "alt"),
             sep = ":", remove = FALSE, convert = TRUE)

  cat(sprintf(" %s variants", format(nrow(chr_data), big.mark = ",")))

  # Add p-values based on standard normal assumption
  # Two-tailed p-value: P(|Z| > |iHS|) where Z ~ N(0,1)
  chr_data <- chr_data %>%
    mutate(
      p_value = 2 * pnorm(-abs(`Std iHS`)),
      neg_log10_p = -log10(p_value)
    )

  # Calculate statistics for this chromosome
  chr_stats_list[[as.character(chr)]] <- tibble(
    chr = chr,
    n_variants = nrow(chr_data),
    mean_iHS = mean(chr_data$iHS, na.rm = TRUE),
    median_iHS = median(chr_data$iHS, na.rm = TRUE),
    sd_iHS = sd(chr_data$iHS, na.rm = TRUE),
    mean_std_iHS = mean(chr_data$`Std iHS`, na.rm = TRUE),
    median_std_iHS = median(chr_data$`Std iHS`, na.rm = TRUE),
    sd_std_iHS = sd(chr_data$`Std iHS`, na.rm = TRUE),
    max_abs_std_iHS = max(abs(chr_data$`Std iHS`), na.rm = TRUE),
    q99 = quantile(abs(chr_data$`Std iHS`), 0.99, na.rm = TRUE),
    q995 = quantile(abs(chr_data$`Std iHS`), 0.995, na.rm = TRUE),
    q999 = quantile(abs(chr_data$`Std iHS`), 0.999, na.rm = TRUE)
  )

  # Extract candidates at each threshold for this chromosome
  for (level in names(thresholds)) {
    thresh <- thresholds[level]

    candidates <- chr_data %>%
      filter(abs(`Std iHS`) >= thresh) %>%
      arrange(desc(abs(`Std iHS`)))

    if (nrow(candidates) > 0) {
      # Save per-chromosome candidates
      outfile <- file.path(output_dir, "candidates/per_chromosome",
                          sprintf("chr%s_candidates_%s.tsv", chr, level))
      write_tsv(candidates, outfile)

      # Store for combined file
      candidates_key <- paste0(level, "_", chr)
      all_candidates[[candidates_key]] <- candidates
    }
  }

  cat(sprintf(" - Done\n"))
}

cat("\n")

# Combine statistics across chromosomes
chr_stats <- bind_rows(chr_stats_list) %>%
  mutate(chr = factor(chr, levels = c(1:22, "X")))

cat("Step 2: Calculating summary statistics...\n")

# Add chromosome type to stats
chr_stats <- chr_stats %>%
  mutate(chr_type = if_else(chr == "X", "X chromosome", "Autosome"))

# Overall statistics by chromosome type
overall_stats <- chr_stats %>%
  group_by(chr_type) %>%
  summarise(
    n_chromosomes = n(),
    total_variants = sum(n_variants),
    mean_variants_per_chr = mean(n_variants),
    mean_std_iHS = mean(mean_std_iHS),
    sd_std_iHS = mean(sd_std_iHS),
    mean_q99 = mean(q99),
    mean_q995 = mean(q995),
    .groups = "drop"
  )

cat("\n=== Distribution Check (Voight et al. 2006) ===\n")
cat("Standardized iHS should have mean ~0, variance ~1 under neutrality\n\n")

cat("Overall Statistics by Chromosome Type:\n")
print(overall_stats, n = Inf)

cat("\n=== Empirical Percentiles (|Std iHS|) ===\n")
cat("These show the empirical thresholds for top 1%, 0.5%, 0.1% of SNPs\n")
chr_percentiles <- chr_stats %>%
  select(chr, chr_type, q99, q995, q999) %>%
  arrange(chr)
print(chr_percentiles, n = Inf)

cat("\n=== Full Per-chromosome Statistics ===\n")
print(chr_stats, n = Inf)

write_tsv(chr_stats, file.path(output_dir, "candidates/chromosome_statistics.tsv"))
cat("\n  Saved: candidates/chromosome_statistics.tsv\n")

# Save empirical percentiles summary
percentile_summary <- chr_stats %>%
  select(chr, chr_type, n_variants, q99, q995, q999) %>%
  arrange(chr)
write_tsv(percentile_summary,
          file.path(output_dir, "candidates/empirical_percentiles.tsv"))
cat("  Saved: candidates/empirical_percentiles.tsv\n\n")

# Check if distributions are approximately normal
cat("=== Normality Check ===\n")
cat("For standard normal, mean should be ~0, SD should be ~1\n")
normality_check <- chr_stats %>%
  select(chr, chr_type, mean_std_iHS, sd_std_iHS) %>%
  mutate(
    mean_deviation = abs(mean_std_iHS),
    sd_deviation = abs(sd_std_iHS - 1)
  )
print(normality_check, n = Inf)
cat("\nNote: Large deviations suggest non-normal distribution (Salazar-Tortosa et al. 2023)\n")
cat("      Consider using empirical percentiles rather than fixed thresholds\n\n")

cat("Step 3: Combining candidates across chromosomes...\n")

# Combine candidates at each threshold level
candidates_list <- map(names(thresholds), function(level) {
  thresh <- thresholds[level]

  # Get all candidates for this level across all chromosomes
  level_candidates <- map_df(chroms, function(chr) {
    key <- paste0(level, "_", chr)
    if (key %in% names(all_candidates)) {
      return(all_candidates[[key]])
    }
    return(NULL)
  })

  # Add chromosome type
  level_candidates <- level_candidates %>%
    mutate(chr_type = if_else(chr == "X", "X chromosome", "Autosome"))

  cat(sprintf("  %s threshold (|Std iHS| >= %.1f): %s variants\n",
              str_to_title(level), thresh,
              format(nrow(level_candidates), big.mark = ",")))

  # Count by chromosome type
  by_type <- level_candidates %>%
    count(chr_type) %>%
    mutate(pct = n / sum(n) * 100)

  if (nrow(by_type) > 0) {
    auto_row <- by_type %>% filter(chr_type == "Autosome")
    x_row <- by_type %>% filter(chr_type == "X chromosome")

    if (nrow(auto_row) > 0) {
      cat(sprintf("    - Autosomes: %s (%.1f%%)\n",
                  format(auto_row$n, big.mark = ","),
                  auto_row$pct))
    }
    if (nrow(x_row) > 0) {
      cat(sprintf("    - X chromosome: %s (%.1f%%)\n",
                  format(x_row$n, big.mark = ","),
                  x_row$pct))
    }
  }

  # Save combined candidates
  outfile <- file.path(output_dir, "candidates",
                       sprintf("ALL_candidates_%s.tsv", level))
  write_tsv(level_candidates, outfile)
  cat(sprintf("    Saved: candidates/ALL_candidates_%s.tsv\n\n", level))

  return(level_candidates)
})

names(candidates_list) <- names(thresholds)

cat("Step 4: Testing for X chromosome enrichment...\n")

# Test if X chromosome is enriched for selection signals
# Using moderate threshold (|Std iHS| >= 2.5)
candidates_moderate <- candidates_list$moderate

# Contingency table
n_candidates_x <- sum(candidates_moderate$chr == "X")
n_candidates_auto <- sum(candidates_moderate$chr != "X")
n_total_x <- chr_stats %>% filter(chr == "X") %>% pull(n_variants) %>% sum()
n_total_auto <- chr_stats %>% filter(chr != "X") %>% pull(n_variants) %>% sum()

contingency <- matrix(c(
  n_candidates_x, n_total_x - n_candidates_x,
  n_candidates_auto, n_total_auto - n_candidates_auto
), nrow = 2, byrow = TRUE,
dimnames = list(c("X", "Autosomes"), c("Candidate", "Non-candidate")))

# Fisher's exact test
fisher_result <- fisher.test(contingency)

cat(sprintf("  X chromosome candidates: %s / %s (%.2f%%)\n",
            n_candidates_x, n_total_x,
            n_candidates_x / n_total_x * 100))
cat(sprintf("  Autosome candidates: %s / %s (%.2f%%)\n",
            n_candidates_auto, n_total_auto,
            n_candidates_auto / n_total_auto * 100))
cat(sprintf("  Odds ratio: %.3f\n", fisher_result$estimate))
cat(sprintf("  P-value: %.2e\n\n", fisher_result$p.value))

# Save test results
enrichment_results <- tibble(
  test = "Fisher's exact test",
  threshold = "moderate (|Std iHS| >= 2.5)",
  x_candidates = n_candidates_x,
  x_total = n_total_x,
  x_pct = n_candidates_x / n_total_x * 100,
  auto_candidates = n_candidates_auto,
  auto_total = n_total_auto,
  auto_pct = n_candidates_auto / n_total_auto * 100,
  odds_ratio = as.numeric(fisher_result$estimate),
  p_value = fisher_result$p.value
)

write_tsv(enrichment_results,
          file.path(output_dir, "candidates/x_enrichment_test.tsv"))
cat("  Saved: candidates/x_enrichment_test.tsv\n\n")

cat("Step 5: Creating visualizations...\n")

# For visualization, we'll use the combined moderate candidates
# and sample from full data if needed for density plots

# Plot 1: Candidate counts by chromosome
cat("  Plotting candidate counts by chromosome...\n")

candidates_by_chr <- candidates_list$moderate %>%
  count(chr, chr_type) %>%
  mutate(chr = factor(chr, levels = c(1:22, "X")))

p1 <- ggplot(candidates_by_chr, aes(x = chr, y = n, fill = chr_type)) +
  geom_col() +
  scale_fill_manual(values = c("Autosome" = "gray50", "X chromosome" = "#E41A1C")) +
  labs(
    title = "Number of Selection Candidates by Chromosome",
    subtitle = "Moderate threshold (|Std iHS| >= 2.5)",
    x = "Chromosome",
    y = "Number of Candidates",
    fill = "Chromosome Type"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

ggsave(file.path(output_dir, "plots/candidates_by_chromosome.png"),
       p1, width = 12, height = 6, dpi = 300)
cat("    Saved: plots/candidates_by_chromosome.png\n")

# Plot 2: Candidate rate by chromosome
cat("  Plotting candidate rate by chromosome...\n")

candidate_rates <- candidates_by_chr %>%
  left_join(chr_stats %>% select(chr, n_variants), by = "chr") %>%
  mutate(rate = n / n_variants * 100)

p2 <- ggplot(candidate_rates, aes(x = chr, y = rate, fill = chr_type)) +
  geom_col() +
  scale_fill_manual(values = c("Autosome" = "gray50", "X chromosome" = "#E41A1C")) +
  labs(
    title = "Selection Candidate Rate by Chromosome",
    subtitle = "Moderate threshold (|Std iHS| >= 2.5)",
    x = "Chromosome",
    y = "Candidate Rate (%)",
    fill = "Chromosome Type"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

ggsave(file.path(output_dir, "plots/candidate_rate_by_chromosome.png"),
       p2, width = 12, height = 6, dpi = 300)
cat("    Saved: plots/candidate_rate_by_chromosome.png\n")

# Plot 3: Distribution of Std iHS in candidates (X vs autosomes)
cat("  Plotting Std iHS distribution in candidates...\n")

p3 <- ggplot(candidates_list$moderate, aes(x = `Std iHS`, fill = chr_type)) +
  geom_histogram(alpha = 0.6, bins = 50, position = "identity") +
  scale_fill_manual(values = c("Autosome" = "gray50", "X chromosome" = "#E41A1C")) +
  labs(
    title = "Distribution of Standardized iHS in Selection Candidates",
    subtitle = "Moderate threshold (|Std iHS| >= 2.5)",
    x = "Standardized iHS",
    y = "Count",
    fill = "Chromosome Type"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

ggsave(file.path(output_dir, "plots/candidates_std_ihs_distribution.png"),
       p3, width = 10, height = 6, dpi = 300)
cat("    Saved: plots/candidates_std_ihs_distribution.png\n\n")

cat("===========================================\n")
cat("Analysis Complete!\n")
cat("===========================================\n\n")

cat("Output files:\n")
cat("\nPer-chromosome candidates:\n")
cat("  - candidates/per_chromosome/chr*_candidates_*.tsv (69 files)\n")
cat("    Each file includes: Location, chr, pos, ref, alt, iHH_0, iHH_1,\n")
cat("                        iHS, Std_iHS, p_value, neg_log10_p\n")
cat("\nCombined files:\n")
cat("  - candidates/chromosome_statistics.tsv (includes empirical percentiles)\n")
cat("  - candidates/empirical_percentiles.tsv (99th, 99.5th, 99.9th percentiles)\n")
cat("  - candidates/ALL_candidates_liberal.tsv (|Std iHS| >= 2.0)\n")
cat("  - candidates/ALL_candidates_moderate.tsv (|Std iHS| >= 2.5)\n")
cat("  - candidates/ALL_candidates_stringent.tsv (|Std iHS| >= 3.0)\n")
cat("  - candidates/x_enrichment_test.tsv\n")
cat("\nPlots:\n")
cat("  - plots/candidates_by_chromosome.png\n")
cat("  - plots/candidate_rate_by_chromosome.png\n")
cat("  - plots/candidates_std_ihs_distribution.png\n\n")

cat("=== Threshold Interpretation (Literature-based) ===\n")
cat("Based on Voight et al. 2006, Paul et al. 2024, and others:\n\n")
cat("Fixed thresholds (assuming Std iHS ~ N(0,1)):\n")
cat("  |Std iHS| >= 2.0: ~2.3% of genome (liberal, Voight et al. 2006)\n")
cat("  |Std iHS| >= 2.5: ~0.6% of genome (moderate, Paul et al. 2024)\n")
cat("  |Std iHS| >= 3.0: ~0.13% of genome (stringent)\n\n")
cat("Equivalent p-value thresholds:\n")
cat("  |Std iHS| >= 2.0: p < 0.046, -log10(p) > 1.34\n")
cat("  |Std iHS| >= 2.5: p < 0.012, -log10(p) > 1.91\n")
cat("  |Std iHS| >= 3.0: p < 0.0027, -log10(p) > 2.57\n")
cat("  Paul et al. 2024 used: -log10(p) > 4 (p < 0.0001)\n\n")
cat("IMPORTANT: Check empirical percentiles in your data!\n")
cat("  If distributions deviate from normality, use empirical percentiles\n")
cat("  rather than fixed thresholds (Salazar-Tortosa et al. 2023)\n\n")

cat("Next steps:\n")
cat("  1. Review the normality check and empirical percentiles\n")
cat("  2. Check if X chromosome distribution differs from autosomes\n")
cat("  3. Decide on threshold: fixed (2.0/2.5/3.0) or empirical (top 1%/0.5%)\n")
cat("  4. Review plots to understand selection patterns\n")
cat("  5. Run 02_define_selection_regions.R to cluster candidates into peaks\n")
cat("      (Following literature: require multiple SNPs per region, not single outliers)\n")
cat("  6. Annotate peaks with genes\n")
cat("  7. Perform enrichment analysis for reproductive/immune genes\n\n")
