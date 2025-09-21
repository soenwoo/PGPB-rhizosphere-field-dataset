######################################################################################################################
####### Field-scale 16S rRNA dataset on rhizosphere microbial community under Rhodobacter capsulatus treatment ####### 
######################################################################################################################
# Date of creation: 2025-07-13

# ================================================================================================
# IMPORTANT:
# Before running this script, please make sure to execute "Integrated_preprocessing_rhizosphere_data.R".
# That script performs the essential preprocessing steps and generates the input files required here.
# ================================================================================================

# Fig 6b beta diversity, yield association heatmap
# =================================================================================================
# Title: Association between β-diversity differences and yield across farms
# Purpose:
#   Summarize whether farms showing significant β-diversity differences (PGPB vs. control)
#   also tend to show increased or decreased crop yield at two time points (1st, 2nd).
#   The heatmap tiles display the COUNT of farms in each combination:
#     - β-diversity: Difference vs. NO_difference (by p ≤ 0.05)
#     - Yield: Increased vs. Decreased (by yield_change > 0)
#   Note: β-diversity p-values (beta_p_value) were computed via PERMANOVA
#         (vegan::adonis2) on Bray–Curtis distances (e.g., from CSS-/library-normalized counts).
#   This mirrors the main-figure legend: left panel = 1st time point, right panel = 2nd time point.
#
# Workflow:
#   1) Read per–time-point CSVs: beta_yield_1st.csv, beta_yield_2nd.csv.
#   2) Derive two flags per farm:
#        - yield_direction = Increased / Decreased (from yield_change)
#        - beta_diff_flag  = Difference / NO_difference (from beta_p_value threshold 0.05)
#   3) Cross-tabulate counts for (beta_diff_flag × yield_direction) within each time point.
#   4) Bind counts from both time points; plot a heatmap with tile labels (counts), faceted by Time point.
# =================================================================================================

# Set the working directory to the folder that contains the data files.
# Make sure "beta_yield_1st.csv" and "beta_yield_2nd.csv" are inside the "Beta_yield" folder.
setwd(".././Beta_yield/")

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)


# Read per–time-point CSVs (1st, 2nd)
beta_yield_1st = read.csv("beta_yield_1st.csv")
beta_yield_2nd = read.csv("beta_yield_2nd.csv")


# Derive flags per farm (yield direction & β-diversity difference) 
annot_1st <- beta_yield_1st %>%
  mutate(
    Yield_dir = if_else(yield_change  > 0, "Increased", "Decreased"),
    beta_diff_flag = if_else(beta_p_value  <= 0.05, "Difference", "NO_difference")
  )

annot_2nd <- beta_yield_2nd %>%
  mutate(
    Yield_dir = if_else(yield_change  > 0, "Increased", "Decreased"),
    beta_diff_flag = if_else(beta_p_value   <= 0.05, "Difference", "NO_difference")
  )

# Cross-tabulate counts: (β-diversity change × yield direction) per time point 
counts_1st <- annot_1st %>%
  dplyr::count(beta_diff_flag, Yield_dir) %>%
  complete(beta_diff_flag, Yield_dir, fill = list(n = 0)) %>%
  dplyr::mutate(Sampling = "1st") %>%
  dplyr::rename(Beta = beta_diff_flag)

counts_1st

counts_2nd <- annot_2nd %>%
  dplyr::count(beta_diff_flag, Yield_dir) %>%
  complete(beta_diff_flag, Yield_dir, fill = list(n = 0)) %>%
  dplyr::mutate(Sampling = "2nd") %>%
  dplyr::rename(Beta = beta_diff_flag)

counts_2nd

# Combine time points
counts_all <- bind_rows(counts_1st, counts_2nd)

### Plot heatmap with count labels
p = ggplot(counts_all, aes(x = Yield_dir, y = Beta, fill = n)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n), size = 6) +
  scale_fill_gradient(low = "white", high = "blue") +
  facet_wrap(~ Sampling, nrow = 1) +
  labs(
    x = "Yield direction",
    y = expression(beta*"-diversity change"),
    fill = "Count"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),   # x축 라벨 45도
    axis.text.y   = element_text(angle = 90, vjust = 0.5, hjust = 0.5,
                                 margin = margin(r = 6))
  )+ 
  theme(legend.position = "none") 

p
