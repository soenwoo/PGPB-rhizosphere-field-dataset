######################################################################################################################
####### Field-scale 16S rRNA dataset on rhizosphere microbial community under Rhodobacter capsulatus treatment ####### 
######################################################################################################################
# Date of creation: 2025-07-13


# ================================================================================================
# IMPORTANT:
# Before running this script, please make sure to execute "Integrated_preprocessing_rhizosphere_data.R".
# That script performs the essential preprocessing steps and generates the input files required here.
# ================================================================================================

# Fig 5 alpha diversity
# =======================================================================================================================
# This script produces per-site alpha diversity boxplots (Observed & Shannon) and overlays per-pair significance marks.
# It builds Ctr/Trt × S1/S2 group labels, checks normality per group, then auto-selects t-test or Wilcoxon per site.
# Boxplots with jittered points are shown; p-value stars (***, **, *, NS) are drawn between paired groups.
#
# Workflow:
#   1) Read alpha-diversity table (e.g., "output/alpha.csv") 
#   2) Create compact group labels: control/Treatment × sampling(1/2) → Ctr-S1, Trt-S1, Ctr-S2, Trt-S2, and fix order.
#   3) For each site_ID:
#        - Run Shapiro–Wilk tests per group to assess normality.
#        - If both groups ~ normal → use t.test; otherwise → use wilcox.test.
#        - Compute p-values and significance symbols; set y-positions for annotations.
#   4) Draw faceted boxplots (one facet per site) with jittered points and significance stars.
# =======================================================================================================================

# =======================================================================================================================
# Directory Setup for Farm-Specific Analysis
# This script requires farm-specific analysis. 
# For each farm, make sure to set the working directory correctly.
# For example, if analyzing data for "EK1", the working directory should be set as follows:
# setwd("./../Kale/EK1/") 
# Change the directory based on the farm you are analyzing.
# Ensure that the relevant data files are in the corresponding farm directory.
# =======================================================================================================================

# Load necessary libraries
library(ggplot2)  # boxplot, jitter, facet_wrap
library(dplyr)    # mutate, group_by, group_modify, filter
library(tibble)   # tibble(), column helpers
library(ggpubr)   # stat_pvalue_manual 
library(ggsci)    # scale_fill_lancet 
library(purrr)    # map_dfr 

# Load alpha diversity data 
alpha_data = read.csv("output/alpha.csv")

# Build compact group labels: control/Treatment × sampling(1/2) → Ctr-S1/Trt-S1/Ctr-S2/Trt-S2
names_alpha = alpha_data %>%
  mutate(control = ifelse(Sample.site == "control" & sampling == "1",
                          "Ctr-S1",
                          ifelse(Sample.site == "control" & sampling == "2",
                                 "Ctr-S2",
                                 ifelse(Sample.site == "Treatment" & sampling == "1",
                                        "Trt-S1",
                                        "Trt-S2"))))



# Fix plotting order on x-axis
newLabel = c("Ctr-S1", "Trt-S1", "Ctr-S2", "Trt-S2")
names_alpha$control = factor(names_alpha$control, levels = newLabel)
head(names_alpha)


# Per-site stats for Observed:
#  - Shapiro–Wilk normality test per group
#  - If both groups ~ normal → t.test; else → wilcox.test
comp_list = list(c("Ctr-S1","Trt-S1"), c("Ctr-S2","Trt-S2"))
stat_obs = names_alpha %>%
  group_by(site_ID) %>%
  group_modify(~ {
    map_dfr(comp_list, function(cp) {
      sub = .x %>% filter(control %in% cp)
      p1  = shapiro.test(sub$Observed[sub$control==cp[1]])$p.value
      p2  = shapiro.test(sub$Observed[sub$control==cp[2]])$p.value
      test_name = if (p1 > 0.05 & p2 > 0.05) "t.test" else "wilcox.test"
      res = if (test_name == "t.test") {
        t.test(Observed ~ control, data=sub)
      } else {
        wilcox.test(Observed ~ control, data=sub)
      }
      tibble(
        group1     = cp[1],
        group2     = cp[2],
        p.value    = res$p.value,
        p.signif   = symnum(res$p.value, corr = FALSE,
                            cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                            symbols   = c("***","**","*","NS")),
        y.position = max(sub$Observed, na.rm=TRUE) * 1.05,
        test = test_name
      )
    })
  }) %>% ungroup()
stat_obs

# Observed richness — faceted boxplot with jitter and significance stars
library(ggsci)
ob_plot_dynamic = ggplot(names_alpha, aes(x = control, y = Observed, fill = control)) +
  geom_boxplot(width = 0.4, outlier.shape = NA) +
  geom_jitter(size = 1, position = position_jitter(width = 0.1)) +
  facet_wrap(~site_ID, scales = "free_y") +
  stat_pvalue_manual(stat_obs,
                     inherit.aes = FALSE,
                     label      = "p.signif",
                     xmin       = "group1",
                     xmax       = "group2",
                     y.position = "y.position",
                     tip.length = 0.02) +
  scale_fill_lancet() +
  theme_bw(base_size = 16) +
  theme(
    strip.text       = element_text(size = 14, face = "bold"),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.title.x     = element_blank(),
    axis.title.y = element_blank()
  ) 
ob_plot_dynamic
no_legend_ob = ob_plot_dynamic +  theme(legend.position = "none")

# Repeat stats for Shannon index (same testing logic)
stat_shannon = names_alpha %>%
  group_by(site_ID) %>%
  group_modify(~ {
    map_dfr(comp_list, function(cp) {
      sub = .x %>% filter(control %in% cp)
      p1  = shapiro.test(sub$Shannon[sub$control==cp[1]])$p.value
      p2  = shapiro.test(sub$Shannon[sub$control==cp[2]])$p.value
      test_name = if (p1 > 0.05 & p2 > 0.05) "t.test" else "wilcox.test"
      res = if (test_name == "t.test") {
        t.test(Shannon ~ control, data=sub)
      } else {
        wilcox.test(Shannon ~ control, data=sub)
      }
      tibble(
        group1     = cp[1],
        group2     = cp[2],
        p.value    = res$p.value,
        p.signif   = symnum(res$p.value, corr = FALSE,
                            cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                            symbols   = c("***","**","*","NS")),
        y.position = max(sub$Shannon, na.rm=TRUE) * 1.05,
        test = test_name
      )
    })
  }) %>% ungroup()
stat_shannon

# Shannon — faceted boxplot with jitter and significance stars
sh_plot_dynamic = ggplot(names_alpha, aes(x = control, y = Shannon, fill = control)) +
  geom_boxplot(width = 0.4, outlier.shape = NA) +
  geom_jitter(size = 1, position = position_jitter(width = 0.1)) +
  facet_wrap(~site_ID, scales = "free_y") +
  stat_pvalue_manual(stat_shannon,
                     inherit.aes = FALSE,
                     label      = "p.signif",
                     xmin       = "group1",
                     xmax       = "group2",
                     y.position = "y.position",
                     tip.length = 0.02) +
  scale_fill_lancet() +
  theme_bw(base_size = 16) +
  theme(
    strip.text       = element_text(size = 14, face = "bold"),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.title.x     = element_blank(),
    axis.title.y = element_blank()
  ) 

sh_plot_dynamic
no_legend_sh = sh_plot_dynamic + theme(legend.position = "none")
