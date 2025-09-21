######################################################################################################################
####### Field-scale 16S rRNA dataset on rhizosphere microbial community under Rhodobacter capsulatus treatment ####### 
######################################################################################################################
# Date of creation: 2025-07-13

# ================================================================================================
# IMPORTANT:
# Before running this script, please make sure to execute "Integrated_preprocessing_rhizosphere_data.R".
# That script performs the essential preprocessing steps and generates the input files required here.
# ================================================================================================

# Fig 4a Relatice abundance 
# ================================================================================================
# This script generates stacked bar plots of relative abundance (RA) at the phylum level.
# The analysis covers multiple experimental sites and different sampling times (1st and 2nd).
# The workflow includes:
#   1. Importing phylum-level RA data from multiple experimental groups.
#   2. Preprocessing and merging the datasets.
#   3. Plotting stacked bar charts of microbial community composition 
#      for each sampling event (1st and 2nd).
# ================================================================================================

# ================================================================================================
# Working directory: set to your Kale project folder
# If you run this script from a subfolder, set the project root like below:
# setwd("./../Kale/")     # relative path example
# Or use an absolute path:
# setwd("/path/to/Kale")  # replace with your actual path
# Ensure the input files exist under Kale/, e.g. AE/output/phylum_RA.csv, CK1/output/...
# ================================================================================================

# Load necessary libraries
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggh4x)

## read genus_RA
AE = read.csv("AE/output/phylum_RA.csv")
CK1 = read.csv("CK1/output/phylum_RA.csv")
EK1 = read.csv("EK1/output/phylum_RA.csv")
EK2 = read.csv("EK2_rep/output/phylum_RA.csv")
EK3 = read.csv("EK3/output/phylum_RA.csv")
EK4 = read.csv("EK4/output/phylum_RA.csv")

## remove unnecessary columns for each dataset
sel_AE = AE %>% select(-c(1,2,6,9))
sel_CK1 = CK1 %>% select(-c(1,2,6,9))
sel_EK1 = EK1 %>% select(-c(1,2,6,9))
sel_EK2 = EK2%>% select(-c(1,2,6,9))
sel_EK3 = EK3 %>% select(-c(1,2,6,9))
sel_EK4 = EK4 %>% select(-c(1,2,6,9))

## merge all datasets into one
merge_data = rbind(sel_AE, sel_CK1, sel_EK1, sel_EK2, sel_EK3, sel_EK4)
head(merge_data)

#########################################################
##                  1st sampling                       ##
#########################################################
sampling_1st = merge_data %>% 
  filter(sampling == "1") %>% 
  mutate(new_col = case_when(
    Sample.site == "control"   ~ "Ctrl-S1",
    Sample.site == "Treatment" ~ "Trt-S1",
    TRUE             ~ NA_character_
  ))

## Color palette settings for plotting
getPalette = colorRampPalette(brewer.pal(11, "Spectral"))


## Phylum-level stacked barplot (relative abundance)
colourCount_P_1 = length(unique(sampling_1st$Phylum))
plot_1 = ggplot(sampling_1st, aes(x = Sample, y = Abundance, fill = Phylum)) +
    geom_bar(stat = "identity", position = "fill", width = 1) +
    facet_nested(~ site_ID + new_col,
                 scales = "free_x") +
    scale_x_discrete(labels = NULL, breaks = NULL) +   
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
    scale_fill_manual(values = getPalette(colourCount_P_1)) +
  theme(
    axis.title.x   = element_blank(),
    axis.text.x    = element_blank(),
    axis.ticks.x   = element_blank(),
    axis.title.y = element_blank(),
    strip.text.x.top    = element_text(face = "bold", size = 8, margin = margin(b = 8)),  
    strip.background = element_blank(),
    panel.spacing.x = unit(0,"line"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_blank(),
    legend.title = element_text(size = 8, face = "bold"),  
    legend.text  = element_text(size = 8),                 
    legend.key.size = unit(0.5, "cm")                      
  ) +
  guides(fill = guide_legend(ncol = 1)) 

plot_1

plot_1_nlegend = plot_1 + theme(legend.position = "none")


#########################################################
##                  2nd sampling                       ##
#########################################################
sampling_2nd = merge_data %>% 
  filter(sampling == "2") %>% 
  mutate(new_col = case_when(
    Sample.site == "control"   ~ "Ctrl-S2",
    Sample.site == "Treatment" ~ "Trt-S2",
    TRUE             ~ NA_character_
  ))

## Color palette settings for plotting
getPalette = colorRampPalette(brewer.pal(11, "Spectral")) # plots Color Palette Settings


## Phylum-level stacked barplot (relative abundance)
colourCount_P_2 = length(unique(sampling_2nd$Phylum))
plot_2 = ggplot(sampling_2nd, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill", width = 1) +
  facet_nested(~ site_ID + new_col,
               #nest_line = TRUE,
               scales = "free_x") +
  scale_x_discrete(labels = NULL, breaks = NULL) +   
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = getPalette(colourCount_P_2)) +
  theme(
    axis.title.x   = element_blank(),
    axis.text.x    = element_blank(),
    axis.ticks.x   = element_blank(),
    axis.title.y = element_blank(),
    strip.text.x.top    = element_text(face = "bold", size = 8, margin = margin(b = 8)), 
    strip.background = element_blank(),
    panel.spacing.x = unit(0,"line"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_blank(),
    legend.title = element_text(size = 8, face = "bold"),  
    legend.text  = element_text(size = 8),                 
    legend.key.size = unit(0.5, "cm")                      
  ) +
  guides(fill = guide_legend(ncol = 1)) 

plot_2

plot_2_nlegend = plot_2 + theme(legend.position = "none")
