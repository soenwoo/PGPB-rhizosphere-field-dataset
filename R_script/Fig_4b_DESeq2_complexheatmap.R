######################################################################################################################
####### Field-scale 16S rRNA dataset on rhizosphere microbial community under Rhodobacter capsulatus treatment ####### 
######################################################################################################################
# Date of creation: 2025-07-13

# ================================================================================================
# IMPORTANT:
# Before running this script, please make sure to execute "Integrated_preprocessing_rhizosphere_data.R".
# That script performs the essential preprocessing steps and generates the input files required here.
# ================================================================================================

# Fig 4b DESeq2 heatmap
# =======================================================================================================================
# This script builds per-site heatmaps of DESeq2-identified differentially abundant genera for two sampling times
# (1st, 2nd). Counts are DESeq2-normalized and log-scaled upstream; Phylum is used as a bottom annotation to keep
# taxonomy visible across panels. Uncultured/ambiguous taxa are excluded in prior preprocessing.

# Workflow:
#   1) Read per-site Genus-level DESeq2 result tables (1st/2nd) and Phylum annotation tables.
#   2) Harmonize sample headers to "control"/"Treatment" (N → control; T → Treatment) and set Genus as rownames.
#   3) Build a unified Phylum color palette from the union of Phyla (keeps colors consistent across sites/times).
#   4) Z-score by Genus (transpose → center/scale) to create heatmap matrices.
#   5) Create site-level heatmaps with Phylum bottom annotations (and optional control/Treatment row annotations).
# =======================================================================================================================

# Working directory:
#   # setwd("./../Kale/")   

# Load necessary libraries
library(ComplexHeatmap)
library(dplyr)
library(tibble)
library(RColorBrewer)
library(circlize)

######################################
###          1st sampling          ###
######################################

# Read Genus-level DESeq2 result tables (one per site) for the 1st sampling
otu_AE_1st = read.csv("AE/output/1st_sampling_Genus_DESeq2_results.csv")
otu_CK1_1st = read.csv("CK1/output/1st_sampling_Genus_DESeq2_results.csv")
otu_EK1_1st = read.csv("EK1/output/1st_sampling_Genus_DESeq2_results.csv")
otu_EK2_1st = read.csv("EK2/output/1st_sampling_Genus_DESeq2_results.csv")
otu_EK3_1st = read.csv("EK3/output/1st_sampling_Genus_DESeq2_results.csv")
otu_EK4_1st = read.csv("EK4/output/1st_sampling_Genus_DESeq2_results.csv")

# Helper to standardize sample headers → "control"/"Treatment" based on N/T patterns
rename_control_treatment <- function(df,
                                     control_label   = "control",
                                     treatment_label = "Treatment",
                                     pattern_control = "N",
                                     pattern_treat   = "T"){

  orig <- colnames(df)
  new <- ifelse(
    grepl(pattern_control, orig),
    sub(paste0("^.*?", pattern_control),
        control_label,
        orig),
    ifelse(
      grepl(pattern_treat, orig),
      sub(paste0("^.*?", pattern_treat),
          treatment_label,
          orig),
      orig
    )
  )
  colnames(df) <- new
  return(df)
}

# Apply header renaming; move Genus column "X" to rownames to form a matrix-like object
AE_1st = rename_control_treatment(otu_AE_1st) %>%
  column_to_rownames(var = "X")
CK1_1st = rename_control_treatment(otu_CK1_1st) %>%
  column_to_rownames(var = "X")
EK1_1st = rename_control_treatment(otu_EK1_1st) %>%
  column_to_rownames(var = "X")
EK2_1st = rename_control_treatment(otu_EK2_1st) %>%
  column_to_rownames(var = "X")
EK3_1st = rename_control_treatment(otu_EK3_1st) %>%
  column_to_rownames(var = "X")
EK4_1st = rename_control_treatment(otu_EK4_1st) %>%
  column_to_rownames(var = "X")

### Read Phylum annotation tables (1st & 2nd) to build consistent color mapping across panels
# phylum data read (1st)
phy_AE_1st = read.csv("AE/output/1st_phylum_annotation_DESeq2.csv")
phy_CK1_1st = read.csv("CK1/output/1st_phylum_annotation_DESeq2.csv")
phy_EK1_1st = read.csv("EK1/output/1st_phylum_annotation_DESeq2.csv")
phy_EK2_1st = read.csv("EK2/output/1st_phylum_annotation_DESeq2.csv")
phy_EK3_1st = read.csv("EK3/output/1st_phylum_annotation_DESeq2.csv")
phy_EK4_1st = read.csv("EK4/output/1st_phylum_annotation_DESeq2.csv")

# phylum data read (1st)
phy_AE_2nd = read.csv("AE/output/2nd_phylum_annotation_DESeq2.csv")
phy_CK1_2nd = read.csv("CK1/output/2nd_phylum_annotation_DESeq2.csv")
phy_EK1_2nd = read.csv("EK1/output/2nd_phylum_annotation_DESeq2.csv")
phy_EK2_2nd = read.csv("EK2/output/2nd_phylum_annotation_DESeq2.csv")
phy_EK3_2nd = read.csv("EK3/output/2nd_phylum_annotation_DESeq2.csv")
phy_EK4_2nd = read.csv("EK4/output/2nd_phylum_annotation_DESeq2.csv")



# Build a global Phylum palette from the union of all Phyla (1st + 2nd); keeps colors consistent
combined_unique <- bind_rows(phy_AE_1st, phy_CK1_1st, phy_EK1_1st, phy_EK2_1st,
                             phy_EK3_1st, phy_EK4_1st,
                             phy_AE_2nd, phy_CK1_2nd,phy_EK1_2nd, phy_EK2_2nd, 
                             phy_EK3_2nd, phy_EK4_2nd) %>% distinct()

combined_unique_1st =  bind_rows(phy_AE_1st, phy_CK1_1st, phy_EK1_1st, 
                                 phy_EK2_1st, phy_EK3_1st, phy_EK4_1st) %>% distinct()


phylum_levels_total <- sort(unique(combined_unique$Phylum))
phylum_levels_1st = sort(unique(combined_unique_1st$Phylum))

phylum_colors <- setNames(
  colorRampPalette(brewer.pal(8, "Set3"))(length(phylum_levels_total)),
  phylum_levels_total
)


### AE_heatmap
## AE (1st) — scale by Genus and prepare annotations
# Z-score by Genus: transpose → scale rows (center/scale) → matrix for heatmap
AE_col_scaled_1 = scale( t(as.matrix(AE_1st)),
                         center = TRUE,
                         scale  = TRUE )
AE_col_scaled_1


# Optional row annotation for sample groups (control vs Treatment) — adjust counts to your design
AE_row_ant = data.frame(
  x = colnames(AE_1st),
  category = c(rep("control",8), rep("Treatment", 8))
)

AE_row_ha = rowAnnotation(
  df = AE_row_ant %>%
    column_to_rownames(var = "x"),
  col = list(
    category = c("control" = "lightgray", "Treatment" = "orange")
  ),
  show_annotation_name = F,
  show_legend = F,
  annotation_name_side = "bottom",
  annotation_name_gp = gpar(fontsize = 10),
  width = unit(3, "mm"),
  simple_anno_size_adjust = TRUE 
)


# Bottom annotation: Phylum per Genus (factor-levels fixed for consistent colors)
AE_col_ant_1 = phy_AE_1st %>%
  column_to_rownames(var = "X")

AE_col_ant_1$Phylum <- factor(AE_col_ant_1$Phylum, levels = phylum_levels_1st)

phylum_colors

AE_col_ha_1 = HeatmapAnnotation(
  df = AE_col_ant_1,
  col    = list(Phylum = phylum_colors),
  show_annotation_name = FALSE,
  show_legend = TRUE,
  annotation_name_side = "right",
  annotation_legend_param = list(
    Phylum =list(
      title = "Phylum",
      at = phylum_levels_1st,
      labels = phylum_levels_1st
    )
  )
)

# Build the AE (1st) heatmap 
AE_ht_1 = Heatmap(
  AE_col_scaled_1,
  name            = "AE",
  column_title = "AE", 
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  col = colorRamp2(c(-4, 0, 4), c("#2471A3", "white", "#C0392B")),
  bottom_annotation  = AE_col_ha_1,
  left_annotation= AE_row_ha,
  show_row_names  = FALSE,
  show_column_names = FALSE,
  show_heatmap_legend = TRUE,
  cluster_rows    = FALSE,
  cluster_columns = TRUE,
  heatmap_legend_param = list(
    title = NULL,
    legend_direction = "vertical",  
    legend_height = unit(6, "cm"), 
    grid_width = unit(4, "mm")
  )
)
AE_ht_1

##  Repeat the same pattern for CK1, EK1, EK2, EK3, EK4 (1st sampling)
## CK1_heatmap
# Z-score by Genus: transpose → scale rows (center/scale) 
CK1_col_scaled_1 = scale( t(as.matrix(CK1_1st)),
                        center = TRUE, 
                        scale  = TRUE )

# Bottom annotation: Phylum per Genus (factor-levels fixed for consistent colors)
CK1_col_ant_1 = phy_CK1_1st %>%
  column_to_rownames(var = "X")

CK1_col_ant_1$Phylum <- factor(CK1_col_ant_1$Phylum, levels = phylum_levels_1st)

phylum_colors

CK1_col_ha_1 = HeatmapAnnotation(
  df = CK1_col_ant_1,
  col    = list(Phylum = phylum_colors),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  annotation_name_side = "right",
  annotation_legend_param = list(
    Phylum =list(
      title = "Phylum",
      at = phylum_levels_1st,
      labels = phylum_levels_1st
    )
  )
)

#Build the CK1 (1st) heatmap
CK1_ht_1 = Heatmap(
  CK1_col_scaled_1,
  name            = "CK1",
  column_title = "CK1",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  col = colorRamp2(c(-4, 0, 4), c("#2471A3", "white", "#C0392B")),
  bottom_annotation= CK1_col_ha_1,
  show_row_names  = FALSE,
  show_column_names = FALSE,
  show_heatmap_legend = FALSE,
  cluster_rows    = FALSE,
  cluster_columns = TRUE,
  heatmap_legend_param = list(
    title = NULL,
    legend_direction = "vertical",  
    legend_height = unit(6, "cm"), 
    grid_width = unit(4, "mm")
  )
)
CK1_ht_1

##################################################################
##EK1_heatmap
# Z-score by Genus: transpose → scale rows (center/scale) 
EK1_col_scaled_1 = scale( t(as.matrix(EK1_1st)),
                        center = TRUE,   
                        scale  = TRUE )

# Bottom annotation: Phylum per Genus (factor-levels fixed for consistent colors)
EK1_col_ant_1 = phy_EK1_1st %>%
  column_to_rownames(var = "X")

EK1_col_ant_1$Phylum <- factor(EK1_col_ant_1$Phylum, levels = phylum_levels_1st)

phylum_colors

EK1_col_ha_1 = HeatmapAnnotation(
  df = EK1_col_ant_1,
  col    = list(Phylum = phylum_colors),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  annotation_name_side = "right",
  annotation_legend_param = list(
    Phylum =list(
      title = "Phylum",
      at = phylum_levels_1st,
      labels = phylum_levels_1st
    )
  )
)

#Build the EK1 (1st) heatmap
EK1_ht_1 = Heatmap(
  EK1_col_scaled_1,
  name            = "EK1",
  column_title = "EK1",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  col = colorRamp2(c(-4, 0, 4), c("#2471A3", "white", "#C0392B")),
  bottom_annotation= EK1_col_ha_1,
  show_row_names  = FALSE,
  show_column_names = FALSE,
  show_heatmap_legend = FALSE,
  cluster_rows    = FALSE,
  cluster_columns = TRUE,
  heatmap_legend_param = list(
    title = NULL,
    legend_direction = "vertical",
    legend_height = unit(6, "cm"),
    grid_width = unit(4, "mm")
  )
)
EK1_ht_1

#####################################################################
##EK2_heatmap
# Z-score by Genus: transpose → scale rows (center/scale) 
EK2_col_scaled_1 = scale( t(as.matrix(EK2_1st)),
                        center = TRUE,
                        scale  = TRUE )

# Bottom annotation: Phylum per Genus (factor-levels fixed for consistent colors)
EK2_col_ant_1 = phy_EK2_1st %>%
  column_to_rownames(var = "X")

EK2_col_ant_1$Phylum <- factor(EK2_col_ant_1$Phylum, levels = phylum_levels_1st)

phylum_colors

EK2_col_ha_1 = HeatmapAnnotation(
  df = EK2_col_ant_1,
  col    = list(Phylum = phylum_colors),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  annotation_name_side = "right",
  annotation_legend_param = list(
    Phylum =list(
      title = "Phylum",
      at = phylum_levels_1st,
      labels = phylum_levels_1st
    )
  )
)

#Build the EK2 (1st) heatmap
EK2_ht_1 = Heatmap(
  EK2_col_scaled_1,
  name            = "EK2",
  column_title = "EK2",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  col = colorRamp2(c(-4, 0, 4), c("#2471A3", "white", "#C0392B")),
  bottom_annotation= EK2_col_ha_1,
  show_row_names  = FALSE,
  show_column_names = FALSE,
  show_heatmap_legend = FALSE,
  cluster_rows    = FALSE,
  cluster_columns = TRUE,
  heatmap_legend_param = list(
    title = NULL,
    legend_direction = "vertical",
    legend_height = unit(6, "cm"),
    grid_width = unit(4, "mm")
  )
)
EK2_ht


##EK3_heatmap
# Z-score by Genus: transpose → scale rows (center/scale) 
EK3_col_scaled_1 = scale( t(as.matrix(EK3_1st)),
                        center = TRUE, 
                        scale  = TRUE )

# Bottom annotation: Phylum per Genus (factor-levels fixed for consistent colors)
EK3_col_ant_1 = phy_EK3_1st %>%
  column_to_rownames(var = "X")

EK3_col_ant_1$Phylum <- factor(EK3_col_ant_1$Phylum, levels = phylum_levels_1st)

phylum_colors

EK3_col_ha_1 = HeatmapAnnotation(
  df = EK3_col_ant_1,
  col    = list(Phylum = phylum_colors),
  show_annotation_name = FALSE,
  show_legend = TRUE,
  annotation_name_side = "right",
  annotation_legend_param = list(
    Phylum =list(
      title = "Phylum",
      at = phylum_levels_1st,
      labels = phylum_levels_1st
    )
  )
)


#Build the EK3 (1st) heatmap
EK3_ht_1 = Heatmap(
  EK3_col_scaled_1,
  name            = "EK3",
  column_title = "EK3", 
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  col = colorRamp2(c(-4, 0, 4), c("#2471A3", "white", "#C0392B")),
  bottom_annotation= EK3_col_ha_1,
  show_row_names  = FALSE,
  show_column_names = FALSE,
  show_heatmap_legend = FALSE,
  cluster_rows    = FALSE,
  cluster_columns = TRUE,
  heatmap_legend_param = list(
    title = NULL,
    legend_direction = "vertical",
    legend_height = unit(6, "cm"),
    grid_width = unit(4, "mm")
  )
)
EK3_ht_1


##EK4_heatmap
# Z-score by Genus: transpose → scale rows (center/scale) 
EK4_col_scaled_1 = scale( t(as.matrix(EK4_1st)),
                        center = TRUE,
                        scale  = TRUE )

# Bottom annotation: Phylum per Genus (factor-levels fixed for consistent colors)
EK4_col_ant_1 = phy_EK4_1st %>%
  column_to_rownames(var = "X")

EK4_col_ant_1$Phylum <- factor(EK4_col_ant_1$Phylum, levels = phylum_levels_1st)

phylum_colors

EK4_col_ha_1 = HeatmapAnnotation(
  df = EK4_col_ant_1,
  col    = list(Phylum = phylum_colors),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  annotation_name_side = "right",
  annotation_legend_param = list(
    Phylum =list(
      title = "Phylum",
      at = phylum_levels_1st,
      labels = phylum_levels_1st
    )
  )
)


#Build the EK4 (1st) heatmap
EK4_ht_1 = Heatmap(
  EK4_col_scaled_1,
  name            = "EK4",
  column_title = "EK4",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  col = colorRamp2(c(-4, 0, 4), c("#2471A3", "white", "#C0392B")),
  bottom_annotation= EK4_col_ha_1,
  show_row_names  = FALSE,
  show_column_names = FALSE,
  show_heatmap_legend = FALSE,
  cluster_rows    = FALSE,
  cluster_columns = TRUE,
  heatmap_legend_param = list(
    title = NULL,
    legend_direction = "vertical",
    legend_height = unit(6, "cm"),
    grid_width = unit(4, "mm")
  )
)
EK4_ht_1

# Combine site panels for the 1st sampling and draw
ht_1_list <-AE_ht_1 + CK1_ht_1 + EK1_ht_1 + EK2_ht_1 + EK3_ht_1 + EK4_ht_1
draw(ht_1_list)

######################################
###          2nd sampling          ###
######################################
# Repeat the pipeline for the 2nd sampling (read → rename → scale → annotate → heatmaps)

### Relative abundance data read
otu_AE_2nd = read.csv("AE/output/2nd_sampling_Genus_DESeq2_results.csv")
otu_CK1_2nd = read.csv("CK1/output/2nd_sampling_Genus_DESeq2_results.csv")
otu_EK1_2nd = read.csv("EK1/output/2nd_sampling_Genus_DESeq2_results.csv")
otu_EK2_2nd = read.csv("EK2/output/2nd_sampling_Genus_DESeq2_results.csv")
otu_EK3_2nd = read.csv("EK3/output/2nd_sampling_Genus_DESeq2_results.csv")
otu_EK4_2nd = read.csv("EK4/output/2nd_sampling_Genus_DESeq2_results.csv")

# Apply header renaming; move Genus column "X" to rownames to form a matrix-like object
AE_2nd = rename_control_treatment(otu_AE_2nd) %>%
  column_to_rownames(var = "X")
CK1_2nd = rename_control_treatment(otu_CK1_2nd) %>%
  column_to_rownames(var = "X")
EK1_2nd = rename_control_treatment(otu_EK1_2nd) %>%
  column_to_rownames(var = "X")
EK2_2nd = rename_control_treatment(otu_EK2_2nd) %>%
  column_to_rownames(var = "X")
EK3_2nd = rename_control_treatment(otu_EK3_2nd) %>%
  column_to_rownames(var = "X")
EK4_2nd = rename_control_treatment(otu_EK4_2nd) %>%
  column_to_rownames(var = "X")


# Build a global Phylum palette from the union of all Phyla (1st + 2nd); keeps colors consistent
combined_unique <- bind_rows(phy_AE_1st, phy_CK1_1st, phy_EK1_1st, phy_EK2_1st,
                             phy_EK3_1st, phy_EK4_1st,
                             phy_AE_2nd, phy_CK1_2nd,phy_EK1_2nd, phy_EK2_2nd, 
                             phy_EK3_2nd, phy_EK4_2nd) %>% 
  distinct()        

combined_unique_2nd =  bind_rows(phy_AE_2nd, phy_CK1_2nd, phy_EK1_2nd, 
                                 phy_EK2_2nd, phy_EK3_2nd, phy_EK4_2nd) %>% 
  distinct()        

phylum_levels_total <- sort(unique(combined_unique$Phylum))
phylum_levels_2nd = sort(unique(combined_unique_2nd$Phylum))

phylum_colors <- setNames(
  colorRampPalette(brewer.pal(8, "Set3"))(length(phylum_levels_total)),
  phylum_levels_total
)

### sampling replication 10 (n=10)
##AE_heatmap
##col scale(AE) 
AE_col_scaled_2 = scale( t(as.matrix(AE_2nd)),
                         center = TRUE,
                         scale  = TRUE )
AE_col_scaled_2


##row_annotataion (AE)
AE_row_ant_2 = data.frame(
  x = colnames(AE_2nd),
  category = c(rep("control",10), rep("Treatment", 10))
)

AE_row_ha_2 = rowAnnotation(
  df = AE_row_ant_2 %>%
    column_to_rownames(var = "x"),
  col = list(
    category = c("control" = "lightgray", "Treatment" = "orange")
  ),
  show_annotation_name = F,
  show_legend = F,
  annotation_name_side = "bottom",
  annotation_name_gp = gpar(fontsize = 10),
  width = unit(3, "mm"),
  simple_anno_size_adjust = TRUE
)


##col_annotation (AE)
AE_col_ant_2 = phy_AE_2nd %>%
  column_to_rownames(var = "X")

AE_col_ant_2$Phylum <- factor(AE_col_ant_2$Phylum, levels = phylum_levels_2nd)

phylum_colors

AE_col_ha_2 = HeatmapAnnotation(
  df = AE_col_ant_2,
  col    = list(Phylum = phylum_colors),
  show_annotation_name = FALSE,
  show_legend = TRUE,
  annotation_name_side = "right",
  annotation_legend_param = list(
    Phylum =list(
      title = "Phylum",
      at = phylum_levels_2nd,
      labels = phylum_levels_2nd
    )
  )
)

#plot
AE_ht_2 = Heatmap(
  AE_col_scaled_2,
  name            = "AE",
  column_title = "AE",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  col = colorRamp2(c(-4, 0, 4), c("#2471A3", "white", "#C0392B")),
  bottom_annotation  = AE_col_ha_2,
  left_annotation= AE_row_ha_2,
  show_row_names  = FALSE,
  show_column_names = FALSE,
  show_heatmap_legend = TRUE,
  cluster_rows    = FALSE,
  cluster_columns = TRUE,
  heatmap_legend_param = list(
    title = NULL,
    legend_direction = "vertical",
    legend_height = unit(6, "cm"),
    grid_width = unit(4, "mm")
  )
)
AE_ht_2

#####################################################################
##CK1_heatmap (n = 10)
##col scale(CK1)
CK1_col_scaled_2 = scale( t(as.matrix(CK1_2nd)),
                          center = TRUE,
                          scale  = TRUE )

##col_annotation (CK1)
CK1_col_ant_2 = phy_CK1_2nd %>%
  column_to_rownames(var = "X")

CK1_col_ant_2$Phylum <- factor(CK1_col_ant_2$Phylum, levels = phylum_levels_2nd)

phylum_colors

CK1_col_ha_2 = HeatmapAnnotation(
  df = CK1_col_ant_2,
  col    = list(Phylum = phylum_colors),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  annotation_name_side = "right",
  annotation_legend_param = list(
    Phylum =list(
      title = "Phylum",
      at = phylum_levels_2nd,
      labels = phylum_levels_2nd
    )
  )
)

#plot
CK1_ht_2 = Heatmap(
  CK1_col_scaled_2,
  name            = "CK1",
  column_title = "CK1",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  col = colorRamp2(c(-4, 0, 4), c("#2471A3", "white", "#C0392B")),
  bottom_annotation= CK1_col_ha_2,
  show_row_names  = FALSE,
  show_column_names = FALSE,
  show_heatmap_legend = FALSE,
  cluster_rows    = FALSE,
  cluster_columns = TRUE,
  heatmap_legend_param = list(
    title = NULL,
    legend_direction = "vertical",
    legend_height = unit(6, "cm"),
    grid_width = unit(4, "mm")
  )
)
CK1_ht_2

##################################################################
##EK4_heatmap (n=10)
##col scale(EK4)
EK4_col_scaled_2 = scale( t(as.matrix(EK4_2nd)),
                          center = TRUE,
                          scale  = TRUE )

##col_annotation (EK4)
EK4_col_ant_2 = phy_EK4_2nd %>%
  column_to_rownames(var = "X")

EK4_col_ant_2$Phylum <- factor(EK4_col_ant_2$Phylum, levels = phylum_levels_2nd)

phylum_colors

EK4_col_ha_2 = HeatmapAnnotation(
  df = EK4_col_ant_2,
  col    = list(Phylum = phylum_colors),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  annotation_name_side = "right",
  annotation_legend_param = list(
    Phylum =list(
      title = "Phylum",
      at = phylum_levels_2nd,
      labels = phylum_levels_2nd
    )
  )
)


#plot
EK4_ht_2 = Heatmap(
  EK4_col_scaled_2,
  name            = "EK4",
  column_title = "EK4",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  col = colorRamp2(c(-4, 0, 4), c("#2471A3", "white", "#C0392B")),
  bottom_annotation= EK4_col_ha_2,
  show_row_names  = FALSE,
  show_column_names = FALSE,
  show_heatmap_legend = FALSE,
  cluster_rows    = FALSE,
  cluster_columns = TRUE,
  heatmap_legend_param = list(
    title = NULL,
    legend_direction = "vertical",
    legend_height = unit(6, "cm"),
    grid_width = unit(4, "mm")
  )
)
EK4_ht_2

##################################################################
##EK1_heatmap (n=8)
##col scale(EK1)
EK1_col_scaled_2 = scale( t(as.matrix(EK1_2nd)),
                          center = TRUE,
                          scale  = TRUE )

##row_annotataion (EK1)
EK1_row_ant_2 = data.frame(
  x = colnames(EK1_2nd),
  category = c(rep("control",8), rep("Treatment", 8))
)

EK1_row_ha_2 = rowAnnotation(
  df = EK1_row_ant_2 %>%
    column_to_rownames(var = "x"),
  col = list(
    category = c("control" = "lightgray", "Treatment" = "orange")
  ),
  show_annotation_name = F,
  show_legend = F,
  annotation_name_side = "bottom",
  annotation_name_gp = gpar(fontsize = 10),
  width = unit(3, "mm"),
  simple_anno_size_adjust = TRUE 
)

##col_annotation (EK1)
EK1_col_ant_2 = phy_EK1_2nd %>%
  column_to_rownames(var = "X")

EK1_col_ant_2$Phylum <- factor(EK1_col_ant_2$Phylum, levels = phylum_levels_2nd)

phylum_colors

EK1_col_ha_2 = HeatmapAnnotation(
  df = EK1_col_ant_2,
  col    = list(Phylum = phylum_colors),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  annotation_name_side = "right",
  annotation_legend_param = list(
    Phylum =list(
      title = "Phylum",
      at = phylum_levels_2nd,
      labels = phylum_levels_2nd
    )
  )
)

#plot
EK1_ht_2 = Heatmap(
  EK1_col_scaled_2,
  name            = "EK1",
  column_title = "EK1",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  col = colorRamp2(c(-4, 0, 4), c("#2471A3", "white", "#C0392B")),
  bottom_annotation= EK1_col_ha_2,
  left_annotation= EK1_row_ha_2,
  show_row_names  = FALSE,
  show_column_names = FALSE,
  show_heatmap_legend = TRUE,
  cluster_rows    = FALSE,
  cluster_columns = TRUE,
  heatmap_legend_param = list(
    title = NULL,
    legend_direction = "vertical", 
    legend_height = unit(6, "cm"),
    grid_width = unit(4, "mm")
  )
)
EK1_ht_2

#####################################################################
##EK2_heatmap
##col scale(EK2)
EK2_col_scaled_2 = scale( t(as.matrix(EK2_2nd)),
                          center = TRUE,
                          scale  = TRUE ) 

##col_annotation (EK2)
EK2_col_ant_2 = phy_EK2_2nd %>%
  column_to_rownames(var = "X")

EK2_col_ant_2$Phylum <- factor(EK2_col_ant_2$Phylum, levels = phylum_levels_2nd)

phylum_colors

EK2_col_ha_2 = HeatmapAnnotation(
  df = EK2_col_ant_2,
  col    = list(Phylum = phylum_colors),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  annotation_name_side = "right",
  annotation_legend_param = list(
    Phylum =list(
      title = "Phylum",
      at = phylum_levels_2nd,
      labels = phylum_levels_2nd
    )
  )
)

#plot
EK2_ht_2 = Heatmap(
  EK2_col_scaled_2,
  name            = "EK2",
  column_title = "EK2",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  col = colorRamp2(c(-4, 0, 4), c("#2471A3", "white", "#C0392B")),
  bottom_annotation= EK2_col_ha_2,
  show_row_names  = FALSE,
  show_column_names = FALSE,
  show_heatmap_legend = FALSE,
  cluster_rows    = FALSE,
  cluster_columns = TRUE,
  heatmap_legend_param = list(
    title = NULL,
    legend_direction = "vertical",
    legend_height = unit(6, "cm"),
    grid_width = unit(4, "mm")
  )
)
EK2_ht_2

##################################################################
##EK3_heatmap
##col scale(EK3)
EK3_col_scaled_2 = scale( t(as.matrix(EK3_2nd)),
                          center = TRUE,
                          scale  = TRUE ) 

##col_annotation (EK3)
EK3_col_ant_2 = phy_EK3_2nd %>%
  column_to_rownames(var = "X")

EK3_col_ant_2$Phylum <- factor(EK3_col_ant_2$Phylum, levels = phylum_levels_2nd)

phylum_colors

EK3_col_ha_2 = HeatmapAnnotation(
  df = EK3_col_ant_2,
  col    = list(Phylum = phylum_colors),
  show_annotation_name = FALSE,
  show_legend = TRUE,
  annotation_name_side = "right",
  annotation_legend_param = list(
    Phylum =list(
      title = "Phylum",
      at = phylum_levels_2nd,
      labels = phylum_levels_2nd
    )
  )
)


#plot
EK3_ht_2 = Heatmap(
  EK3_col_scaled_2,
  name            = "EK3",
  column_title = "EK3",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  col = colorRamp2(c(-4, 0, 4), c("#2471A3", "white", "#C0392B")),
  bottom_annotation= EK3_col_ha_2,
  show_row_names  = FALSE,
  show_column_names = FALSE,
  show_heatmap_legend = FALSE,
  cluster_rows    = FALSE,
  cluster_columns = TRUE,
  heatmap_legend_param = list(
    title = NULL,
    legend_direction = "vertical",
    legend_height = unit(6, "cm"),
    grid_width = unit(4, "mm")
  )
)
EK3_ht_2


# Combine 2nd-sampling panels by replicate count (grouped by n):
# n = 10 → AE, CK1, EK4 | n = 8 → EK1, EK2, EK3
ht_2_list_n10 <-AE_ht_2 + CK1_ht_2 + EK4_ht_2
ht_2_list_n8 = EK1_ht_2 + EK2_ht_2 + EK3_ht_2

draw(ht_2_list_n10)
draw(ht_2_list_n8)