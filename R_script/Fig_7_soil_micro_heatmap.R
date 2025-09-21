######################################################################################################################
####### Field-scale 16S rRNA dataset on rhizosphere microbial community under Rhodobacter capsulatus treatment ####### 
######################################################################################################################
# Date of creation: 2025-07-13

# ================================================================================================
# IMPORTANT:
# Before running this script, please make sure to execute "Integrated_preprocessing_rhizosphere_data.R".
# That script performs the essential preprocessing steps and generates the input files required here.
# ================================================================================================

# Fig 7 soil physicochemical - rhizosphere heatmaps
# =======================================================================================
# This script creates heatmaps (Fig. 7) to show soil physicochemical 
# responses and microbial community shifts after R. capsulatus treatment.
#
# Workflow:
#   1) Load beta diversity and soil data for each crop  
#   2) Separate data into 1st and 2nd samplings  
#   3) Perform independent t-tests and adjust p-values (FDR)  
#   4) Generate heatmaps with significance stars and annotations (crop, beta diversity)  
#   5) Output final heatmaps for each sampling time
# =======================================================================================

### Load and combine beta diversity results by crop

## Read beta diversity distance tables from each site and add crop label
# Mallow (sites: AB, EC6, EC7, EC8)
AB = read.csv("mallow/AB/output/beta_diversity_distance.csv")
EC6 = read.csv("mallow/EC6/output/beta_diversity_distance.csv")
EC7 = read.csv("mallow/EC7/output/beta_diversity_distance.csv")
EC8 = read.csv("mallow/EC8/output/beta_diversity_distance.csv")
mallow = bind_rows(AB, EC6, EC7, EC8) %>%
  mutate(plant = "mallow")

# Kale (sites: AE, CK1, EK1–EK4)
AE = read.csv("kale/AE/output/beta_diversity_distance.csv")
CK1 = read.csv("kale/CK1/output/beta_diversity_distance.csv")
EK1 = read.csv("kale/EK1/output/beta_diversity_distance.csv")
EK2 = read.csv("kale/EK2/output/beta_diversity_distance.csv")
EK3 = read.csv("kale/EK3/output/beta_diversity_distance.csv")
EK4 = read.csv("kale/EK4/output/beta_diversity_distance.csv")
kale = bind_rows(AE, CK1, EK1, EK2, EK3, EK4) %>%
  mutate(plant = "kale")

# Green lettuce (sites: AC, CL1, EL4, EL6, EL10, EL11)
AC = read.csv("green_lettuce/AC/output/beta_diversity_distance.csv")
CL1 = read.csv("green_lettuce/CL1/output/beta_diversity_distance.csv")
EL4 = read.csv("green_lettuce/EL4/output/beta_diversity_distance.csv")
EL6 = read.csv("green_lettuce/EL6/output/beta_diversity_distance.csv")
EL10 = read.csv("green_lettuce/EL10/output/beta_diversity_distance.csv")
EL11 = read.csv("green_lettuce/EL11/output/beta_diversity_distance.csv")
green = bind_rows(AC, CL1, EL4, EL6, EL10, EL11) %>%
  mutate(plant = "green lettuce")

# Red lettuce (sites: EL1, EL5, EL7, EL8, EL9)
EL1 = read.csv("red_lettuce/EL1/output/beta_diversity_distance.csv")
EL5 = read.csv("red_lettuce/EL5/output/beta_diversity_distance.csv")
EL7 = read.csv("red_lettuce/EL7/output/beta_diversity_distance.csv")
EL8 = read.csv("red_lettuce/EL8/output/beta_diversity_distance.csv")
EL9 = read.csv("red_lettuce/EL9/output/beta_diversity_distance.csv")
red = bind_rows(EL1, EL5, EL7, EL8, EL9) %>%
  mutate(plant = "red lettuce")

# Romaine (sites: AD, EL3)
AD = read.csv("romain/AD/output/beta_diversity_distance.csv")
EL3 = read.csv("romain/EL3/output/beta_diversity_distance.csv")
romain = bind_rows(AD, EL3) %>%
  mutate(plant = "romain")

## Merge all crops into one table
# Final dataset combines all sites from all crops
crop_total = bind_rows(mallow, kale, green, red, romain)

######################################
###          1st sampling          ###
######################################

### beta_disance 1st sampling
beta_s1 = crop_total %>% 
  filter(sampling == 1) %>% 
  select(crop = plant,
         site_ID,
         diff = distance,
         p_value = p.value) %>% 
  mutate(variable = "Beta")

beta_s2 = crop_total %>% 
  filter(sampling == 2) %>% 
  select(crop = plant,
         site_ID,
         diff = distance,
         p_value = p.value) %>% 
  mutate(variable = "Beta")

# soil 1st sampling
soil_s1 <- read.csv("soil_physicochemical/total_soil_1st_v3.csv", stringsAsFactors = FALSE)


vars <- c("pH","EC","OM","P2O5","NO3.N","K","Ca","Mg","Na")
options(scipen = 999)

### Run independent t-tests (control vs treatment) for each variable and site_ID (farm)
soil_s1_ttest <- map_dfr(vars, function(v) {
  soil_s1 %>%
    split(.$site_ID) %>%
    map_dfr(function(x) {
      ctrl <- x[[v]][x$Sample.site == "control"]   # values for control group
      trt  <- x[[v]][x$Sample.site == "treatment"] # values for treatment group
      tt   <- t.test(trt, ctrl, paired = FALSE)    # independent two-sample t-test
      tibble(
        crop           = unique(x$Crop),
        site_ID        = unique(x$site_ID),
        variable       = v,
        mean_control   = mean(ctrl, na.rm = TRUE),
        mean_treatment = mean(trt,  na.rm = TRUE),
        diff           = (mean(trt, na.rm = TRUE) - mean(ctrl, na.rm = TRUE)) /
                          mean(ctrl, na.rm = TRUE) * 100,   # relative difference (%)
        p.value        = tt$p.value
      )
    })
})


# Adjust p-values for multiple testing (Benjamini–Hochberg FDR) 
# and assign significance stars
soil_s1_tadj = soil_s1_ttest %>% 
  mutate(
    p_adjust = p.adjust(p.value, method = "BH"),  # adjusted p-value
    star = cut(p_adjust,
               breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),  # significance thresholds
               labels = c("***", "**", "*", "")    # significance stars
               )
    )

# Create a final summary table (only key variables)
soil_s1_tsum = soil_s1_tadj %>% 
  select(crop, site_ID, variable, diff,
         p_adjust, star) %>%
  remove_rownames() 

### Heatmap data (sampling 1)

# Create a heatmap-ready matrix from soil physicochemical data
hmap_s1w = soil_s1_tsum %>%  
  select(site_ID, variable, diff) %>%     # select only site_ID, variable, and diff
  pivot_wider(
    names_from = variable, 
    values_from = diff) %>%  # reshape data: long → wide (variables as columns)
  arrange(site_ID) %>%  # order rows by site_ID
  column_to_rownames(var = "site_ID") # set site_ID as row names 

# Column-wise z-score normalization of soil variables
hmap_s1_z = scale(hmap_s1w)

### Right annotation bar for the heatmap
# Use beta diversity results (relative difference, %) from sampling 1
bar_s1_df <- beta_s1 %>%
  select(site_ID, diff)

# Convert the two-column data frame (site_ID → diff) into a named vector
bar_s1_vec = deframe(bar_s1_df)
# Reorder vector values to match the rownames (site_ID) of the heatmap matrix
# This ensures correct alignment between the annotation bar and the heatmap rows
bar_s1_vec = bar_s1_vec[rownames(hmap_s1_z)]

## Right-side annotation bar (beta diversity)
a1_r <- rowAnnotation(
  "Community\ndistance" = anno_barplot(
    bar_s1_vec,                             # named vector (site_ID → relative difference)
    gp = gpar(fill = "#7FBF7B", col = NA), # bar color and no border
    border = TRUE,                         # draw border around bars
    axis = TRUE,                           # display axis
    ylim = c(0, max(bar_s1_vec) * 1.2),       # expand y-axis range for better spacing
    width = unit(3.5, "cm"),               # set annotation bar width
    axis_param = list(                     # customize axis appearance
      gp = gpar(
        fontsize = 10,                     # axis label font size
        lwd = 2)                            # axis line width
      )
    )
  )

### Left annotation strip for the heatmap
# Use crop categories (sampling 1) to show site-specific crop information

# Create a mapping between site_ID and crop (one crop per site)
crop_s1_vec = soil_s1_tsum %>% 
  distinct(site_ID, crop) %>% 
  deframe()


# Reorder to match the heatmap row order (important for alignment)
crop_s1_vec = crop_s1_vec[rownames(hmap_s1_z)]

# Build color palette for crop categories
lvls_s1 <- unique(crop_s1_vec)   
cols_s1 <- as.character(paletteer::paletteer_d("RColorBrewer::Set2", length(lvls_s1)))
pal_crop_s1 <- setNames(cols_s1, lvls_s1)

# Define left annotation bar (categorical: crop)
a1_l <- rowAnnotation(
  crop = crop_s1_vec,               # categorical annotation
  col  = list(crop = pal_crop_s1),   # assign colors to crop categories
  simple_anno_size = unit(3, "mm"),     # thickness
  annotation_legend_param = list(crop = list(title = "Crop"))  # legend settings
)

###Heatmap with Statistical Significance Annotations 
# Create a star matrix (site_ID × variable)  
# Each cell stores the significance level symbol (*, **, ***)
soil_s1_star <- soil_s1_tsum %>%
  distinct(site_ID, variable, star) %>%                          
  pivot_wider(names_from = variable, values_from = star) %>%
  column_to_rownames(var = "site_ID") %>% 
  as.matrix()


# Reorder matrix rows and columns to match the heatmap (hmap_scale)
# This ensures stars are correctly overlaid in the heatmap cells
soil_s1_star = soil_s1_star[rownames(hmap_s1_z), colnames(hmap_s1_z), drop = FALSE]

### plot soil heatmap
# Helper function: convert variable names into plotmath expressions 
# (e.g., "P2O5" → "P[2]*O[5]", "NO3.N" → "NO[3]-N")
pretty_var <- function(x) {
  map <- c(
    "pH"    = "pH",
    "EC"    = "EC",
    "OM"    = "OM",
    "P2O5"  = "P[2]*O[5]",
    "NO3.N" = "NO[3]-N",
    "K"     = "K",
    "Ca"    = "Ca",
    "Mg"    = "Mg",
    "Na"    = "Na"
  )
  as.expression(parse(text = map[x]))
}

col_s1_labs <- pretty_var(colnames(hmap_s1_z))


## Draw the heatmap
soil_s1_ht <- Heatmap(
  hmap_s1_z,
  name = "Relative\nchange (%)",  # legend title
  column_labels = col_s1_labs,    # nicer chemical labels
  right_annotation = a1_r,        # beta diversity bar annotation
  left_annotation = a1_l,         # crop category annotation
  show_row_names  = T,
  cluster_rows = TRUE,    
  cluster_columns = TRUE, 
  # Overlay significance stars in each heatmap cell
  cell_fun = function(j, i, x, y, w, h, fill) {
    lab <- soil_s1_star[i, j]
    if (!is.na(lab) && nzchar(lab)) {     # "", NA는 출력 안 함
      grid.text(lab, x = x, y = y,
                gp = gpar(fontsize = 12, fontface = "bold"))
    }
  }
)

soil_s1_ht

### Add significance stars to the beta diversity distance bar plot
# based on PERMANOVA p-values from sampling 1

# Assign significance levels (*, **, ***) from PERMANOVA p-values
bar_s1_star = beta_s1 %>% 
  mutate(star = cut(p_value,
                    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                    labels = c("***", "**", "*", "")
  ))

# Create a named vector: site_ID -> star symbol
sig_s1_vec = bar_s1_star %>% 
  select(site_ID, star) %>% 
  deframe()

# Reorder to match heatmap row order (important for alignment)
sig_s1_vec = sig_s1_vec[rownames(hmap_s1_z)]

# Draw the heatmap to fix the actual row order (after clustering)
dr_s1_ht <- draw(soil_s1_ht)
ord_s1 <- row_order(dr_s1_ht) # extract row order from heatmap object

# Get star labels only for significant rows
labs_s1 <- as.character(sig_s1_vec[ord_s1])


# Add stars next to bar plot (right annotation)
decorate_annotation("Community\ndistance", slice = 1, {
  grid.text(
    label = rev(labs_s1),
    x = unit(rev(bar_s1_vec[ord_s1]), "native") + unit(4, "mm"),
    y = unit(seq_along(ord_s1), "native"),
    gp = gpar(fontsize = 12, fontface = "bold")
  )
})

######################################
###          2nd sampling          ###
######################################
# The following code is identical in structure to the 1st sampling.
# keep the same workflow and logic.

# soil 2nd sampling
soil_s2 <- read.csv("soil_physicochemical/total_soil_2nd_v3.csv", stringsAsFactors = FALSE)


vars <- c("pH","EC","OM","P2O5","NO3.N","K","Ca","Mg","Na")
options(scipen = 999)

### Run independent t-tests (control vs treatment) for each variable and site_ID (farm)
soil_s2_ttest <- map_dfr(vars, function(v) {
  soil_s2 %>%
    split(.$site_ID) %>%
    map_dfr(function(x) {
      ctrl <- x[[v]][x$Sample.site == "control"]   # values for control group
      trt  <- x[[v]][x$Sample.site == "treatment"] # values for treatment group
      tt   <- t.test(trt, ctrl, paired = FALSE)    # independent two-sample t-test
      tibble(
        crop           = unique(x$Crop),
        site_ID        = unique(x$site_ID),
        variable       = v,
        mean_control   = mean(ctrl, na.rm = TRUE),
        mean_treatment = mean(trt,  na.rm = TRUE),
        diff           = (mean(trt, na.rm = TRUE) - mean(ctrl, na.rm = TRUE)) /
          mean(ctrl, na.rm = TRUE) * 100,   # relative difference (%)
        p.value        = tt$p.value
      )
    })
})


# Adjust p-values for multiple testing (Benjamini–Hochberg FDR) 
# and assign significance stars
soil_s2_tadj = soil_s2_ttest %>% 
  mutate(
    p_adjust = p.adjust(p.value, method = "BH"),  # adjusted p-value
    star = cut(p_adjust,
               breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),  # significance thresholds
               labels = c("***", "**", "*", "")    # significance stars
    )
  )

# Create a final summary table (only key variables)
soil_s2_tsum = soil_s2_tadj %>% 
  select(crop, site_ID, variable, diff,
         p_adjust, star) %>%
  remove_rownames() 

### Heatmap data (sampling 2)

# Create a heatmap-ready matrix from soil physicochemical data
hmap_s2w = soil_s2_tsum %>%  
  select(site_ID, variable, diff) %>%     # select only site_ID, variable, and diff
  pivot_wider(
    names_from = variable, 
    values_from = diff) %>%  # reshape data: long → wide (variables as columns)
  arrange(site_ID) %>%  # order rows by site_ID
  column_to_rownames(var = "site_ID") # set site_ID as row names 

# Column-wise z-score normalization of soil variables
hmap_s2_z = scale(hmap_s2w)

### Right annotation bar for the heatmap
# Use beta diversity results (relative difference, %) from sampling 2
bar_s2_df <- beta_s2 %>%
  select(site_ID, diff)

# Convert the two-column data frame (site_ID → diff) into a named vector
bar_s2_vec = deframe(bar_s2_df)
# Reorder vector values to match the rownames (site_ID) of the heatmap matrix
# This ensures correct alignment between the annotation bar and the heatmap rows
bar_s2_vec = bar_s2_vec[rownames(hmap_s2_z)]

## Right-side annotation bar (beta diversity)
a2_r <- rowAnnotation(
  "Community\ndistance" = anno_barplot(
    bar_s2_vec,                             # named vector (site_ID → relative difference)
    gp = gpar(fill = "#7FBF7B", col = NA), # bar color and no border
    border = TRUE,                         # draw border around bars
    axis = TRUE,                           # display axis
    ylim = c(0, max(bar_s2_vec) * 1.2),       # expand y-axis range for better spacing
    width = unit(3.5, "cm"),               # set annotation bar width
    axis_param = list(                     # customize axis appearance
      gp = gpar(
        fontsize = 10,                     # axis label font size
        lwd = 2)                            # axis line width
    )
  )
)

### Left annotation strip for the heatmap
# Use crop categories (sampling 2) to show site-specific crop information

# Create a mapping between site_ID and crop (one crop per site)
crop_s2_vec = soil_s2_tsum %>% 
  distinct(site_ID, crop) %>% 
  deframe()


# Reorder to match the heatmap row order (important for alignment)
crop_s2_vec = crop_s2_vec[rownames(hmap_s2_z)]

# Build color palette for crop categories
lvls_s2 <- unique(crop_s2_vec)   
cols_s2 <- as.character(paletteer::paletteer_d("RColorBrewer::Set2", length(lvls_s2)))
pal_crop_s2 <- setNames(cols_s2, lvls_s2)

# Define left annotation bar (categorical: crop)
a2_l <- rowAnnotation(
  crop = crop_s2_vec,               # categorical annotation
  col  = list(crop = pal_crop_s2),   # assign colors to crop categories
  simple_anno_size = unit(3, "mm"),     # thickness
  annotation_legend_param = list(crop = list(title = "Crop"))  # legend settings
)

###Heatmap with Statistical Significance Annotations 
# Create a star matrix (site_ID × variable)  
# Each cell stores the significance level symbol (*, **, ***)
soil_s2_star <- soil_s2_tsum %>%
  distinct(site_ID, variable, star) %>%                          
  pivot_wider(names_from = variable, values_from = star) %>%
  column_to_rownames(var = "site_ID") %>% 
  as.matrix()


# Reorder matrix rows and columns to match the heatmap (hmap_scale)
# This ensures stars are correctly overlaid in the heatmap cells
soil_s2_star = soil_s2_star[rownames(hmap_s2_z), colnames(hmap_s2_z), drop = FALSE]

### plot soil heatmap
# Helper function: convert variable names into plotmath expressions 
# (e.g., "P2O5" → "P[2]*O[5]", "NO3.N" → "NO[3]-N")
pretty_var <- function(x) {
  map <- c(
    "pH"    = "pH",
    "EC"    = "EC",
    "OM"    = "OM",
    "P2O5"  = "P[2]*O[5]",
    "NO3.N" = "NO[3]-N",
    "K"     = "K",
    "Ca"    = "Ca",
    "Mg"    = "Mg",
    "Na"    = "Na"
  )
  as.expression(parse(text = map[x]))
}

col_s2_labs <- pretty_var(colnames(hmap_s2_z))


## Draw the heatmap
soil_s2_ht <- Heatmap(
  hmap_s2_z,
  name = "Relative\nchange (%)",  # legend title
  column_labels = col_s2_labs,    # nicer chemical labels
  right_annotation = a2_r,        # beta diversity bar annotation
  left_annotation = a2_l,         # crop category annotation
  show_row_names  = T,
  cluster_rows = TRUE,    
  cluster_columns = TRUE, 
  # Overlay significance stars in each heatmap cell
  cell_fun = function(j, i, x, y, w, h, fill) {
    lab <- soil_s2_star[i, j]
    if (!is.na(lab) && nzchar(lab)) {     # "", NA는 출력 안 함
      grid.text(lab, x = x, y = y,
                gp = gpar(fontsize = 12, fontface = "bold"))
    }
  }
)

soil_s2_ht

### Add significance stars to the beta diversity distance bar plot
# based on PERMANOVA p-values from sampling 2

# Assign significance levels (*, **, ***) from PERMANOVA p-values
bar_s2_star = beta_s2 %>% 
  mutate(star = cut(p_value,
                    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                    labels = c("***", "**", "*", "")
  ))

# Create a named vector: site_ID -> star symbol
sig_s2_vec = bar_s2_star %>% 
  select(site_ID, star) %>% 
  deframe()

# Reorder to match heatmap row order (important for alignment)
sig_s2_vec = sig_s2_vec[rownames(hmap_s2_z)]

# Draw the heatmap to fix the actual row order (after clustering)
dr_s2_ht <- draw(soil_s2_ht)
ord_s2 <- row_order(dr_s2_ht) # extract row order from heatmap object

# Get star labels only for significant rows
labs_s2 <- as.character(sig_s2_vec[ord_s2])


# Add stars next to bar plot (right annotation)
decorate_annotation("Community\ndistance", slice = 1, {
  grid.text(
    label = rev(labs_s2),
    x = unit(rev(bar_s2_vec[ord_s2]), "native") + unit(4, "mm"),
    y = unit(seq_along(ord_s2), "native"),
    gp = gpar(fontsize = 12, fontface = "bold")
  )
})
