######################################################################################################################
####### Field-scale 16S rRNA dataset on rhizosphere microbial community under Rhodobacter capsulatus treatment ####### 
######################################################################################################################
# Date of creation: 2025-08-13

# =========================================================
# Integrated Preprocessing Code
# This code should be run before any analyses such as fig4, fig5, and fig7.
# It preprocesses and filters the data to prepare it for further analysis.
# Please run this code first before starting any specific analysis.
# =========================================================

# =========================================================
# Directory Setup for Farm-Specific Analysis
# This script requires farm-specific analysis. 
# For each farm, make sure to set the working directory correctly.
# For example, if analyzing data for "EK1", the working directory should be set as follows:
# setwd("./../Kale/EK1/") 
# Change the directory based on the farm you are analyzing.
# Ensure that the relevant data files are in the corresponding farm directory.
# =========================================================

# Example for setting working directory for a specific farm:
# setwd("./../Kale/EK1/")  # Set the working directory to the EK1 folder
# setwd("./../green_lettuce/EL4/")  # Set the working directory to the EL4 folder

# Load necessary libraries
library(phyloseq)
library(biomformat)
library(tidyverse)
library(tibble)
library(vegan)

## Creates directories to save data if they don't exist
if (any(list.files() == "output")) {
  print("The 'output' directory already exists.")
} else {
  print("Creating 'output' directory.")
  dir.create("output")
}

### READ THE REQUIRED OTU, TAXONOMY, META DATA FILES 
biom_file = import_biom(biomformat::read_biom(biom_file = "table-w-tax-meta.biom"))
biom_file

tree_file = read_tree("rooted_tree.nwk")
ps_data = merge_phyloseq(biom_file, tree_file)

sample_data(ps_data)

### Change the Taxonomy Ranks
colnames(tax_table(ps_data)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

### Remove the patterns of taxonomy classification names [k__/p__/c__/o__/f__/g__/s__]
tax_table(ps_data)[, colnames(tax_table(ps_data))] =
  gsub(tax_table(ps_data)[, colnames(tax_table(ps_data))],
       pattern = "[A-Z]_[0-7]__", replacement = "")

### Data preprocessing steps
## If you wish to apply any filters at any taxonomy level
## then please uncomment respective lines below and run accordingly
tax_filt = ps_data %>%
  subset_taxa(Kingdom != "Archaea") %>%         # Remove Archaea
  subset_taxa(Kingdom != "Unassigned") %>%
  subset_taxa(Phylum != "NA") %>%               # Remove NA (Not assigned)
  subset_taxa(Order != "Chloroplast") %>%       # Remove Chloroplast reads
  subset_taxa(Family != "Mitochondria")         # Remove Mitochondria reads

tax_filt

## Normalize to Even depth
median = median(sample_sums(tax_filt))
standf = function(x, t = median) round(t * (x / sum(x)))
norm = transform_sample_counts(tax_filt, standf)

## Remove taxa not seen more than 3 times in at least 20% of the samples.
## This protects against an OTU with small mean & trivially large C.V.
abundance = filter_taxa(tax_filt, function(x) sum(x > 3) > (0.2*length(x)), TRUE)

# Add sampling information
abundance_meta = sample_data(abundance)
abundance_meta$sampling= ifelse(abundance_meta$week == "4", "1", "2")
abundance_meta
sample_data(abundance) = abundance_meta


# =========================================================
# Fig 4A Preprocessing
# =========================================================
# Preprocessing for Fig 4A: Relative Abundance at Phylum Level
# This code processes the relative abundance (RA) data at the Phylum level.
# It normalizes the abundance, groups the data at the Phylum level, 
# and filters out taxa with RA below 0.01.
# The processed data is saved as "phylum_RA.csv" for further use in Fig 4A.

#### Relative Abundance
phylum_RA = transform_sample_counts(abundance, function(x) x/sum(x)) %>%
  tax_glom(taxrank="Phylum") %>%    # Phylum Level
  psmelt() %>%
  filter(Abundance > 0.01)   # Filter RA above 0.01

write.csv(phylum_RA, file="output/phylum_RA.csv")

# =========================================================
# Fig 4B Preprocessing
# =========================================================
# Preprocessing for Fig 4B: Analysis for 1st Sampling and 2nd Sampling
# This code performs preprocessing for Fig 4B, where we separately analyze 
# data from the first and second sampling events.
# The first part handles the preprocessing for the 1st sampling, 
# including normalization, filtering, DESeq2 analysis, and significant taxa selection.
# The processed data is saved for further use in the analysis.

# Load necessary libraries
library("DESeq2") 

# Perform initial preprocessing for the first sampling
abundance

# Group data by Genus level
tax_level = tax_glom(abundance, "Genus")

######################################
###          1st sampling          ###
######################################

# Filter the data to include only the 1st sampling group
tax_1st = subset_samples(tax_level, sampling == "1")

sample_data(tax_1st)

# Prepare data for DESeq2 analysis
ps_ds_1st  = phyloseq_to_deseq2(tax_1st, ~ Sample.site )

# Estimate size factors and perform DESeq2 analysis using Wald test
dds_1st = estimateSizeFactors(ps_ds_1st , type="poscounts")
dds_1st = DESeq(dds_1st, test="Wald", fitType="parametric")

# Extract results
res_1st = results(dds_1st, cooksCutoff=T)

# Map OTU names to genus names in the taxonomy table
res_df_1st = as.data.frame(res_1st) %>% 
  rownames_to_column(var = "OTU")

# Retrieve the genus, phylum, and family information for OTUs
tax_name_1st = tax_table(tax_level) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "OTU") %>% 
  filter(OTU %in% res_df_1st$OTU) %>% 
  select(OTU, Genus, Phylum, Family)

# Merge DESeq2 results with taxonomy information
merge_df_1st = merge(res_df_1st, tax_name_1st, by = "OTU")


# Filter out "uncultured" and ambiguous taxa
merge_fl_1st = merge_df_1st[!grepl("uncultured", merge_df_1st$Genus) &
                          !grepl("Ambiguous_taxa", merge_df_1st$Genus) &
                          !grepl("metagenome", merge_df_1st$Genus), ]

# Select significant taxa based on p-value and log2FoldChange thresholds
alpha = 0.05
log2fc_threshold = 1
merge_fl_1st$log2FoldChange = as.numeric(merge_fl_1st$log2FoldChange)
merge_fl_1st$padj = as.numeric(merge_fl_1st$padj)

# Filter for significant taxa
t_sigtab_1st = merge_fl_1st[which(merge_fl_1st$padj < alpha & 
                                    abs(merge_fl_1st$log2FoldChange) > log2fc_threshold), ]

# Extract normalized counts for significant taxa
norm_1st = as.data.frame(counts(dds_1st, normalized = T)) %>% 
  rownames_to_column(var = "OTU")

# Merge normalized counts with significant taxa
sig_1st = t_sigtab_1st %>%
  left_join(norm_1st, by = c("OTU"))

# Convert to matrix for heatmap generation
hmap_data_1st = sig_1st %>% 
  column_to_rownames(var = "Genus") %>% 
  select(-c(Phylum, Family, OTU,
            baseMean, log2FoldChange, lfcSE, stat, pvalue, padj))

# Convert to matrix format for heatmap plotting
mat_hmap_1st = as.matrix(hmap_data_1st)

# Log-transform for visualization
log_hmap_1st = log(mat_hmap_1st + 1)
rownames(log_hmap_1st)

#### Phylum, Family annotation transformation
tax_annotaion_1st = sig_1st %>% 
  column_to_rownames(var = "Genus") %>% 
  select(Phylum)


# Save the results for the first sampling
write.csv(tax_annotaion_1st, file = "output/1st_phylum_annotation_DESeq2.csv", row.names = TRUE,
          quote = FALSE)
write.csv(log_hmap_1st, file = "output/1st_sampling_Genus_DESeq2_results.csv", row.names = TRUE,
          quote = FALSE)

######################################
###          2nd sampling          ###
######################################
# Run this section only if "2" exists in the sampling metadata
if ("2" %in% unique(sample_data(abundance)$sampling)) {
  # Filter the data to include only the 1st sampling group
  tax_2nd = subset_samples(tax_level, sampling == "2")
  
  sample_data(tax_2nd)
  
  # Prepare data for DESeq2 analysis
  ps_ds_2nd  = phyloseq_to_deseq2(tax_2nd, ~ Sample.site )
  
  # Estimate size factors and perform DESeq2 analysis using Wald test
  dds_2nd = estimateSizeFactors(ps_ds_2nd , type="poscounts")
  dds_2nd = DESeq(dds_2nd, test="Wald", fitType="parametric")
  
  # Extract results
  res_2nd = results(dds_2nd, cooksCutoff=T)
  
  # Map OTU names to genus names in the taxonomy table
  res_df_2nd = as.data.frame(res_2nd) %>% 
    rownames_to_column(var = "OTU")
  
  # Retrieve the genus, phylum, and family information for OTUs
  tax_name_2nd = tax_table(tax_level) %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "OTU") %>% 
    filter(OTU %in% res_df_2nd$OTU) %>% 
    select(OTU, Genus, Phylum, Family)
  
  # Merge DESeq2 results with taxonomy information
  merge_df_2nd = merge(res_df_2nd, tax_name_2nd, by = "OTU")
  
  
  # Filter out "uncultured" and ambiguous taxa
  merge_fl_2nd = merge_df_2nd[!grepl("uncultured", merge_df_2nd$Genus) &
                                !grepl("Ambiguous_taxa", merge_df_2nd$Genus) &
                                !grepl("metagenome", merge_df_2nd$Genus), ]
  
  # Select significant taxa based on p-value and log2FoldChange thresholds
  alpha = 0.05
  log2fc_threshold = 1
  merge_fl_2nd$log2FoldChange = as.numeric(merge_fl_2nd$log2FoldChange)
  merge_fl_2nd$padj = as.numeric(merge_fl_2nd$padj)
  
  # Filter for significant taxa
  t_sigtab_2nd = merge_fl_2nd[which(merge_fl_2nd$padj < alpha & 
                                      abs(merge_fl_2nd$log2FoldChange) > log2fc_threshold), ]
  
  # Extract normalized counts for significant taxa
  norm_2nd = as.data.frame(counts(dds_2nd, normalized = T)) %>% 
    rownames_to_column(var = "OTU")
  
  # Merge normalized counts with significant taxa
  sig_2nd = t_sigtab_2nd %>%
    left_join(norm_2nd, by = c("OTU"))
  
  # Convert to matrix for heatmap generation
  hmap_data_2nd = sig_2nd %>% 
    column_to_rownames(var = "Genus") %>% 
    select(-c(Phylum, Family, OTU,
              baseMean, log2FoldChange, lfcSE, stat, pvalue, padj))
  
  # Convert to matrix format for heatmap plotting
  mat_hmap_2nd = as.matrix(hmap_data_2nd)
  
  # Log-transform for visualization
  log_hmap_2nd = log(mat_hmap_2nd + 1)
  rownames(log_hmap_2nd)
  
  #### Phylum, Family annotation transformation
  tax_annotaion_2nd = sig_2nd %>% 
    column_to_rownames(var = "Genus") %>% 
    select(Phylum)
  
  # Save the results for the first sampling
  write.csv(tax_annotaion_2nd, file = "output/2nd_phylum_annotation_DESeq2.csv", row.names = TRUE,
            quote = FALSE)
  write.csv(log_hmap_2nd, file = "output/2nd_sampling_Genus_DESeq2_results.csv", row.names = TRUE,
            quote = FALSE)
} else {
  message("No 2nd sampling data detected. Skipping analysis.")
}

# =========================================================
# Fig 5 Preprocessing, Analysis
# =========================================================
# Fig 5 Analysis: Alpha Diversity (Richness Estimation)
# The results are then saved as a CSV file for further use.

### Rarefaction
# Rarefaction is applied to even out the depth of sequencing across samples.
# This process will randomly subsample the OTU tables to a common depth.
# `rarefy_even_depth` ensures all samples have the same sequencing depth (90% of the smallest sample's reads).
# A random seed is set to ensure reproducibility.
rare.alpha = rarefy_even_depth(abundance, rngseed=1, sample.size=0.9*min(sample_sums(abundance)), replace=F)
sample_sums(rare.alpha)

### Richness Estimation
alpha = estimate_richness(
  rare.alpha,
  split=TRUE)  # TRUE -> separate richness estimates for each sample


head(alpha)

meta_data = data.frame(sample_data(rare.alpha)) %>% 
  rownames_to_column(var = "sample_ID") %>% 
  select(sample_ID,site_ID, Sample.site, replicate, sampling)

df_alpha = alpha %>% 
  rownames_to_column(var = "sample_ID") 

merge_alpha = merge(meta_data, df_alpha, by = "sample_ID")

write.csv(merge_alpha, file="output/alpha.csv")

# =========================================================
# Fig 7 Preprocessing, Analysis
# =========================================================
# Fig 7 Analysis : Beta Diversity Analysis(Distance Matrix & PCoA)
# This script performs beta diversity analysis including:
#   (1) CSS normalization of OTU counts
#   (2) Re-construction of the phyloseq object
#   (3) Bray-Curtis distance calculation and PCoA ordination
#   (4) Extraction of PC1/PC2 coordinates and merging with metadata
#   (5) Calculation of ΔPCoA distances between Control and Treatment
#   (6) PERMANOVA testing for each sampling event
#   (7) Export of final results (ΔPCoA distances + p-values)

### Distance Matrix & PCoA
# CSS Normalization 
Bacteria_css_norm = phyloseq::phyloseq_to_metagenomeSeq(abundance)
Bacteria_css_otu_norm = metagenomeSeq::MRcounts(
  Bacteria_css_norm, norm = TRUE, log = TRUE
)

# Attach normalized OTU table back to Phyloseq object
Bacteria_css = phyloseq(
  sample_data(abundance),
  phy_tree(abundance),
  tax_table(abundance),
  otu_table(Bacteria_css_otu_norm, taxa_are_rows = TRUE)
)

# Bray-Curtis dissimilarity & PCoA ordination
wunifrac_dist = phyloseq::distance(Bacteria_css, method = "bray")
ordination = ordinate(Bacteria_css, method = "PCoA", distance = wunifrac_dist)

### PCoA Distance Calculation
# Extract PC1 and PC2 scores
pc1_2_scores = as.data.frame(ordination$vectors[, 1:2]) %>%
  tibble::rownames_to_column(var = "rownames")

# Merge sample metadata with PC1, PC2 scores
df_samples = as(sample_data(abundance), "data.frame") 

df_merge = df_samples %>%
  tibble::rownames_to_column(var = "rownames") %>%
  left_join(pc1_2_scores, by = "rownames") %>%
  dplyr::select(site_ID, sampling, Sample.site, replicate,
                pc1 = Axis.1, pc2 = Axis.2)

# Average PC1 and PC2 values within groups
df_avg = df_merge %>%
  group_by(site_ID, sampling, Sample.site) %>%
  summarise(avg_pc1 = mean(pc1),
            avg_pc2 = mean(pc2), .groups = "drop")

# Calculate PCoA distance between Control and Treatment
df_avg_dist = df_avg %>%
  pivot_wider(
    id_cols     = c(site_ID, sampling),
    names_from  = Sample.site,
    values_from = c(avg_pc1, avg_pc2),
    names_sep   = "_"
  ) %>%
  mutate(distance = sqrt(
    (avg_pc1_Treatment - avg_pc1_control)^2 +
      (avg_pc2_Treatment - avg_pc2_control)^2
  ))

### PERMANOVA Testing
# Run PERMANOVA for each sampling level
meta_df = data.frame(sample_data(Bacteria_css))
sampling_levels = sort(unique(sample_data(Bacteria_css)$sampling))
pvalues = list()

if (1 %in% sampling_levels) {
  Bacteria_css_sub1 = subset_samples(Bacteria_css, sampling == 1)
  dist_sub1 = phyloseq::distance(Bacteria_css_sub1, method = "bray", weighted = TRUE)
  df_sub1 = data.frame(sample_data(Bacteria_css_sub1))
  
  adonis_sub1 = adonis2(dist_sub1 ~ Sample.site, data = df_sub1, permutations = 999)
  pvalues[["1"]] = adonis_sub1$`Pr(>F)`[1]
}

if (2 %in% sampling_levels) {
  Bacteria_css_sub2 = subset_samples(Bacteria_css, sampling == 2)
  dist_sub2 = phyloseq::distance(Bacteria_css_sub2, method = "bray", weighted = TRUE)
  df_sub2 = data.frame(sample_data(Bacteria_css_sub2))
  
  adonis_sub2 = adonis2(dist_sub2 ~ Sample.site, data = df_sub2, permutations = 999)
  pvalues[["2"]] = adonis_sub2$`Pr(>F)`[1]
}

# Combine p-values with ΔPCoA distances
df_pvalue = data.frame(
  site_ID  = unique(meta_df$site_ID),
  sampling = as.numeric(names(pvalues)),
  p.value  = unlist(pvalues)
)

df_merg = merge(df_avg_dist, df_pvalue, by = c("site_ID", "sampling"))

### Save results
df_merg
write.csv(df_merg, file = "output/beta_diversity_distance.csv", row.names = FALSE)

