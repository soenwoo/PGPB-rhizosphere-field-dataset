######################################################################################################################
####### Field-scale 16S rRNA dataset on rhizosphere microbial community under Rhodobacter capsulatus treatment ####### 
######################################################################################################################
# Date of creation: 2025-07-13

# Fig_6a beta diversity
# ================================================================================================
# Beta diversity (PCoA) with PERMANOVA and inline p-values
# Workflow
#   1) CSS normalization (metagenomeSeq) → Bray–Curtis distance → PCoA ordination.
#   2) PERMANOVA per sampling round (1st / 2nd) against Sample.site; convert p to star-formatted labels.
#   3) Build PCoA (fill = group, shape = sampling), facet by site_ID, and annotate p-values:
#      - If both rounds exist → print “1st …” and “2nd …”.
#      - If only 1st exists → print the 1st line only.
# ================================================================================================

# ================================================================================================
# Directory Setup for Farm-Specific Analysis
# This script requires farm-specific analysis. 
# For each farm, make sure to set the working directory correctly.
# For example, if analyzing data for "EK1", the working directory should be set as follows:
# setwd("./../Kale/EK1/") 
# Change the directory based on the farm you are analyzing.
# Ensure that the relevant data files are in the corresponding farm directory.
# ================================================================================================


# Load necessary libraries
library(phyloseq)
library(biomformat)
library(tidyverse)
library(vegan) 
library(metagenomeSeq)

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


#### beta diversity ####
# CSS Normalization 
Bacteria_css_norm = phyloseq::phyloseq_to_metagenomeSeq(abundance)
Bacteria_css_otu_norm = metagenomeSeq::MRcounts(
  Bacteria_css_norm, norm = TRUE, log = TRUE)

# Attach normalized OTU table back to Phyloseq object
Bacteria_css = phyloseq(
  sample_data(abundance),
  phy_tree(abundance),
  tax_table(abundance),
  otu_table(Bacteria_css_otu_norm, taxa_are_rows = TRUE))

# Bray-Curtis dissimilarity & PCoA ordination
wunifrac_dist = phyloseq::distance(Bacteria_css, method = "bray")
ordination = ordinate(Bacteria_css, method = "PCoA", distance = wunifrac_dist)


### PERMANOVA (adonis2) per sampling (1st / 2nd) against Group (control vs Treatment)
a = data.frame(sample_data(Bacteria_css))
sampling_levels = sort(unique(sample_data(Bacteria_css)$sampling))
pvalues = list()

# make compact p-value labels with stars
star_for_p = function(pval) {
  if (pval <= 0.001) {
    return(paste0("p<0.001 ***"))
  } else if (pval < 0.01) {
    return(paste0("p=", round(pval, 3), " **"))
  } else if (pval < 0.05) {
    return(paste0("p=", round(pval, 3), "  *"))
  } else {
    return(paste0("p=", round(pval, 3), "   "))
  }
}

# Check & run PERMANOVA for 1st sampling
# If "1" exists in sampling_levels: subset to sampling==1, compute Bray–Curtis distance,
# run PERMANOVA against Sample.site (control vs Treatment), then store a star-formatted p-value.
if (1 %in% sampling_levels) {
  Bacteria_css_sub1 = subset_samples(Bacteria_css, sampling == 1)   # keep only 1st sampling
  dist_sub1 = phyloseq::distance(Bacteria_css_sub1, method = "bray")
  df_sub1 = data.frame(sample_data(Bacteria_css_sub1))              # metadata as data.frame
  
  adonis_sub1 = adonis2(dist_sub1 ~ Sample.site, data = df_sub1, permutations = 999)  # PERMANOVA by group
  pvalues[["1"]] = star_for_p(adonis_sub1$`Pr(>F)`[1])  # convert p-value to compact label with stars
}

# Check & run PERMANOVA for 2nd sampling
# Same logic for sampling==2: distance → PERMANOVA → star label.
if (2 %in% sampling_levels) {
  Bacteria_css_sub2 = subset_samples(Bacteria_css, sampling == 2)   # keep only 2nd sampling
  dist_sub2 = phyloseq::distance(Bacteria_css_sub2, method = "bray")
  df_sub2 = data.frame(sample_data(Bacteria_css_sub2))              # metadata as data.frame
  
  adonis_sub2 = adonis2(dist_sub2 ~ Sample.site, data = df_sub2, permutations = 999)  # PERMANOVA by group
  pvalues[["2"]] = star_for_p(adonis_sub2$`Pr(>F)`[1])  # convert p-value to compact label with stars
}
pvalues

# Prepare data for PCoA scatter; merge ordination axes with metadata
df_meta =  data.frame(sample_data(Bacteria_css)) %>%
  dplyr::rename(Sample_site = Sample.site)
df_ord = data.frame(scores(ordination$vectors)) %>%
  select(1:2)
df_merge = merge(x = df_ord, y = df_meta, by = "row.names", all.x = TRUE)
explained = ordination$values$Relative_eig # change the "axis" title

### Draw PCoA plot
pcoa_plot = ggplot(df_merge, aes(x = Axis.1, y = Axis.2)) +
  geom_point(aes(fill = Sample_site, shape = sampling), 
             col ="gray30", size=4, alpha = 0.8, stroke = 0.3) +
  theme_test(base_size = 8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray10", linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray10", linewidth = 0.4) +
  labs(x = paste("PCoA 1 [", round(explained[1] * 100, 2), "%]", sep = ""),
       y = paste("PCoA 2 [", round(explained[2] * 100, 2), "%]", sep = "")) +
  facet_wrap(.~site_ID, scales = "free_y") +
  theme(strip.text = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 9, face = "bold"),
        axis.title.y = element_text(size = 9, face = "bold")) + 
  scale_fill_manual(name = "Group",
                    values = c("control" = "#A0E515","Treatment" = "#B084FF"),
                    labels = c("control" = "Control","Treatment" = "Treatment")) +
  scale_shape_manual(name = "Sampling",
                     values = c("1" = 21,"2" = 22),
                     labels = c("1" = "1st","2" = "2nd")) +
  guides(fill = guide_legend(override.aes = list(shape = 21, color = "black",
                                                 size = 4, stroke = 0.3)))

pcoa_plot


# Add inline PERMANOVA p-values to the plot
# If both "1" and "2" exist → print two lines (1st & 2nd)
# If only "1" exists → print the 1st line only
if (2 %in% sampling_levels) {
  pcoa_plot_annot = pcoa_plot +
    annotate(
      "text",
      x = -Inf, y = Inf,
      label = paste0("1st sampling: ", pvalues[["1"]], "\n",
                     "2nd sampling: ", pvalues[["2"]]),
      hjust = -0.05, vjust = 1.5, 
      size = 4
    )
} else {
  pcoa_plot_annot = pcoa_plot +
    annotate(
      "text",
      x = -Inf, y = Inf,
      label = paste0("1st sampling: ", pvalues[1]),
      hjust = -0.05, vjust = 1.5,  
      size = 4
    )
  }

pcoa_plot_annot



