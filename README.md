# Rhizosphere microbiome responses to plant growth-promoting bacteria (PGPB) treatment in leafy vegetable farms

## Introduction
This repository contains processed **16S rRNA gene (V3–V4) sequencing data**, **soil physicochemical measurements**, and **IAA production data**, along with **R scripts for downstream analysis**.  
The dataset was generated from field-scale experiments assessing rhizosphere microbial communities under *Rhodobacter capsulatus* treatment.

- **Title**: Rhizosphere microbiome responses to plant growth-promoting bacteria (PGPB) treatment in leafy vegetable farms  
- **Authors**: Seonwoo Choi, Jin Young Kim, Ji-Young Moon, and Junhyun Jeon  

---

## Contents

### Available Datasets

#### **Beta_yield.csv**
- **Description**: Summary table linking β-diversity differences (PERMANOVA p-values) with crop yield changes across farms.  
- **Columns**:  
  - `Site_ID`: Farm identifier  
  - `Plant`: Crop type (mallow, kale, romaine, etc.)  
  - `Beta_p_value`: PERMANOVA p-value (adonis2, Bray–Curtis distances, CSS-/library-normalized counts)  
  - `yield_change`: Relative change in crop yield compared to control  

#### **IAA_production.csv**
- **Description**: Indole-3-acetic acid (IAA) production data measured in isolated *Rhodobacter* strains.  
- **Columns**:  
  - `isolate_ID`: Isolate identifier (e.g., 21PW6, 21PW7, …)  
  - `rep1`, `rep2`, `rep3`: IAA production values (mg/mL) from replicate experiments  

#### **table-w-tax-meta.biom**
- **Description**: `.biom` file generated from the DADA2 pipeline on QIIME2, provided separately for each farm.  
- **Components**:  
  - **OTU table**: Rows = OTUs, Columns = Sample IDs, Values = Abundance counts  
  - **Sample metadata**: `sample-id`, `site_ID`, `week`, `region`, `Sample-site`, `replicate`  
  - **Taxonomy annotation**: Hierarchical classification of OTUs  

#### **rooted_tree.nwk / unrooted_tree.nwk**
- **Description**: Phylogenetic tree files (Newick format), generated using QIIME2 *q2-fragment-insertion*.  

#### **soil_physicochemical**  
- **Files**: `total_soil_1st.csv`, `total_soil_2nd.csv`  
- **Description**: Soil physicochemical properties measured at 1st and 2nd sampling events.  
- **Columns**: Crop type, site ID, treatment (control/treatment), replicate, and soil physicochemical variables (`pH`, `EC`, `OM`, `P<sub>2</sub>O<sub>5</sub>`,`NO<sub>3</sub>-N`,`K`, `Ca`, `Mg`, `Na`)

---

## Analysis Codes
These R scripts were used to generate the figures presented in the paper.  
**Note**: Run `Integrated_preprocessing_rhizosphere_data.R` first before any figure analysis.

### **Preprocessing**
- **Integrated_preprocessing_rhizosphere_data.R**  
  - Must be run before all analyses  
  - Generates outputs: `alpha.csv`, `phylum_RA.csv`, `DESeq2_results.csv`, `beta_diversity_distance.csv`  

### **Figure 1**
- **Fig_1_IAA_production_bar_plot.R**  
  - Visualizes IAA production of *Rhodobacter* isolates  
  - Bar plots with error bars and threshold line  

### **Figure 4**
- **Fig_4a_Relative_abundance_kale.R**  
  - Stacked bar plots of phylum-level relative abundance  
- **Fig_4b_DESeq2_complexheatmap.R**  
  - Heatmaps of differentially abundant genera (DESeq2)  
  - Log-transformed normalized counts, phylum annotation  

### **Figure 5**
- **Fig_5_alpha_diversity.R**  
  - Boxplots of Observed & Shannon diversity  
  - Statistical testing (t-test / Wilcoxon) with significance stars  

### **Figure 6**
- **Fig_6a_beta_diversity.R**  
  - PCoA plots (Bray–Curtis distances)  
  - PERMANOVA p-values annotated  
- **Fig_6b_beta_diversity_yield_association_heatmap.R**  
  - Heatmap linking β-diversity differences (p ≤ 0.05) to yield changes  
  - Counts of farms categorized by Difference/No difference × Increased/Decreased yield  

### **Figure 7**
- **Fig_7_soil_micro_heatmap.R**  
  - Heatmaps integrating soil physicochemical variables with microbial β-diversity changes  
  - Includes significance testing (t-tests, FDR adjustment)  
  - Visualized separately for 1st and 2nd sampling  
