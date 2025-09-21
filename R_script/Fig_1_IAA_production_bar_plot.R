######################################################################################################################
####### Field-scale 16S rRNA dataset on rhizosphere microbial community under Rhodobacter capsulatus treatment ####### 
######################################################################################################################
# Script created by: Seonwoo Choi
# Date of creation: 2025-08-13

# Fig 1: Indole-3-acetic acid (IAA) production by 66 Rhodobacter isolates -------------------------

# Load required libraries
library(ggplot2)
library(dplyr)
library(tibble)

# Set the working directory to the folder where the data file is located
# Make sure the "IAA_production.csv" file is inside the "IAA_production" folder
setwd(".././IAA_production")

#Read CSV, and convert row names to a column
IAA_raw_data = read.csv("IAA_production.csv") %>%
  rownames_to_column(var = "num")

# Calculate the mean and standard error for IAA production (rep1, rep2, rep3)
IAA_raw_data$mean = rowMeans(IAA_raw_data[, 3:5], na.rm = TRUE)
IAA_raw_data$se = apply(IAA_raw_data[, c("rep1", "rep2", "rep3")], 1, function(x) {
  sd(x) / sqrt(length(x))
})

# Sort data by mean values and ensure isolate_ID is ordered correctly for the plot
IAA_raw_data = IAA_raw_data %>%
  arrange(mean) %>%
  mutate(isolate_ID = factor(isolate_ID, levels = isolate_ID))


# Create a bar plot with error bars and a threshold line
ggplot(IAA_raw_data, aes(x = isolate_ID, y = mean)) +
  geom_bar(stat = "identity", fill = "grey40", width = 0.8, alpha = 0.9) + # Bar plot
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.5, size = 0.3, color = "black") + # Error bars
  geom_hline(yintercept = 16.6, color = "red", linetype = "dashed", linewidth = 0.5) +  # Threshold line
  labs(x = "Isolate ID", 
       y = "IAA production (mg/mL)") + # Labels
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    axis.text.y = element_markdown(size = 8),
    axis.title = element_text(size = 9, face = "bold"),
    axis.ticks = element_blank()) + # Customize axis appearance
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     breaks = c(seq(0, 25, 5), 16.6),
                     labels = c("0","5","10","15","20","25",
                                "<span style='color:red;'>16.6</span>")) # y-axis scale


# Save the plot as a PNG file
ggsave("IAA_plot.png", width = 20, height = 8, units = "cm", dpi = 600)

