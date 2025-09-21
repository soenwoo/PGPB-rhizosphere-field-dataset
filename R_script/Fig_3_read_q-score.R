library(ShortRead)
library(ggplot2)
library(data.table)

# Set your base directory here
base_dir <- "Q-mean"

# Recursively find all FASTQ or FASTQ.GZ files
fastq_files <- list.files(base_dir, pattern = "\\.fastq(\\.gz)?$", recursive = TRUE, full.names = TRUE)

# Function to compute mean Q-scores
get_mean_qscores <- function(file) {
  fq <- readFastq(file)
  quals <- as(quality(fq), "matrix")
  means <- rowMeans(quals)
  return(data.frame(
    file = file,
    folder = basename(dirname(file)),
    mean_qscore = means
  ))
}

#Process all files and combine

all_data <- rbindlist(lapply(fastq_files, get_mean_qscores))

# Plot using ggplot2
all <- ggplot(all_data, aes(x = mean_qscore, fill = folder)) +
  geom_density(alpha = 0.5) +
  theme_bw(base_size = 16) +
  scale_x_continuous(
    breaks = seq(0,40, by =5),
    limits = c(1,40)
  )+
  labs(
    title = "Distribution of Mean Q-scores (Overall)",
    x = "Mean Q-score of reads",
    y = "Proportion",
    fill = "Crops"
  )
all

ggsave("mean_qscore_all_reads.png", all, width = 28, height = 20.2, units = "cm", dpi = 600)

