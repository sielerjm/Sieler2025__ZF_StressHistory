# plot_significant_taxa.R
# This script creates a bar chart showing the number of significantly abundant taxa
# that are up or down regulated across different experimental effects.
# Input: TSV files containing MaAsLin2 differential abundance results
# Output: Bar chart showing up/down regulated taxa counts

# Load required libraries
library(tidyverse)
library(ggplot2)

# Set seed for reproducibility
set.seed(42)

# Define file paths
file_paths <- c(
  "Code/Analysis/DiffAbund/MaAsLin2__Parasite_effect/significant_results.tsv",
  "Code/Analysis/DiffAbund/MaAsLin2__Antibiotics_parasite/significant_results.tsv",
  "Code/Analysis/DiffAbund/MaAsLin2__Temperature_parasite/significant_results.tsv",
  "Code/Analysis/DiffAbund/MaAsLin2__Combined_parasite/significant_results.tsv"
)

# Define effect names
effect_names <- c("Parasite Effect", "Antibiotics Effect", "Temperature Effect", "Combined Effect")

# Create empty list to store results
results_list <- list()

# Read and process each file
for (i in seq_along(file_paths)) {
  # Read the TSV file
  df <- read_tsv(file_paths[i], show_col_types = FALSE)
  
  # Count up and down regulated taxa
  up_regulated <- sum(df$coef > 0)
  down_regulated <- sum(df$coef < 0)
  
  # Store results
  results_list[[i]] <- data.frame(
    Effect = effect_names[i],
    Regulation = c("Up-regulated", "Down-regulated"),
    Count = c(up_regulated, down_regulated)
  )
}

# Combine all results
results_df <- do.call(rbind, results_list)

# Convert Effect to factor with specified order
results_df$Effect <- factor(results_df$Effect, 
                           levels = c("Parasite Effect", 
                                    "Antibiotics Effect", 
                                    "Temperature Effect", 
                                    "Combined Effect"))

# Create the bar chart
p <- ggplot(results_df, aes(x = Effect, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Up-regulated" = "red", "Down-regulated" = "blue")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  ) +
  labs(
    title = "Number of Significantly Abundant Taxa by Effect",
    x = "Effect",
    y = "Number of Taxa",
    fill = "Regulation"
  )

# Save the plot
ggsave("Code/Analysis/DiffAbund/significant_taxa_barplot.pdf", p, width = 10, height = 6) 