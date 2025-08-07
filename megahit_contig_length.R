library(tidyverse)
library(RColorBrewer)

# Load the TSV file
file_path <- "/Users/spencertoth/Desktop/biocode/mega_hit/M997_CORAL1286/quast_combined_2.tsv"
data <- read_tsv(file_path)

# Define contig length categories
contig_metrics <- c(
  "# contigs",
  "# contigs (>= 0 bp)",
  "# contigs (>= 1000 bp)",
  "# contigs (>= 5000 bp)",
  "# contigs (>= 10000 bp)",
  "# contigs (>= 25000 bp)",
  "# contigs (>= 50000 bp)"
)

# Filter relevant rows
data_filtered <- data %>%
  filter(Metric %in% contig_metrics)

# Pivot longer to long format
data_long <- data_filtered %>%
  pivot_longer(
    cols = -Metric,
    names_to = "Sample",
    values_to = "Count"
  ) %>%
  mutate(Count = as.numeric(Count))

# Pivot wider for easier calculation
data_wide <- data_long %>%
  pivot_wider(names_from = Metric, values_from = Count)

# Calculate contig count bins by successive differences
data_bins <- data_wide %>%
  mutate(
    `>= 0 bp` = `# contigs (>= 0 bp)` - `# contigs (>= 1000 bp)`,
    `>= 1000 bp` = `# contigs (>= 1000 bp)` - `# contigs (>= 5000 bp)`,
    `>= 5000 bp` = `# contigs (>= 5000 bp)` - `# contigs (>= 10000 bp)`,
    `>= 10000 bp` = `# contigs (>= 10000 bp)` - `# contigs (>= 25000 bp)`,
    `>= 25000 bp` = `# contigs (>= 25000 bp)` - `# contigs (>= 50000 bp)`,
    `>= 50000 bp` = `# contigs (>= 50000 bp)`
  ) %>%
  select(Sample, `# contigs`, `>= 0 bp`, `>= 1000 bp`, `>= 5000 bp`, `>= 10000 bp`, `>= 25000 bp`, `>= 50000 bp`)

# Reshape to long format for plotting
data_bins_long <- data_bins %>%
  pivot_longer(
    cols = -c(Sample, `# contigs`),
    names_to = "Length_Category",
    values_to = "Count"
  )

# Order Samples by total contigs descending (largest to smallest)
sample_order <- data_bins %>%
  arrange((`# contigs`)) %>%
  pull(Sample)
data_bins_long$Sample <- factor(data_bins_long$Sample, levels = sample_order)

# Set category factor levels ordered shortest to longest contig lengths (bottom to top stacking)
category_order <- c(">= 0 bp", ">= 1000 bp", ">= 5000 bp", ">= 10000 bp", ">= 25000 bp", ">= 50000 bp")
data_bins_long$Length_Category <- factor(data_bins_long$Length_Category, levels = category_order)

# Use a colorblind-friendly palette from RColorBrewer ("Set2")
colors <- brewer.pal(6, "Set2")

# Plot horizontal stacked bar chart
ggplot(data_bins_long, aes(x = Count, y = Sample, fill = Length_Category)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Contig Length Distribution per Sample",
    x = "# contigs (bar length)",
    y = "Sample",
    fill = "Contig Length Category"
  ) +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 6),
    legend.position = "right"
  )

# Categories to classify samples by contig length thresholds
# Define the categories in descending order of length threshold
library(dplyr)

# Define categories ordered longest to shortest (priority order for testing zero)
categories <- c(
  "# contigs (>= 50000 bp)",
  "# contigs (>= 25000 bp)",
  "# contigs (>= 10000 bp)",
  "# contigs (>= 5000 bp)",
  "# contigs (>= 1000 bp)",
  "# contigs (>= 0 bp)"
)

assign_category <- function(row) {
  # Check categories from longest to shortest
  if(!is.na(row[["# contigs (>= 50000 bp)"]]) && row[["# contigs (>= 50000 bp)"]] > 0) {
    return("# contigs (>= 50000 bp)")
  } else if(!is.na(row[["# contigs (>= 50000 bp)"]]) && row[["# contigs (>= 50000 bp)"]] == 0) {
    if(!is.na(row[["# contigs (>= 25000 bp)"]]) && row[["# contigs (>= 25000 bp)"]] > 0) {
      return("# contigs (>= 25000 bp)")
    }
  }
  
  if(!is.na(row[["# contigs (>= 25000 bp)"]]) && row[["# contigs (>= 25000 bp)"]] == 0) {
    if(!is.na(row[["# contigs (>= 10000 bp)"]]) && row[["# contigs (>= 10000 bp)"]] > 0) {
      return("# contigs (>= 10000 bp)")
    }
  }
  
  if(!is.na(row[["# contigs (>= 10000 bp)"]]) && row[["# contigs (>= 10000 bp)"]] == 0) {
    if(!is.na(row[["# contigs (>= 5000 bp)"]]) && row[["# contigs (>= 5000 bp)"]] > 0) {
      return("# contigs (>= 5000 bp)")
    }
  }
  
  if(!is.na(row[["# contigs (>= 5000 bp)"]]) && row[["# contigs (>= 5000 bp)"]] == 0) {
    if(!is.na(row[["# contigs (>= 1000 bp)"]]) && row[["# contigs (>= 1000 bp)"]] > 0) {
      return("# contigs (>= 1000 bp)")
    }
  }
  
  if(!is.na(row[["# contigs (>= 1000 bp)"]]) && row[["# contigs (>= 1000 bp)"]] == 0) {
    if(!is.na(row[["# contigs (>= 0 bp)"]]) && row[["# contigs (>= 0 bp)"]] > 0) {
      return("# contigs (>= 0 bp)")
    }
  }
  
  # Finally, try to assign the longest category with positive count available
  categories_rev <- rev(categories)
  for(cat in categories_rev) {
    if(!is.na(row[[cat]]) && row[[cat]] > 0) {
      return(cat)
    }
  }
  
  # If no positive counts, assign NA
  return(NA_character_)
}


# Apply category assignment row-wise
data_wide <- data_wide %>%
  rowwise() %>%
  mutate(Category = assign_category(cur_data_all())) %>%
  ungroup()

# Check category assignments summary
print(table(data_wide$Category))

library(dplyr)
library(readr)

# Starting from data_wide with samples and a Category column assigned
# e.g. data_wide has columns: Sample, Category

# Read sample names from text file
txt_samples <- readLines("/Users/spencertoth/Desktop/biocode/mega_hit/M997_CORAL1286/M997_CORAL1286_basenames.txt")

# Samples in TSV data
tsv_samples <- data_wide$Sample

# Create dataframe of TSV samples and categories
df_tsv_samples <- data_wide %>%
  select(Sample, Category)

# Find samples in text file but not in TSV
samples_not_in_tsv <- setdiff(txt_samples, tsv_samples)

# Create dataframe for samples not in TSV
df_not_in_tsv <- data.frame(Sample = samples_not_in_tsv, Category = "Not in TSV", stringsAsFactors = FALSE)

# Combine all samples and categories
all_samples <- bind_rows(df_tsv_samples, df_not_in_tsv)

# Sort combined list (optional)
all_samples <- all_samples %>%
  arrange(Category, Sample)

# Write full sample-category mapping to TSV
write_tsv(all_samples, "/Users/spencertoth/Desktop/biocode/mega_hit/M997_CORAL1286/samples_sorted_by_category.tsv")

cat("All samples with categories written to samples_sorted_by_category.tsv\n")
