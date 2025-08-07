library(dplyr)
library(ggplot2)
library(scales)

# File paths
unmapped_fp <- "/Users/spencertoth/Desktop/biocode/barrnap_phyloflash_comparison/Moorea_Corals/CORAL1286_M997_SSU_RRNA - all_barrnap_not_mapped_phyloflash.csv"
mapped_fp <- "/Users/spencertoth/Desktop/biocode/barrnap_phyloflash_comparison/Moorea_Corals/CORAL1286_M997_SSU_RRNA - all_barrnap_phyloflash_mapping.csv"

# Read CSV files
mapped_df <- read.csv(mapped_fp, stringsAsFactors = FALSE)
unmapped_df <- read.csv(unmapped_fp, stringsAsFactors = FALSE)

# Function to count unique qseqid by superkingdom category
count_unique_by_superkingdom <- function(df) {
  df %>%
    group_by(superkingdom) %>%
    summarize(unique_qseqid = n_distinct(qseqid), .groups = "drop")
}

# Summarize counts for mapped ("True 16S") and unmapped ("False 16S")
mapped_counts <- count_unique_by_superkingdom(mapped_df) %>%
  mutate(Category = "YES")

unmapped_counts <- count_unique_by_superkingdom(unmapped_df) %>%
  mutate(Category = "NO")

# Combine results
combined_counts <- bind_rows(mapped_counts, unmapped_counts)

# Calculate percentages for labeling within each Category
combined_counts <- combined_counts %>%
  group_by(Category) %>%
  mutate(
    total = sum(unique_qseqid),
    pct = unique_qseqid / total,
    pct_label = paste0(round(pct * 100, 1), "%")
  ) %>%
  ungroup()

# Create the plot
ggplot(combined_counts, aes(x = Category, y = unique_qseqid, fill = superkingdom)) +
  geom_bar(stat = "identity", color = "black") +       # stacked bars with counts on y-axis
  geom_text(aes(label = pct_label),                     # percentage labels inside segments
            position = position_stack(vjust = 0.5),
            color = "black", size = 4) +
  scale_fill_manual(values = c("Bacteria" = "purple", "Eukaryota" = "green", "N/A" = "gray")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + # little space above bars
  labs(
    x = "MEGAHIT-Barrnap Map to Phyloflash Sample Assembly?",
    y = "Number of MEGAHIT-Barrnap contigs (count)",
    fill = "Superkingdom",
    title = "Moorea Corals (CORAL1286 and M997) 16S SSU rRNA: MEGAHIT-Barrnap Mapped to Phyloflash Contigs"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, face = "bold"))
 