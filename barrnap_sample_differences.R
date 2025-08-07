library(ggplot2)
library(dplyr)
library(forcats)

# Path to your input file (update as needed)
file_path <- "/Users/spencertoth/Desktop/biocode/phyloflash/M997_CORAL1286/CORAL1286_M997_SSU_RRNA - 16S_bac_barrnap_assemblies.csv"
#file_path <- "/Users/spencertoth/Desktop/biocode/phyloflash/M997_CORAL1286/phyloflash_M997_CORAL1286_assemblies_revised_06.csv"
data <- read.csv(file_path, stringsAsFactors = FALSE)

# Filter to keep only rows where superkingdom is "Bacteria"
data_bacteria <- data %>%
  filter(superkingdom == "Bacteria")

# Summarize counts of unique query_ids by Sample_ID (using filtered data)
summary_df <- data_bacteria %>%
  group_by(Sample_ID) %>%
  summarise(contig_count = n_distinct(query_id), .groups = "drop")

# Assign bar color based on Sample_ID
summary_df <- summary_df %>%
  mutate(
    bar_color = case_when(
      grepl("_M997$", Sample_ID) ~ "orange",
      grepl("CORAL1286$", Sample_ID) ~ "green",
      TRUE ~ "#CCCCCC" # fallback for all other samples
    )
  )

# Order samples by count
summary_df$Sample_ID <- factor(
  summary_df$Sample_ID,
  levels = summary_df %>% arrange(contig_count) %>% pull(Sample_ID)
)

color_map <- c("orange" = "orange", "green" = "green", "#CCCCCC" = "#CCCCCC")

ggplot(summary_df, aes(x = contig_count, y = Sample_ID, fill = bar_color)) +
  geom_bar(stat = "identity", color = "black", width = 0.8) +
  scale_fill_manual(
    name = "Sample Group",                # Legend title
    values = color_map,                   # Named vector of colors
    breaks = c("orange", "green", "#CCCCCC"),        # Order of legend keys
    labels = c("M997", "CORAL_1286", "Other")  # Friendly legend labels
  ) +
  scale_x_continuous(n.breaks = 20) +  # Increase number of ticks here
  labs(
    title = "MEGA-HIT/Barrnap -- Moorea Coral Samples: Total Unique 16S SSU rRNA Contigs",
    x = "Number of Contigs",
    y = "Moorea Coral Samples"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 5),
    legend.position = "right"
  )

# Summarize total contigs and count of unique Sample_IDs per sample color group
contigs_and_samples_per_group <- summary_df %>%
  group_by(bar_color) %>%
  summarise(
    total_contigs = sum(contig_count),
    unique_samples = n_distinct(Sample_ID),
    .groups = "drop"
  )

# Print the results
cat("Total contigs and unique samples per sample group:\n")
print(contigs_and_samples_per_group)


