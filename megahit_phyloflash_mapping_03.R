library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)  # for fct_reorder

file_path <- "/Users/spencertoth/Desktop/biocode/barrnap_phyloflash_comparison/Moorea_Corals/CORAL1286_M997_SSU_RRNA - 16S_barrnap_phyloflash_mapping.csv"

# Read CSV
data <- read.csv(file_path, stringsAsFactors = FALSE)

# Filter for pident between 99 and 100
data_100 <- data %>% filter(between(pident, 99, 100))

# Filter for contigs with query length (qlen) >= 1000
data_100_long <- data_100 %>% filter(qlen >= 1400)

# Calculate percent overlap = (alignment length / qlen) * 100
data_100_long <- data_100_long %>%
  mutate(percent_overlap = (length / qlen) * 100)

# Categorize lengths based on 'length' column around 1000 ± 10 bp
data_100_long <- data_100_long %>%
  mutate(length_category = case_when(
    length < (1400 - 10) ~ "less than 1000 bp",
    length >= (1400 - 10) & length <= (1000 + 10) ~ "1000 bp (±10)",
    length > (1400 + 10) ~ "greater than 1000 bp"
  ))

# Summarize counts and average percent overlap per Sample_ID and length_category
summary_df <- data_100_long %>%
  group_by(Sample_ID, length_category) %>%
  summarize(
    megahit_contig_count = n_distinct(qseqid),
    avg_percent_overlap = mean(percent_overlap, na.rm = TRUE),
    .groups = "drop"
  )

# Prepare all sample IDs and categories for complete framing
all_samples <- unique(summary_df$Sample_ID)
all_categories <- c("less than 1000 bp", "1000 bp (±10)", "greater than 1000 bp")

# Complete data frame to include zero counts for categories missing in summary
complete_df <- expand.grid(Sample_ID = all_samples, length_category = all_categories,
                           stringsAsFactors = FALSE) %>%
  left_join(summary_df, by = c("Sample_ID", "length_category")) %>%
  mutate(
    megahit_contig_count = ifelse(is.na(megahit_contig_count), 0, megahit_contig_count),
    avg_percent_overlap = ifelse(is.na(avg_percent_overlap), NA, avg_percent_overlap)
  )

# Factor for ordering in plot
complete_df$Sample_ID <- factor(complete_df$Sample_ID, levels = all_samples)
complete_df$length_category <- factor(complete_df$length_category, levels = all_categories)

# Compute total counts per sample for ordering
sample_totals <- complete_df %>%
  group_by(Sample_ID) %>%
  summarise(total_count = sum(megahit_contig_count, na.rm = TRUE), .groups = "drop")

# Reorder factors by total counts ascending (lowest top)
complete_df <- complete_df %>%
  left_join(sample_totals, by = "Sample_ID") %>%
  mutate(Sample_ID = fct_reorder(Sample_ID, total_count, .desc = FALSE))

max_count <- max(complete_df$megahit_contig_count, na.rm = TRUE)
# Example approach: for each Sample_ID and length_category, assign avg_percent_overlap as fill

ggplot(filter(complete_df, megahit_contig_count > 0), 
       aes(x = megahit_contig_count, y = Sample_ID, fill = avg_percent_overlap)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  
  geom_text(aes(label = megahit_contig_count),
            position = position_stack(vjust = 0.5),
            size = 3,
            color = "black") +
  
  labs(
    title = "MEGAHIT-Barrnap 16S SSU rRNA contigs with ≥99% identity to Phyloflash contigs (qlen ≥ 1400)",
    y = "Coral Samples (sorted by total sequence count, highest top)",
    x = "Number of MEGA-HIT contigs",
    fill = "Percent Query Coverage"
  ) +
  
  # Color gradient from low to high coverage
  scale_fill_gradient(low = "blue", high = "red", na.value = "grey80") +
  
  scale_x_continuous(
    breaks = seq(0, max_count + 5, by = 1),
    limits = c(0, max_count + 5)
  ) +
  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0))
