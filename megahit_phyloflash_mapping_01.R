library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)  # for fct_reorder

file_path <- "/Users/spencertoth/Desktop/biocode/barrnap_phyloflash_comparison/Moorea_Corals/CORAL1286_M997_SSU_RRNA - 16S_barrnap_phyloflash_mapping.csv"

data <- read.csv(file_path, stringsAsFactors = FALSE)

data_100 <- data %>% filter(between(pident, 95, 100))

data_100 <- data_100 %>%
  mutate(length_category = case_when(
    qlen < (slen - 10) ~ "Shorter",
    qlen >= (slen - 10) & qlen <= (slen + 10) ~ "Equal (±10 bp)",
    qlen > (slen + 10) ~ "Longer"
  ))

summary_df <- data_100 %>%
  group_by(Sample_ID, length_category) %>%
  summarize(megahit_contig_count = n_distinct(qseqid), .groups = "drop")

all_samples <- unique(summary_df$Sample_ID)
all_categories <- c("Shorter", "Equal (±10 bp)", "Longer")

complete_df <- expand.grid(Sample_ID = all_samples, length_category = all_categories,
                           stringsAsFactors = FALSE) %>%
  left_join(summary_df, by = c("Sample_ID", "length_category")) %>%
  mutate(megahit_contig_count = ifelse(is.na(megahit_contig_count), 0, megahit_contig_count))

complete_df$Sample_ID <- factor(complete_df$Sample_ID, levels = all_samples)
complete_df$length_category <- factor(complete_df$length_category, levels = all_categories)

# Prepare zero bars for black lines
zero_bars <- complete_df %>% filter(megahit_contig_count == 0)

dodge_width <- 0.9
n_groups <- length(all_categories)

# Position adjustments for zero bars and rectangles
zero_bars <- zero_bars %>%
  mutate(
    sample_num = as.numeric(Sample_ID),
    category_num = as.numeric(length_category),
    x_pos = sample_num - dodge_width/2 + (category_num - 0.5) * (dodge_width / n_groups)
  )

# Compute xmin and xmax for sample boxes (a little padding around group)
rect_df <- data.frame(
  sample_num = seq_along(all_samples),
  xmin = seq_along(all_samples),
  xmax = seq_along(all_samples) + 0.5
)
max_count <- max(complete_df$megahit_contig_count, na.rm = TRUE)

# Compute total counts per sample
sample_totals <- complete_df %>%
  group_by(Sample_ID) %>%
  summarise(total_count = sum(megahit_contig_count, na.rm = TRUE))

# Reorder Sample_ID factor by descending total_count
complete_df <- complete_df %>%
  left_join(sample_totals, by = "Sample_ID") %>%
  mutate(Sample_ID = fct_reorder(Sample_ID, total_count, .desc = FALSE))

max_count <- max(complete_df$megahit_contig_count, na.rm = TRUE)

# Check for NA in megahit_contig_count, just in case
na_rows <- complete_df[is.na(complete_df$megahit_contig_count), ]

# Check values outside the x-axis limits
x_limits <- c(0, max(complete_df$megahit_contig_count, na.rm = TRUE) + 1)
outside_limit_rows <- complete_df[
  complete_df$megahit_contig_count < x_limits[1] | complete_df$megahit_contig_count > x_limits[2],
]

ggplot(filter(complete_df, megahit_contig_count > 0), aes(x = megahit_contig_count, y = Sample_ID, fill = length_category)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  
  # Add segment labels centered vertically within each stacked bar segment
  geom_text(aes(label = megahit_contig_count),
            position = position_stack(vjust = 0.5),
            size = 3,     # adjust text size as needed
            color = "black") +
  
  labs(
    title = "MEGAHIT-Barrnap 16S SSU rRNA contigs with ≥95% identity to Phyloflash contigs",
    y = "Coral Samples (sorted by total sequence count, highest top)",
    x = "Number of MEGA-HIT contigs",
    fill = "Length category"
  ) +
  scale_x_continuous(
    breaks = seq(0, max_count + 3, by = 1),
    limits = c(0, max_count + 3)
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0))

# --- ADD THIS TO THE END OF YOUR SCRIPT ---

# Arrange sample_totals in descending order
sample_totals_sorted <- sample_totals %>%
  arrange(desc(total_count))

# Join total_count info into summary_table and sort descending by total
summary_table_sorted <- summary_table %>%
  left_join(sample_totals, by = "Sample_ID") %>%
  arrange(desc(total_count)) %>%
  select(-total_count)  # Remove extra column if you only want categories

cat("\nSummary Table (sorted by total sequences, highest to lowest):\n")
print(summary_table_sorted, n = Inf)

cat("\nTotal number of contigs with ≥99% identity per sample (sorted):\n")
print(sample_totals_sorted, n = Inf)


