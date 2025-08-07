library(ggplot2)
library(dplyr)
library(forcats)

# Path to your input file (update as needed)
file_path <- "/Users/spencertoth/Desktop/biocode/phyloflash/M997_CORAL1286/CORAL1286_M997_SSU_RRNA - 16S_bac_barrnap_assemblies.csv"
#file_path <- "/Users/spencertoth/Desktop/biocode/barrnap_phyloflash_comparison/Moorea_Corals/CORAL1286_M997_SSU_RRNA - bac_barrnap_assemblies_all.csv"
# Read CSV
data <- read.csv(file_path, stringsAsFactors = FALSE)

# Summarize counts of unique query_ids by Sample and taxonomy levels
summary_df <- data %>%
  group_by(Sample_ID, superkingdom) %>%
  summarise(contig_count = n_distinct(query_id), .groups = "drop") %>%
  # Combine taxonomy levels for fill assignment
  mutate(taxonomy_full = paste(superkingdom))

custom_colors <- c(
  # Bacteria groups (brighter purples and pinks)
  "Bacteria" = "purple",
  "Bacteria | Actinobacteriota" = "#9e37c8",    # vivid purple
  "Bacteria | Bacteroidota" = "#5f1796",        # bright pink-purple
  "Bacteria | Campylobacterota" = "#a445c7",    # strong medium purple
  "Bacteria | Cyanobacteria" = "#d61a9a",       # bright magenta-pink
  "Bacteria | Firmicutes" = "#ad9efa",          # bright royal purple
  "Bacteria | Patescibacteria" = "#7133ba",     # bright lavender
  "Bacteria | Proteobacteria" = "#e27bdc",      # vivid dark purple
  "Bacteria | Verrucomicrobiota" = "#d5b0e9",   # pale bright purple
  
  # Eukaryota groups (brighter greens and blues)
  "Eukaryota" = "green",
  "Eukaryota | Other" = "#2ca730",               # bright rich green
  "Eukaryota | Cnidaria" = "#7fff00",            # bright lime green
  
  # SAR groups as Eukaryota supergroup members
  "Eukaryota | Apicomplexa sp." = "#53c06a",     # vibrant softer green
  "Eukaryota | Symbiodinium sp." = "#3399ff"      # bright vivid blue
)



# Assign a default color for any missing taxonomy groups
missing_taxa <- setdiff(unique(summary_df$taxonomy_full), names(custom_colors))
if(length(missing_taxa) > 0) {
  custom_colors[missing_taxa] <- "#CCCCCC"  # light grey for unknown
}

# Calculate total contigs per taxonomy_full and superkingdom for ordering legend
taxa_totals <- summary_df %>%
  group_by(taxonomy_full, superkingdom) %>%
  summarise(total = sum(contig_count), .groups = "drop")

# Separate Bacteria and Eukaryota and order by descending total count
bacteria_levels <- taxa_totals %>%
  filter(superkingdom == "Bacteria") %>%
  arrange(desc(total)) %>%
  pull(taxonomy_full)

eukaryota_levels <- taxa_totals %>%
  filter(superkingdom == "Eukaryota") %>%
  arrange(desc(total)) %>%
  pull(taxonomy_full)

# Combine factor levels with Bacteria first then Eukaryota
ordered_levels <- c(bacteria_levels, eukaryota_levels)

# Set taxonomy_full factor with ordered levels for legend and plotting
summary_df$taxonomy_full <- factor(summary_df$taxonomy_full, levels = ordered_levels)

# Set priority order of taxa within stacks
priority_order <- c(
  "Eukaryota | Symbiodinium sp.",
  "Eukaryota | Cnidaria"
)

all_taxa <- if(is.factor(summary_df$taxonomy_full)) levels(summary_df$taxonomy_full) else unique(summary_df$taxonomy_full)
other_taxa <- setdiff(all_taxa, priority_order)
new_levels <- c(sort(other_taxa), priority_order)

summary_df$taxonomy_full <- factor(summary_df$taxonomy_full, levels = new_levels)

# Order samples by total contig counts (ascending for vertical axis)
sample_totals <- summary_df %>%
  group_by(Sample_ID) %>%
  summarise(total_contigs = sum(contig_count), .groups = "drop") %>%
  arrange(total_contigs)

summary_df$Sample_ID <- factor(summary_df$Sample_ID, levels = sample_totals$Sample_ID)

# Calculate total contig counts per taxonomy_full (once, without changing factor)
taxa_totals <- summary_df %>%
  group_by(taxonomy_full) %>%
  summarise(total = sum(contig_count), .groups = "drop") %>%
  arrange(desc(total))

# Define the desired legend order by taxa abundance
ordered_taxa <- taxa_totals$taxonomy_full

# Plot
ggplot(summary_df, aes(x = contig_count, y = Sample_ID, fill = taxonomy_full)) +
  geom_bar(stat = "identity", position = "stack", color = "black") + 
  scale_fill_manual(
    values = custom_colors,
    breaks = ordered_taxa,   # only affects legend and stacking order in plot
    name = "BLAST Taxonomy (Superkingdom)"
  ) +
  scale_x_continuous(
    breaks = seq(0, max(sample_totals$total_contigs), by = 5),
    limits = c(0, max(sample_totals$total_contigs) + 1)
  ) +
  labs(
    title = "Moorea Coral Samples: MEGA-HIT Barrnap 16S SSU rRNA Contigs",
    x = "Number of Contigs",
    y = "Moorea Coral Samples"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 5),
    legend.position = "right"
  )

# Identify samples with the most bacterial contigs
bacteria_summary <- summary_df %>%
  filter(grepl("^Bacteria", taxonomy_full)) %>%
  group_by(Sample_ID) %>%
  summarise(total_bacteria_contigs = sum(contig_count), .groups = "drop") %>%
  arrange(desc(total_bacteria_contigs))

# Print top samples with most bacterial contigs
cat("Samples with the most bacterial contigs:\n")
print(bacteria_summary)

# Identify samples whose contigs contain no Cnidaria contigs
# First, find samples with any Cnidaria contigs
samples_with_cnidaria <- summary_df %>%
  filter(grepl("Eukaryota | Cnidaria", taxonomy_full) & contig_count > 0) %>%
  distinct(Sample_ID) %>%
  pull(Sample_ID)

# Samples without Cnidaria contigs are those not in samples_with_cnidaria
samples_without_cnidaria <- setdiff(levels(summary_df$Sample_ID), samples_with_cnidaria)

cat("\nSamples with NO Cnidaria contigs:\n")
print(samples_without_cnidaria)

# For samples without Cnidaria, list all their contigs
# Assuming 'data' dataframe contains query_id, Sample, and taxonomy info for each contig
if ("query_id" %in% colnames(data)) {
  no_cnidaria_contigs <- data %>%
    filter(Sample_ID %in% samples_without_cnidaria) %>%
    select(Sample_ID, query_id) %>%
    distinct()
  
  cat("\nContigs in samples without any Cnidaria contigs:\n")
  print(no_cnidaria_contigs)
} else {
  cat("\nNo query_id or contig-level data available to list contigs without Cnidaria.\n")
}
# Find samples that have any Bacteria contigs
# Identify samples with any Bacteria contigs
samples_with_bacteria <- summary_df %>%
  filter(superkingdom == "Bacteria" & contig_count > 0) %>%
  distinct(Sample_ID) %>%
  pull(Sample_ID)

# Identify all samples in the dataset
all_samples <- unique(summary_df$Sample_ID)

# Find samples with NO Bacteria contigs
samples_without_bacteria <- setdiff(all_samples, samples_with_bacteria)
cat(samples_without_bacteria, sep = "\n")
cat("\nSamples WITHOUT any Bacteria sequences:\n")
print(samples_without_bacteria)
