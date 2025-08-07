library(ggplot2)
library(dplyr)
library(forcats)

# Path to your input file
# Path to your input file (update as needed)
#file_path <- "/Users/spencertoth/Desktop/biocode/phyloflash/M997_CORAL1286/phyloflash_M997_CORAL1286_assemblies_revised_06.csv"
#file_path <- "/Users/spencertoth/Desktop/biocode/phyloflash/M997_CORAL1286/phyloflash_M997_CORAL1286_assemblies_revised_euk_01.csv"
file_path <- "/Users/spencertoth/Desktop/biocode/barrnap_phyloflash_comparison/Moorea_Corals/CORAL1286_M997_SSU_RRNA - 16S_phyloflash_assemblies.csv"

# Read CSV
data <- read.csv(file_path, stringsAsFactors = FALSE)


# Combine superkingdom and supergroup for coloring/facetting
data <- data %>%
  mutate(
    taxonomy_full = paste(superkingdom, supergroup, sep = " | "),
    contig_id = OTU  # Assuming each OTU is unique per contig
  )

# If you want to color by taxonomy or by sample, adjust as needed
custom_colors <- c(
  # Bacteria groups (brighter purples and pinks)
  "Bacteria | Actinobacteriota" = "#9e37c8",    # vivid purple
  "Bacteria | Bacteroidota" = "#5f1796",        # bright pink-purple
  "Bacteria | Campylobacterota" = "#a445c7",    # strong medium purple
  "Bacteria | Cyanobacteria" = "#d61a9a",       # bright magenta-pink
  "Bacteria | Firmicutes" = "#ad9efa",          # bright royal purple
  "Bacteria | Patescibacteria" = "#7133ba",     # bright lavender
  "Bacteria | Proteobacteria" = "#e27bdc",      # vivid dark purple
  "Bacteria | Verrucomicrobiota" = "#d5b0e9",   # pale bright purple
  
  # Eukaryota groups (brighter greens and blues)
  "Eukaryota | Other" = "#2ca730",               # bright rich green
  "Eukaryota | Cnidaria" = "#7fff00",            # bright lime green
  
  # SAR groups as Eukaryota supergroup members
  "Eukaryota | Apicomplexa sp." = "#53c06a",     # vibrant softer green
  "Eukaryota | Symbiodinium sp." = "#3399ff"      # bright vivid blue
)


# Order samples by total number of contigs (optional)
sample_order <- data %>%
  group_by(Sample) %>%
  summarise(contig_count = n()) %>%
  arrange(desc(contig_count)) %>%
  pull(Sample)

data$Sample <- factor(data$Sample, levels = sample_order)

# Order taxonomy_full by abundance (optional, for color legend)
taxa_order <- data %>%
  group_by(taxonomy_full) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  pull(taxonomy_full)
data$taxonomy_full <- factor(data$taxonomy_full, levels = taxa_order)

# Order contigs within each sample by coverage, descending (optional)
data <- data %>%
  group_by(Sample) %>%
  arrange(desc(coverage), .by_group = TRUE) %>%
  mutate(contig_id = factor(contig_id, levels = unique(contig_id))) %>%
  ungroup()
# Calculate total coverage per Sample
sample_totals <- data %>%
  group_by(Sample) %>%
  summarise(total_coverage = mean(coverage, na.rm = TRUE)) %>%
  arrange(desc(total_coverage))

# Reorder Sample factor levels by descending total coverage
data$Sample <- factor(data$Sample, levels = sample_totals$Sample)


# Now your plotting code uses the reordered Sample factor
ggplot(data, aes(x = contig_id, y = coverage, fill = taxonomy_full)) +
  geom_col(color = "black") +
  facet_wrap(~ Sample, scales = "free_x") +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = "SSU rRNA Contig Coverage per Sample (ordered by total coverage per taxa)",
    x = "Contigs (per Sample)",
    y = "Coverage",
    fill = "Taxonomy"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 4),
    strip.text.x = element_text(size = 5),
    legend.position = "right"
  )
