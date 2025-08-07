library(ggplot2)
library(dplyr)
library(forcats)

# Path to your input file
file_path <- "/Users/spencertoth/Desktop/biocode/barrnap_phyloflash_comparison/Moorea_Corals/CORAL1286_M997_SSU_RRNA - bac_barrnap_assemblies_all.csv"

# Read CSV
data <- read.csv(file_path, stringsAsFactors = FALSE)

# Select only the columns you need
data <- data %>% dplyr::select(superkingdom, Sample_ID, query_id, multi)

# Create taxonomy_full as superkingdom (adjust if you want to combine more taxonomy levels)
data <- data %>%
  mutate(taxonomy_full = superkingdom)

# Custom colors for taxa (adjust as needed)
custom_colors <- c(
  "Bacteria | Actinobacteriota" = "#9e37c8",
  "Bacteria | Bacteroidota" = "#5f1796",
  "Bacteria | Campylobacterota" = "#a445c7",
  "Bacteria | Cyanobacteria" = "#d61a9a",
  "Bacteria | Firmicutes" = "#ad9efa",
  "Bacteria | Patescibacteria" = "#7133ba",
  "Bacteria | Proteobacteria" = "#e27bdc",
  "Bacteria | Verrucomicrobiota" = "#d5b0e9",
  "Eukaryota | Other" = "#2ca730",
  "Eukaryota | Cnidaria" = "#7fff00",
  "Eukaryota | Apicomplexa sp." = "#53c06a",
  "Eukaryota | Symbiodinium sp." = "#3399ff",
  "Bacteria" = "purple",
  "Eukaryota" = "green"
)

# Order samples by total coverage descending
sample_order <- data %>%
  group_by(Sample_ID) %>%
  summarise(total_coverage = sum(multi, na.rm = TRUE)) %>%
  arrange(desc(total_coverage)) %>%
  pull(Sample_ID)

data$Sample_ID <- factor(data$Sample_ID, levels = sample_order)

# Order taxonomy_full by abundance descending
taxa_order <- data %>%
  group_by(taxonomy_full) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  pull(taxonomy_full)

data$taxonomy_full <- factor(data$taxonomy_full, levels = taxa_order)

# Order contigs within samples by coverage descending
data <- data %>%
  group_by(Sample_ID) %>%
  arrange(desc(multi), .by_group = TRUE) %>%
  mutate(contig_id = factor(query_id, levels = unique(query_id))) %>%
  ungroup()

# Plot: coverage (multi) per contig, faceted by sample
ggplot(data, aes(x = contig_id, y = multi, fill = taxonomy_full)) +
  geom_col(color = "black") +
  facet_wrap(~ Sample_ID, scales = "free_x") +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = "MEGA-HIT Barrnap: SSU rRNA Contig Coverage per Sample (ordered by total multi)",
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
