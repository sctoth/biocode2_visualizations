library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)

# Load your CSV file (adjust the path as needed)
dat <- read.csv("/Users/spencertoth/Desktop/biocode/barrnap/M997_CORAL1286/28S/sequence_counts.csv", stringsAsFactors = FALSE)

# Gather all length columns into long format
dat_long <- dat %>%
  pivot_longer(
    cols = starts_with("Length_"),
    names_to = "Sequence",
    values_to = "Length"
  ) %>%
  filter(!is.na(Length) & Length != "")  # Remove empty sequence slots

# Ensure sample_order is unique and sorted descending by Count
sample_order <- dat %>%
  distinct(Samples, Count) %>%       # Remove duplicates if any
  arrange((Count)) %>%           # Sort descending by Count (highest first)
  pull(Samples)

# Apply factor levels to dat_long$Samples for ordering in plot
dat_long$Samples <- factor(dat_long$Samples, levels = sample_order)

ggplot(dat_long, aes(x = Length, y = Samples, fill = Sequence)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  labs(
    title = "Lengths of 28S Sequences per Sample",
    y = "Sample",
    x = "Sequence Length (bp)",
    fill = "Sequence"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0),
    axis.text.y = element_text(size = 5)  # Decrease y-axis font size here
  )
