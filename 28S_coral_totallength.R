library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)

# Load CSV
dat <- read.csv("/Users/spencertoth/Desktop/biocode/barrnap/M997_CORAL1286/28S/sequence_counts.csv", stringsAsFactors = FALSE)

# Reshape to long format
dat_long <- dat %>%
  pivot_longer(
    cols = starts_with("Length_"),
    names_to = "Sequence",
    values_to = "Length"
  ) %>%
  filter(!is.na(Length) & Length != "") %>%  # Remove empty values
  mutate(Length = as.numeric(Length))

# Order samples by total length descending and reverse factor for plotting top-to-bottom
# Order samples by total length descending for left-to-right plot
sample_total_length <- dat_long %>%
  group_by(Samples) %>%
  summarise(total_length = sum(Length, na.rm=TRUE)) %>%
  arrange((total_length))
sample_order <- sample_total_length$Samples
dat_long$Samples <- factor(dat_long$Samples, levels = sample_order)

# Order sequences within samples descending by length for stacking order
dat_long <- dat_long %>%
  group_by(Samples) %>%
  mutate(Sequence = factor(Sequence, levels = unique(Sequence[order(-Length)]))) %>%
  ungroup()



# Define your color palette and repeat to cover all sequences
base_cb_colors <- c(
  "red",
  "orange",
  "yellow",
  "green",
  "blue",
  "purple",
  "lightpink"
)
n_seq <- length(unique(dat_long$Sequence))
cb_colors <- rep(base_cb_colors, length.out = length(unique(dat_long$Sequence)))

# Plot with reversed stacking direction to match factor levels
ggplot(dat_long, aes(x = Length, y = Samples, fill = Sequence)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE), color = "black") +
  scale_fill_manual(values = cb_colors) +
  scale_x_continuous(n.breaks = 15) +   # Increase number of x-axis ticks here
  labs(
    title = "Moorea Coral Samples: MEGAHIT-Barrnap 28S SSU rRNA Contigs",
    y = "Moorea Coral Samples (M997 and CORAL1286)",
    x = "Sequence Length (bp)",
    fill = "Sequence"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0),
    axis.text.y = element_text(size = 5)
  )
