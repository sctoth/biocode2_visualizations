library(ggplot2)

library(dplyr)
library(stringr)
# Adjust the path to your actual file location
file_path <- "/Users/spencertoth/Desktop/biocode/barrnap_phyloflash_comparison/Moorea_Corals/CORAL1286_M997_SSU_RRNA - bacterial_barrnap_not_mapped.csv"

# Read the data
data <- read.csv(file_path, stringsAsFactors = FALSE)


# Suppose your data frame is called data and has a column query_ID like "16S_rRNA::k141_20270:71364-72912(+)"
data <- data %>%
  mutate(
    # Extract the coordinates using regex:
    start = as.numeric(str_extract(qseqid, "(?<=:)\\d+(?=-)")),
    end = as.numeric(str_extract(qseqid, "(?<=-)\\d+")),
    length_calc = end - start
  )

# Bin the lengths
data <- data %>%
  mutate(length_bin = cut(length_calc,
                          breaks = c(-Inf, 1000, 2000, 3000, 4000, 5000, 6000),
                          labels = c("<1000", "1000-2000", "2000-3000", "3000-4000", "4000-5000", "5000-6000"),
                          right = FALSE))


# Get unique sequence IDs per bin
len_summary <- data %>%
  filter(!is.na(length_bin)) %>%
  group_by(length_bin) %>%
  summarise(seq_count = n_distinct(qseqid))

# Plot the bar chart
ggplot(len_summary, aes(x = length_bin, y = seq_count)) +
  geom_bar(stat = "identity", fill = "purple", color = "black") +
  geom_text(aes(label = seq_count), vjust = -0.5, size = 3.5) +  # Add counts above bars
  scale_y_continuous(
    breaks = seq(0, max(len_summary$seq_count) + 5, by = 5) # More ticks from 0 to max count + 5 by 1
  ) +
  labs(
    title = "Bacterial 16S rRNA assemblies from MEGA-HIT that did not map to Phyloflash assemblies",
    x = "Assembly Length (bp)",
    y = "Number of MEGA-HIT contigs"
  ) +
  theme_minimal()

