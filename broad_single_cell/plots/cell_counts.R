# load the brain object 
brain <- readRDS(file.path(objects_dir,"brain_objects",
                           "2024_11_05_brain_v1.5_broad.labels_added.RDS"))

# get a table of the sample - broad.label - condition
meta <- brain@meta.data

# get the cell type proportions within each sample
counts_summary <- meta %>%
  select(broad.labels, sample_name, condition) %>%
  group_by(broad.labels, sample_name, condition) %>%
  summarize(count = n(), .groups = "drop")

counts_summary <- meta %>%
  select(broad.labels, sample_name, condition) %>%
  group_by(broad.labels, sample_name, condition) %>%
  summarize(count = n(), .groups = "drop") %>%
  group_by(sample_name, condition) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  ungroup()


# View the result
counts_summary <- as.data.frame(counts_summary)

#  complex plot -----------------------------------------------------------

# Install ggsignif if not already installed
if (!requireNamespace("ggsignif", quietly = TRUE)) {
  install.packages("ggsignif")
}

# Calculate p-values for each cluster
library(dplyr)
p_values <- counts_summary %>%
  group_by(broad.labels) %>%
  summarize(p_value = t.test(proportion ~ condition)$p.value)

# Add significance labels
p_values <- p_values %>%
  mutate(significance = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ "ns"
  ))

library(ggsignif)

# Merge significance data with plot data
counts_summary <- counts_summary %>%
  left_join(p_values, by = "broad.labels")
######
# Dynamically calculate y_position for significance annotation
y_positions <- counts_summary %>%
  group_by(broad.labels) %>%
  summarize(y_max = max(proportion)) %>%
  mutate(y_position = y_max + 5) # Add padding above max y-value

# Merge y_positions into p_values
p_values <- p_values %>%
  left_join(y_positions, by = "broad.labels")

# Plot with significance (proportion)
p <- ggplot(counts_summary, aes(x = broad.labels, y = proportion, fill = condition)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_jitter(color = "black", alpha = 1, position = position_dodge(width = 0.8)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Cell type proportions per sample",
       x = "Cell Type (Broad Labels)",
       y = "Cell Proportion") +
  scale_fill_manual(values = c("blue", "red")) +
  scale_color_manual(values = c("blue", "red")) +
  geom_signif(data = p_values,
              mapping = aes(xmin = broad.labels, xmax = broad.labels, annotations = significance, y_position = y_position),
              manual = TRUE, inherit.aes = FALSE,
              textsize = 4)

print(p)

# save to plots di r
ggsave(file.path(data_dir, "broad_single_cell", "cell_counts.png"), plot = p, width = 10, height = 6, dpi = 300, bg="white")



