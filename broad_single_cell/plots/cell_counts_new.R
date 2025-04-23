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

# export counts summary 
write.csv(counts_summary, file = file.path(data_dir, "broad_single_cell", "counts_summary_per_sample.csv"), row.names = FALSE)