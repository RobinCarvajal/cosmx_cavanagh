
# Load pseudo_brain -------------------------------------------------------

pseudo_brain <- readRDS(file.path(data_dir,"single_cell","differential_expression",
                                  "broad.labels_pseudo_bulk.RDS"))


# RANDOM EM table ---------------------------------------------------------------

# just take any EM
em_path <- file.path(data_dir,"single_cell","differential_expression",
                     "condition","Astrocytes","em_Astrocytes.tsv")
counts <- read.table(em_path, header = TRUE, sep = "\t", row.names = 1)

# segment meaning ---------------------------------------------------------

# segment 1: mouse (e.g. s1055)
# segment 2: position (e.g. ant)
# segment 3: brain slice (e.g. bs1)
# segment 4: condition (e.g. A)

# COMPLETE SS -------------------------------------------------------------

# generate sample level metadata
ss <- tibble(row.names = colnames(counts))
names(ss)[1] <- "sample"
print(ss)

segments <- c(4) # condition

replacement <- paste0("\\", segments[1])

# modified so that Searchligth2 will accept it
ss$sample_group <- gsub('^([^_]+)_([^_]+)+_([^_]+)+_([^_]+)$',
                             replacement, 
                             colnames(counts))

# print to check
print(ss)

# save path
out<-file.path(data_dir, "single_cell","differential_expression", "condition", "ss")
if (!dir.exists(out)) {dir.create(out, recursive = TRUE)}
# save
write.table(ss,
            file = file.path(out, "ss_universal.tsv"),
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")

