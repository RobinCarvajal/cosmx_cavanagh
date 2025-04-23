
# Load pseudo_brain -------------------------------------------------------

pseudo_brain <- readRDS(file.path(data_dir,"single_cell","differential_expression",
                                  "broad.labels_pseudo_bulk.RDS"))


# RANDOM EM table ---------------------------------------------------------------

# just take any EM
em_path <- file.path(data_dir,"single_cell","differential_expression",
                     "cell_types","Astrocytes","em_Astrocytes.tsv")
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

segments <- c(2, # position
              4) # condition

replacement <- paste0("\\", segments[1], "\\", segments[2])

# modified so that Searchligth2 will accept it
ss$sample_group <- gsub('^([^_]+)_([^_]+)+_([^_]+)+_([^_]+)$',
                        replacement, 
                        colnames(counts))

# print to check
print(ss)


# save path
out<-file.path(data_dir, "single_cell","differential_expression", "position.condition", "ss")
if (!dir.exists(out)) {dir.create(out, recursive = TRUE)}
# save
write.table(ss,
            file = file.path(out, "ss_universal.tsv"),
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")


# antC vs posC ------------------------------------------------------------

# just make a subset of the ss
ss_antC_vs_posC <- subset(ss, sample_group %in% c("antC","posC"))

# write csv to out location
write.table(ss_antC_vs_posC,
          file = file.path(out, "ss_antC_vs_posC.csv"),
          row.names = FALSE, quote = FALSE, sep = "\t")


# antA vs posA ------------------------------------------------------------

# just make a subset of the ss
ss_antA_vs_posA <- subset(ss, sample_group %in% c("antA","posA"))

# write csv to out location
write.table(ss_antA_vs_posA,
          file = file.path(out, "ss_antA_vs_posA.csv"),
          row.names = FALSE, quote = FALSE, sep = "\t")

# antC vs antA ------------------------------------------------------------

# just make a subset of the ss
ss_antC_vs_antA <- subset(ss, sample_group %in% c("antC","antA"))

# write csv to out location
write.table(ss_antC_vs_antA,
          file = file.path(out, "ss_antC_vs_antA.csv"),
          row.names = FALSE, quote = FALSE, sep = "\t")

# posC vs posA ------------------------------------------------------------

# just make a subset of the ss
ss_posC_vs_posA <- subset(ss, sample_group %in% c("posC","posA"))

# write csv to out location
write.table(ss_posC_vs_posA,
          file = file.path(out, "ss_posC_vs_posA.csv"),
          row.names = FALSE, quote = FALSE, sep = "\t")


