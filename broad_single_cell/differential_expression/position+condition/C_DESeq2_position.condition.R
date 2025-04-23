# load brain/pseudo_brain objects
pseudo_brain <- readRDS(file.path(data_dir,"single_cell","differential_expression",
                                  "broad.labels_pseudo_bulk.RDS"))

# Getting DE tables -------------------------------------------------------

GetDE_fn <- function(em_path, celltype, suffix, sample_sheet, reference, contrast, p_threshold, fold_threshold, out){
  # em_list: list of em df
  # celltype : celltype name to be used in the output file name
  # sample_sheet: sample sheet for DESeq2
  # suffix: suffix to add to the file name
  # reference <- reference var
  # contrast <- contrast var
  # fold_threshold: default 1 
  # p_threshold: default 0.5
  # out: output directory where the CSVs will be saved

  suppressMessages({
    # CODE
    
    # load em
    em_df <- read.table(em_path, sep="\t", header = T, row.names = 1)
    # set ID column as rownames 
    #rownames(em_df) <- em_df$ID

    
    # subsetting only samples in ss
    counts <- em_df[, sample_sheet$sample]
    
    # filtering 
    counts <- subset(counts,apply(counts, 1, mean) >= 1)
    counts <-  as.matrix(counts)
    
    # other way of filtering
    # keep <- rowSums(counts >=10) >= 3
    # counts <- dds[counts,]
  
    # Deseq2 analysis
    dds <- DESeqDataSetFromMatrix(countData = counts,
                                          colData = sample_sheet,
                                          design = ~ sample_group)

    # setting  reference condition
    dds$sample_group <- relevel(dds$sample_group, ref = reference)
    
    # uncomment in case you want to check the reference and contrast
    print(paste0("Reference: ",levels(dds$sample_group)[1]))
    print(paste0("Contrast: ",levels(dds$sample_group)[2]))
  

  
    # run DESeq2
    dds <- DESeq(dds)
  
  
    # check dispersion shrinkage
    # plotDispEsts(dds, main="Gene Matrix", cex.lab = 1, cex.main=2)
  
    # Generate results object USING LFC
    resLFC <- lfcShrink(dds,
                        contrast = c("sample_group", reference, contrast),
                        type = "ashr", # ashr or apeglm
                        lfcThreshold = fold_threshold,
                        alpha = p_threshold)
  
  
    # converting to a df
    res_df <- data.frame(resLFC)
  
    # order by padj
    res_df <- res_df[order(res_df$padj),]
    # gene(ID) column
    res_df$ID <- rownames(res_df)
    # Move 'ID' column to the first position
    res_df <- res_df[, c("ID", setdiff(names(res_df), "ID"))]
  
    # subset columns
    res_df <- select(res_df, ID, log2FoldChange, pvalue, padj)
    # rename the columns
    colnames(res_df) <- c("ID", "log2Fold", "P", "P.adj")
  
    # write csv to out location
    write.table(res_df,
                file = file.path(out, paste0("de_", celltype, "_", suffix, ".tsv")),
                sep = "\t",
                row.names = FALSE,
                quote = FALSE)
  
    print(paste("DE table for ", celltype," | ", suffix, "saved"))
    
  })

}


# cell types list  --------------------------------------------------------

celltypes <- names(pseudo_brain)


# sample sheets -----------------------------------------------------------

# antC vs posC
antC_vs_posC_ss <- read.table(file.path(data_dir,"single_cell","differential_expression",
                                        "position.condition","ss","ss_antC_vs_posC.csv"), sep = "\t", header = T)
# antA vs posA
antA_vs_posA_ss <- read.table(file.path(data_dir,"single_cell","differential_expression",
                                        "position.condition","ss","ss_antA_vs_posA.csv"), sep="\t", header = T)
# antC vs antA
antC_vs_antA_ss <- read.table(file.path(data_dir,"single_cell","differential_expression",
                                        "position.condition","ss","ss_antC_vs_antA.csv"), sep="\t", header = T)
# pos C vs pos A
posC_vs_posA_ss <- read.table(file.path(data_dir,"single_cell","differential_expression",
                                        "position.condition","ss","ss_posC_vs_posA.csv"), sep="\t", header = T)

# ALL comparisons  -------------------------------------------------------

# NOTE: I do all comparisons in a loop to avoid repeating code

# in a loop
for (celltype in celltypes){
  
  # current celltype
  print(celltype)
  
  # base path
  base_path <- file.path(data_dir,"single_cell","differential_expression", "position.condition", "cell_types")
  
  # out path 
  celltype_path <- file.path(base_path, celltype)
  # make dir
  if (!dir.exists(celltype_path)) {dir.create(celltype_path, recursive = TRUE)}
  
  # em path 
  em <- file.path(base_path, celltype, paste0("em_", celltype, ".tsv"))
  
  ## position 
  # antC vs posC
  GetDE_fn(em_path = em,
           celltype = celltype,
           suffix = "antC_vs_posC", 
           sample_sheet = antC_vs_posC_ss,
           reference = "antC", 
           contrast = "posC", 
           p_threshold = 0.5, 
           fold_threshold = 1, 
           out = celltype_path )
  
  # antA vs posA
  GetDE_fn(em_path = em, 
           celltype = celltype,
           suffix = "antA_vs_posA", 
           sample_sheet = antA_vs_posA_ss, 
           reference = "antA", 
           contrast = "posA", 
           p_threshold = 0.5, 
           fold_threshold = 1, 
           out = celltype_path)

  
  ## condition
  # antC vs antA
  GetDE_fn(em_path = em, 
           celltype = celltype,
           suffix = "antC_vs_antA", 
           sample_sheet = antC_vs_antA_ss, 
           reference = "antC", 
           contrast = "antA", 
           p_threshold = 0.5, 
           fold_threshold = 1, 
           out = celltype_path)
  
  # posC vs posA
  GetDE_fn(em_path = em, 
           celltype = celltype,
           suffix = "posC_vs_posA", 
           sample_sheet = posC_vs_posA_ss, 
           reference = "posC", 
           contrast = "posA", 
           p_threshold = 0.5, 
           fold_threshold = 1, 
           out = celltype_path)
  
}

