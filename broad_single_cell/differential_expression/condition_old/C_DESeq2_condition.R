# load brain/pseudo_brain objects
pseudo_brain <- readRDS(file.path(objects_dir,"pseudo_objects",
                                  "broad.labels_pseudo_bulk.RDS"))

# sample sheets -----------------------------------------------------------

# Control vs Aldara
C_vs_A_ss <- read.table(file.path(data_dir,"broad_single_cell","differential_expression",
                                        "condition","ss","ss_universal.tsv"), sep = "\t", header = T)

# ALL comparisons  -------------------------------------------------------

# NOTE: I do all comparisons in a loop to avoid repeating code

# in a loop
for (celltype in names(pseudo_brain)){
  
  # current celltype
  print(celltype)
  
  # base path
  base_path <- file.path(data_dir,"broad_single_cell","differential_expression","condition","cell_types")
  
  # out path 
  celltype_path <- file.path(base_path, celltype)
  # make dir
  if (!dir.exists(celltype_path)) {dir.create(celltype_path, recursive = TRUE)}
  
  # em path 
  em <- pseudo_brain[[celltype]]
  
  # Control vs Aldara
  GetDE_fn(em_df = em,
           celltype = celltype,
           suffix = "C_vs_A", 
           sample_sheet = C_vs_A_ss,
           reference = "C", 
           contrast = "A",
           fold_threshold = 1,
           p_threshold = 0.05,
           out = celltype_path )
  
}

