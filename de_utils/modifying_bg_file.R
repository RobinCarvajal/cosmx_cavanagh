bg <- read.csv("F:/cosmx_data/single_cell/differential_expression/mouse.GRCm39_original.csv", sep="\t")

# change the order of the columns 
bg <- bg[,c(2,1,3,4,5,6)]

# rename the columns
names(bg)[1] <- "ID"
names(bg)[2] <- "SYMBOL"

bg$SYMBOL <- bg$ID

# remove duplicates in the first column 
bg <- bg[!duplicated(bg$ID),]

# save 
write.table(bg, file = "F:/cosmx_data/single_cell/differential_expression/mouse.GRCm39.csv", sep = "\t", row.names = FALSE, quote = FALSE)
