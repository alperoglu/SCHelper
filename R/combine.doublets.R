library(tidyverse)

args = commandArgs(TRUE)

lib_file <- args[1]

ids <- read.table(file = lib_file, row.names = NULL, col.names=c("ids", "path"), sep=",")
ids <- ids$ids
for(i in ids){
  # df_doublets <- read.table(paste0("doublets_df_", i, ".txt"))
  # # df_doublets <- gsub(pattern="-1", replacement = paste0("-",i), df_doublets)
  # df_doublets$Sample <- i
  # write.table(df_doublets, file = "df.combined.doublet.cells.txt", quote = F, col.names= F, row.names = F, append=T, sep = "\t")
  # rm(df_doublets)
  
  scr_doublets <- read.table(paste0("doublets/doublets_scrublet_", i, ".csv"))
  # scr_doublets <- gsub(pattern="-1", replacement = paste0("-",i), scr_doublets)
  scr_doublets$Sample <- i
  write.table(scr_doublets, file = "scrublet.combined.doublet.cells.txt", quote = F, col.names= F, row.names = F, append=T, sep = "\t")
  rm(scr_doublets)
  
}