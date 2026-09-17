# Prepare your gene sets text file in R first:
# "Pos_Network GENE1 GENE2 GENE3..."
# "Neg_Network GENE1 GENE2 GENE3..."

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)

prepare_magma_geneset <- function(file_path, set_name) {
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  ensembl_clean <- gsub("\\..*", "", df$Gene)
  suppressMessages({
    mapped <- tryCatch(
      bitr(ensembl_clean, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
      error = function(e) return(data.frame(ENTREZID = character()))
    )
  })
  
  entrez_ids <- unique(mapped$ENTREZID)
  
  if(length(entrez_ids) > 0) {
    magma_line <- paste(c(set_name, entrez_ids), collapse = " ")
    return(magma_line)
  } else {
    return(NULL)
  }
}


line_pos_pre <- prepare_magma_geneset("../C4A_Positive_Network_pre03.csv", "Prenatal_Positive")
line_neg_pre <- prepare_magma_geneset("../C4A_Negative_Network_pre03.csv", "Prenatal_Negative")
line_pos_post <- prepare_magma_geneset("../C4A_Positive_Network_post03.csv", "Postnatal_Positive")
line_neg_post <- prepare_magma_geneset("../C4A_Negative_Network_post03.csv", "Postnatal_Negative")
line_neg_post <- prepare_magma_geneset("../C4A_Negative_Network_post03.csv", "Postnatal_Negative")
downsampling_neg_post <- prepare_magma_geneset("../Postnatal_C4A_Reproducible_N170_Negative.csv", "Postnatal_Negative_D")
downsampling_pos_post <- prepare_magma_geneset("../Postnatal_C4A_Reproducible_N170_Positive.csv", "Postnatal_Positive_D")



lines_to_write <- c(line_pos_pre, line_neg_pre, line_pos_post, line_neg_post)
lines_to_write <- lines_to_write[!sapply(lines_to_write, is.null)]

#lines_to_write <- c(downsampling_neg_post,downsampling_pos_post)
lines_to_write <- lines_to_write[!sapply(lines_to_write, is.null)]

writeLines(lines_to_write, "C4A_MAGMA_0.3RGeneSets.txt")

