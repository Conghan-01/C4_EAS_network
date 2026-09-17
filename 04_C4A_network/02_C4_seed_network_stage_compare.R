library(tidyverse)
setwd("/gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/02_C4_network/downsampling")
pre_neg <- read.csv(
"C4A_Network_tieredFullSet_Prenatal_R0.3_Negative.csv"
)

post_neg <- read.csv(
"C4A_Network_tieredFullSet_Postnatal_R0.3_Negative.csv"
)

pre_genes <- gsub("\\..*", "", pre_neg$Gene)

post_genes <- gsub("\\..*", "", post_neg$Gene)


overlap <- intersect(pre_genes, post_genes)

length(overlap)


pre <- read.csv(
"C4A_Network_Prenatal_N170.csv"
)

post <- read.csv(
"C4A_Resampled_Median_Network_Postnatal_N170.csv"
)


pre$Gene <- gsub("\\..*","",pre$Gene)
post$Gene <- gsub("\\..*","",post$Gene)


common <- intersect(
pre_genes,
post_genes
)


df <- inner_join(
pre,
post,
by="Gene"
)


df_shared <- df %>% 
filter(Gene %in% common)


cor.test(
df_shared$R,
df_shared$Median_R
)