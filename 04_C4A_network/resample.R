# ==============================================================================
# C4A co-expression network robustness analysis
# Postnatal downsampling to prenatal sample size (N=170)
# 100 random resampling
# Reproducible network defined by frequency >= 70%
# ==============================================================================

rm(list=ls())
options(stringsAsFactors = FALSE)

library(data.table)
library(dplyr)
library(tidyverse)

# Parameters
n_iterations <- 100
target_n <- 170
r_cutoff <- 0.3
fdr_cutoff <- 0.05
frequency_cutoff <- 0.7

expr_file <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/bulk_data/postnatal.logTMM.ComBat.txt"
cov_file <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/bulk_data/postnatalcovariatesToUse.PCAforQTLfromINT.txt"

output_prefix <- "V2Postnatal_C4A_Reproducible_N170"


# Load expression matrix
message("Loading expression matrix...")

expr <- fread(expr_file,header = TRUE,sep = "\t",check.names = FALSE)
gene<- read.table("/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/01_data/RNA-seq_bak/postnatal/genes_to_keep.txt")
gene_ids <- expr[[1]]
expr <- as.matrix(expr[, -1, with = FALSE])
rownames(expr) <- gene_ids
expr<- expr[gene$V1,]

# Load covariates
message("Loading covariates...")

cov <- fread(cov_file,header = TRUE,check.names = FALSE)

sample_ids <- cov[[1]]
cov <- as.matrix(cov[, -1, with = FALSE])
rownames(cov) <- sample_ids
cov <- as.data.frame(t(cov))


# Match samples
common_samples <- intersect(colnames(expr),rownames(cov))

expr <- expr[, common_samples]
cov <- cov[common_samples, ]

message("Samples available: ", length(common_samples))


# Identify C4A
c4a_id <- grep("ENSG00000244731",rownames(expr),value = TRUE)[1]
message("C4A ID: ", c4a_id)

# Storage
all_R <- matrix(NA, nrow = nrow(expr),ncol = n_iterations,dimnames = list(rownames(expr), NULL))
all_FDR <- matrix(NA,nrow = nrow(expr),ncol = n_iterations,dimnames = list(rownames(expr), NULL))
all_negative <- matrix(FALSE,nrow = nrow(expr),ncol = n_iterations,dimnames = list(rownames(expr), NULL))
all_positive <- matrix(FALSE,nrow = nrow(expr), ncol = n_iterations,dimnames = list(rownames(expr), NULL))

# Storage - 增加用于记录每次迭代基因数目的向量
n_iter_positive <- numeric(n_iterations)
n_iter_negative <- numeric(n_iterations)

# Resampling
for(i in 1:n_iterations){

  message("Iteration ", i, "/", n_iterations)

  set.seed(i)

  # Randomly select N=170 samples
  selected_samples <- sample(
    common_samples,
    target_n
  )

  sub_expr <- expr[, selected_samples]
  sub_cov <- cov[selected_samples, ]

  # Remove known covariates
  covs_to_regress <- sub_cov %>%
    select(
      sex,
      RIN,
      Age,
      PC1,
      PC2,
      PC3
    )

  fit <- lm(
    t(sub_expr) ~ .,
    data = covs_to_regress
  )

  expr_res <- t(residuals(fit))
  rownames(expr_res) <- rownames(expr)

  # Pearson correlation with C4A
  c4a_exp <- expr_res[c4a_id, ]

  r <- cor(
    t(expr_res),
    c4a_exp,
    method = "pearson"
  )[,1]

  # P value and FDR
  t_stat <- r * sqrt(
    (target_n - 2) / (1 - r^2)
  )

  p <- 2 * pt(
    -abs(t_stat),
    df = target_n - 2
  )

  fdr <- p.adjust(
    p,
    method = "BH"
  )

  # Store results
  all_R[,i] <- r
  all_FDR[,i] <- fdr

  # Network membership
  all_negative[,i] <- (
    r < -r_cutoff &
    fdr < fdr_cutoff
  )

  all_positive[,i] <- (
    r > r_cutoff &
    fdr < fdr_cutoff
  )
  
  # 统计并记录当前迭代满足条件的基因数目（排除 C4A 自身）
  current_neg_genes <- rownames(expr)[all_negative[,i] & rownames(expr) != gsub("\\..*", "", c4a_id)]
  current_pos_genes <- rownames(expr)[all_positive[,i] & rownames(expr) != gsub("\\..*", "", c4a_id)]
  
  n_iter_negative[i] <- length(current_neg_genes)
  n_iter_positive[i] <- length(current_pos_genes)
}

# 将每轮迭代的统计结果保存为表格，方便后续查看或画图
iteration_stats <- data.frame(
  Iteration = 1:n_iterations,
  Positive_Count = n_iter_positive,
  Negative_Count = n_iter_negative
)

write.csv(
  iteration_stats,
  paste0(output_prefix, "_Iteration_Counts.csv"),
  row.names = FALSE
)

# 打印均值和标准差供参考
cat("\n============================\n")
cat("Iteration counts summary:\n")
cat(sprintf("Positive network genes per iteration: Mean = %.1f, SD = %.1f\n", mean(n_iter_positive), sd(n_iter_positive)))
cat(sprintf("Negative network genes per iteration: Mean = %.1f, SD = %.1f\n", mean(n_iter_negative), sd(n_iter_negative)))
cat("============================\n")

# ==============================================================================
# GO enrichment analysis for reproducible C4A co-expression networks
#
# Postnatal N=170 downsampling robustness analysis
#
# Input:
#   reproducible positive / negative C4A networks
#
# Background:
#   all genes tested in co-expression analysis
# ==============================================================================


rm(list=ls())
options(stringsAsFactors = FALSE)


library(data.table)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(stringr)



# ==============================================================================
# File settings
# ==============================================================================


positive_file <- "Postnatal_C4A_Reproducible_N170_Positive.csv"

negative_file <- "Postnatal_C4A_Reproducible_N170_Negative.csv"


# expression matrix used for network analysis
expr_file <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/bulk_data/postnatal.logTMM.ComBat.txt"



# ==============================================================================
# Load reproducible networks
# ==============================================================================


pos_net <- fread(
  positive_file
)

neg_net <- fread(
  negative_file
)


cat("Positive network genes:",
    nrow(pos_net),
    "\n")

cat("Negative network genes:",
    nrow(neg_net),
    "\n")



# ==============================================================================
# Define background genes
# Important:
# use all tested genes, NOT network genes
# ==============================================================================


expr <- fread(
  expr_file,
  header=TRUE,
  sep="\t",
  check.names=FALSE
)


background <- expr[[1]]

background <- gsub(
  "\\..*",
  "",
  background
)

background <- unique(background)



cat(
  "Background genes:",
  length(background),
  "\n"
)



# ==============================================================================
# GO enrichment function
# ==============================================================================


run_GO <- function(
    gene_list,
    background,
    prefix
){

  genes <- unique(
    gsub("\\..*", "", gene_list)
  )

  bg <- unique(
    gsub("\\..*", "", background)
  )


  # ENSEMBL -> ENTREZ

  gene_entrez <- bitr(
    genes,
    fromType="ENSEMBL",
    toType="ENTREZID",
    OrgDb=org.Hs.eg.db
  )


  bg_entrez <- bitr(
    bg,
    fromType="ENSEMBL",
    toType="ENTREZID",
    OrgDb=org.Hs.eg.db
  )


  cat(
    prefix,
    "mapped genes:",
    nrow(gene_entrez),
    "\n"
  )


  ego <- enrichGO(
    gene = gene_entrez$ENTREZID,
    universe = bg_entrez$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )


  if(!is.null(ego) && nrow(ego@result)>0){

    result <- as.data.frame(ego)


    write.table(
      result,
      paste0(
        "C4A_GO_Reproducible_Postnatal_",
        prefix,
        ".txt"
      ),
      sep="\t",
      quote=FALSE,
      row.names=FALSE
    )


  }


  return(ego)

}



# ==============================================================================
# Run GO enrichment
# ==============================================================================


cat("\nRunning Positive network GO...\n")


pos_GO <- run_GO(
  pos_net$Gene,
  background,
  "Positive"
)



cat("\nRunning Negative network GO...\n")


neg_GO <- run_GO(
  neg_net$Gene,
  background,
  "Negative"
)





# ==============================================================================
# Plot function
# ==============================================================================


plot_GO <- function(
    ego,
    title,
    filename
){

  if(is.null(ego)) return(NULL)


  df <- as.data.frame(ego)


  df <- df %>%
    filter(p.adjust <0.05) %>%
    arrange(p.adjust) %>%
    head(10)


  if(nrow(df)==0)
    return(NULL)


  df$Description <- factor(
    df$Description,
    levels=rev(df$Description)
  )


  p <- ggplot(
    df,
    aes(
      x=-log10(p.adjust),
      y=Description
    )
  )+
    geom_col(
      fill="#83bcd7",
      alpha=0.8
    )+
    geom_point(
      aes(size=Count),
      color="#2166ac"
    )+
    theme_classic(
      base_size=14
    )+
    labs(
      title=title,
      x="-log10(FDR)",
      y=NULL,
      size="Gene number"
    )+
    theme(
      plot.title=
        element_text(
          hjust=0.5,
          face="bold"
        )
    )


  ggsave(
    filename,
    p,
    width=7,
    height=5,
    dpi=400
  )


  return(p)

}



# ==============================================================================
# Generate figures
# ==============================================================================


p1 <- plot_GO(
  pos_GO,
  "Reproducible C4A Positive Network (Postnatal N=170)",
  "GO_Reproducible_Postnatal_Positive.png"
)


p2 <- plot_GO(
  neg_GO,
  "Reproducible C4A Negative Network (Postnatal N=170)",
  "GO_Reproducible_Postnatal_Negative.png"
)



cat("\nGO analysis completed\n")