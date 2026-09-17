# 1.  prenatal: C4A seed co-expression network
# 2. Goal: Standardize Postnatal power to match Prenatal (N=170) | downsampling & resampling 
# 3. Goal: Build tiered seed networks based on effect size for prenatal (|R|)
# 4. Goal: Build tiered seed networks based on effect size for postnatal (After downsampling) (|R|)

# ==============================================================================
# 1.  prenatal: C4A seed co-expression network
# ==============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(ggplot2)
library(corrplot)
library(data.table)
library(biomaRt)
expr_file <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/bulk_data/stage1.logTMM.ComBat.txt"
cov_file <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/bulk_data/prenatalcovariatesToUse.PCAforQTLfromINT.txt"
#expr_file <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/bulk_data/postnatal.logTMM.ComBat.txt"
#cov_file <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/bulk_data/postnatalcovariatesToUse.PCAforQTLfromINT.txt"

datExpr <- fread(expr_file, header = TRUE, sep = "\t", check.names = FALSE)
gene_ids <- datExpr[[1]]
datExpr <- datExpr[, -1, with = FALSE]
datExpr <- as.matrix(datExpr)
rownames(datExpr) <- gene_ids
## Or: datExpr<-  read.table(expr_file, header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
datCov <- fread(cov_file, header = TRUE, sep = " ", check.names = FALSE)
# datCov <- fread(cov_file, header = TRUE, sep = "\t", check.names = FALSE)

cov_rownames <- datCov$cov
datCov <- datCov[, -1, with = FALSE]
datCov <- as.matrix(datCov)
rownames(datCov) <- cov_rownames
datCov <- t(datCov)
datCov <- as.data.frame(datCov)

# sample match
common_samples <- intersect(colnames(datExpr), rownames(datCov))
datExpr <- datExpr[, common_samples]
datCov <- datCov[common_samples, ]

# --- 2. Regress out Covariates (Excluding Expression PCs) ---
# We keep sex, RIN, GW, and expression_pcs to remove technical/demographic noise
datCov<- as.data.frame(datCov)
covs_to_regress <- as.data.frame(datCov %>% dplyr::select(sex, RIN, GW, PC1,PC2,PC3))
# For postntal: covs_to_regress <- as.data.frame(datCov %>% dplyr::select(sex, RIN, Age, PC1,PC2,PC3))
# Fast matrix regression
expr<- as.data.frame(datExpr)
fit <- lm(t(expr) ~ ., data = covs_to_regress)
# expr_res <- t(residuals(fit)) + rowMeans(expr)
expr_res <- t(residuals(fit))

# --- 3. Calculate C4A Co-expression ---
# Identify C4A row (update "ENSG00000244731" if your matrix uses symbols)
c4a_id <- grep("ENSG00000244731.9|C4A", rownames(expr_res), value = TRUE)[1]
c4a_expr <- expr_res[c4a_id, ]

# Calculate Pearson r and p-values for all genes simultaneously
r_vals <- cor(t(expr_res), c4a_expr)[, 1]
n <- ncol(expr_res)
t_stats <- r_vals * sqrt((n - 2) / (1 - r_vals^2))
p_vals <- 2 * pt(-abs(t_stats), df = n - 2)
fdr <- p.adjust(p_vals, method = "BH")

# Create Network DataFrame
network_df <- data.frame(Gene = rownames(expr_res), R = r_vals, P = p_vals, FDR = fdr) %>%
  filter(Gene != c4a_id)

# Define Positive and Negative Networks (FDR < 0.05)
pos_net <- network_df %>% filter(FDR < 0.05 & R > 0.3)
neg_net <- network_df %>% filter(FDR < 0.05 & R < -0.3)


write.csv(pos_net, "C4A_Positive_Network_pre03.csv", row.names = FALSE)
write.csv(neg_net, "C4A_Negative_Network_pre03.csv", row.names = FALSE)
write.table(expr_res, "prenatal_expression.residual.txt", col.names=NA,quote=F)
## volcano plot

library(tidyverse)
library(ggplot2)
library(patchwork)

volcano_data <- network_df %>%
  mutate(
    log10_fdr = -log10(FDR),
    group = case_when(
      FDR < 0.05 & R > 0 ~ "Positive",
      FDR < 0.05 & R < 0 ~ "Negative",
      TRUE ~ "NS"
    )
  )


volcano_data$group <- factor(volcano_data$group, levels = c("Positive", "Negative", "NS"))

p_hist <- ggplot(volcano_data, aes(x = R)) +
  geom_histogram(bins = 50, fill = "gray70", color = "black", alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
  theme_classic(base_size = 14) +
  labs(
    title = "C4A Correlation Distribution (Prenatal)",
    x = "Pearson Correlation (R)",
    y = "Frequency"
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


p_volcano <- ggplot(volcano_data, aes(x = R, y = log10_fdr, color = group)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_color_manual(values = c("Positive" = "#d7191c", "Negative" = "#2c7bb6", "NS" = "gray80")) +
  theme_classic(base_size = 14) +
  labs(
    title = "C4A Seed Network Volcano Plot (Prenatal)",
    x = "Pearson Correlation (R)",
    y = expression(-Log[10](FDR)),
    color = "Significance"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  )

final_plot <- p_hist | p_volcano
print(final_plot)
ggsave("C4A_correlation_distribution_and_volcano_pre.png", final_plot, width = 11, height = 5, dpi = 400)

pdf("./C4A_correlation_distribution_and_volcano.pdf", width = 9, height = 5)
print(final_plot)
dev.off()


# Using clusterProfiler, translate Ensembl IDs into Entrez IDs and run the enrichments.

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(stringr)


run_enrichment <- function(gene_list, background, prefix) {

  genes_clean <- gsub("\\..*", "", gene_list)
  bg_clean <- gsub("\\..*", "", background)

  gene_entrez <- tryCatch(
    bitr(genes_clean, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID,
    error = function(e) return(NULL)
  )
  bg_entrez <- tryCatch(
    bitr(bg_clean, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID,
    error = function(e) return(NULL)
  )
  
  if(is.null(gene_entrez) || length(gene_entrez) < 5) return(NULL)
  ego <- enrichGO(gene = gene_entrez, universe = bg_entrez, OrgDb = org.Hs.eg.db, 
                  ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, 
                  qvalueCutoff = 0.2, readable = TRUE)
  if (!is.null(ego) && nrow(ego@result) > 0) {
    write.table(ego@result, paste0("C4A_GO_Result_prenatal_", prefix, ".txt"), 
                quote = FALSE, sep = "\t", row.names = FALSE)
  }
  return(ego)
}

plot_dual_axis_label_on_bar <- function(enrich_obj, color_theme, title_text, top_n = 10, keywords = NULL) {
  
  if (is.null(enrich_obj) || nrow(enrich_obj@result) == 0) return(NULL)
  
  df <- enrich_obj@result %>% filter(p.adjust < 0.05)
  
  if (!is.null(keywords) && length(keywords) > 0) {
    pattern <- paste(keywords, collapse = "|")
    df <- df %>% filter(grepl(pattern, Description, ignore.case = TRUE))
  }
  
  if (nrow(df) == 0) {
    message("No significant pathway was found for: ", title_text)
    return(NULL)
  }

  df <- df %>%
    arrange(desc(Count)) %>% 
    slice_head(n = top_n) %>%
    mutate(
      GeneRatioNum = as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio)) * 100,
      MinusLog10P = -log10(p.adjust),
      Description = str_wrap(Description, width = 50) 
    ) %>%
    arrange(Count) 
  
  df$Description <- factor(df$Description, levels = df$Description)
  max_p <- max(df$MinusLog10P)
  max_gr <- max(df$GeneRatioNum)
  coeff <- if(max_gr == 0) 1 else max_p / max_gr
  
  p <- ggplot(df, aes(y = Description)) +
    geom_col(aes(x = GeneRatioNum * coeff), fill = color_theme, alpha = 0.3, width = 0.8) +
    geom_path(aes(x = MinusLog10P, group = 1), color = "grey50", linewidth = 0.8) +
    geom_point(aes(x = MinusLog10P, size = Count), color = color_theme) +
    geom_text(aes(x = 0.05, label = Description), 
              hjust = 0, size = 4.5, color = "black", lineheight = 0.9) + 
    scale_x_continuous(
      name = expression(-Log[10](p.adjust)),
      sec.axis = sec_axis(~ . / coeff, name = "GeneRatio (%)"),
      expand = expansion(mult = c(0, 0.1)) 
    ) +
    scale_size_continuous(range = c(4, 10)) + 
    labs(title = title_text, y = NULL) +
    theme_classic(base_size = 15) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(), 
      axis.text.x = element_text(color = "black", size = 12),
      axis.title.x.top = element_text(vjust = 2),
      legend.position = "bottom",
      plot.margin = margin(t = 20, r = 20, b = 20, l = 10)
    )
  return(p)
}


bg_genes <- network_df$Gene
pos_go <- run_enrichment(pos_net$Gene, bg_genes, "Positive")
cat("Running GO enrichment for Positive Network...\n")
cat("Running GO enrichment for Negative Network...\n")
neg_go <- run_enrichment(neg_net$Gene, bg_genes, "Negative")


if(!is.null(pos_go)) {
  p1 <- plot_dual_axis_label_on_bar(
    enrich_obj = pos_go, 
    color_theme = "#D6604D",     
    title_text = "C4A Positive Seed Network (Prenatal)", 
    top_n = 8, 
    keywords = NULL              
  )
  
  if(!is.null(p1)) {
    print(p1)
    ggsave("C4A_GO_prenatal_Positive.png", p1, width = 6.5, height = 6, dpi = 400)
    #ggsave("C4A_GO_prenatal_Positive.pdf", p1, width = 7, height = 6)
  }
}

if(!is.null(neg_go)) {
  p2 <- plot_dual_axis_label_on_bar(
    enrich_obj = neg_go, 
    color_theme = "#D6604D",     
    title_text = "C4A Negative Seed Network (Prenatal)", 
    top_n = 8, 
    keywords = NULL             
  )
 #  #4393C3 postantal  ; #D6604D #prenatal
 #  #83bcd7 postantal  ; #e89a92 #prenatal
  if(!is.null(p2)) {
    print(p2)
    ggsave("C4A_GO_Prenatal_Negative.png", p2, width = 6.5, height = 6, dpi = 400)
    #ggsave("C4A_GO_prenatal_Negative.pdf", p2, width = 7, height = 6)
  }
}


# ==============================================================================
# 2a.  postnatal: C4A seed co-expression network
# ==============================================================================

expr_file <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/bulk_data/postnatal.logTMM.ComBat.txt"
cov_file <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/bulk_data/postnatalcovariatesToUse.PCAforQTLfromINT.txt"
gene<- read.table("/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/01_data/RNA-seq_bak/postnatal/genes_to_keep.txt")
datExpr <- fread(expr_file, header = TRUE, sep = "\t", check.names = FALSE)
gene_ids <- datExpr[[1]]
datExpr <- datExpr[, -1, with = FALSE]
datExpr <- as.matrix(datExpr)
rownames(datExpr) <- gene_ids
#prenatal_genes<- rownames(prenatal_exp)
datExpr<- datExpr[gene$V1,]
## Or: datExpr<-  read.table(expr_file, header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
datCov <- fread(cov_file, header = TRUE, sep = "\t", check.names = FALSE)
# datCov <- fread(cov_file, header = TRUE, sep = "\t", check.names = FALSE)

cov_rownames <- datCov$cov
datCov <- datCov[, -1, with = FALSE]
datCov <- as.matrix(datCov)
rownames(datCov) <- cov_rownames
datCov <- t(datCov)
datCov <- as.data.frame(datCov)

# sample match
common_samples <- intersect(colnames(datExpr), rownames(datCov))
datExpr <- datExpr[, common_samples]
datCov <- datCov[common_samples, ]

# --- 2. Regress out Covariates (Excluding Expression PCs) ---
# We keep sex, RIN, GW, and expression_pcs to remove technical/demographic noise
datCov<- as.data.frame(datCov)
covs_to_regress <- as.data.frame(datCov %>% dplyr::select(sex, RIN, Age, PC1,PC2,PC3))
# For postntal: covs_to_regress <- as.data.frame(datCov %>% dplyr::select(sex, RIN, Age, PC1,PC2,PC3))
# Fast matrix regression
expr<- as.data.frame(datExpr)
fit <- lm(t(expr) ~ ., data = covs_to_regress)
expr_res <- t(residuals(fit)) + rowMeans(expr)

# --- 3. Calculate C4A Co-expression ---
# Identify C4A row (update "ENSG00000244731" if your matrix uses symbols)
c4a_id <- grep("ENSG00000244731.9|C4A", rownames(expr_res), value = TRUE)[1]
c4a_expr <- expr_res[c4a_id, ]

# Calculate Pearson r and p-values for all genes simultaneously
r_vals <- cor(t(expr_res), c4a_expr)[, 1]
n <- ncol(expr_res)
t_stats <- r_vals * sqrt((n - 2) / (1 - r_vals^2))
p_vals <- 2 * pt(-abs(t_stats), df = n - 2)
fdr <- p.adjust(p_vals, method = "BH")

# Create Network DataFrame
network_df <- data.frame(Gene = rownames(expr_res), R = r_vals, P = p_vals, FDR = fdr) %>%
  filter(Gene != c4a_id)

# Define Positive and Negative Networks (FDR < 0.05)
pos_net <- network_df %>% filter(FDR < 0.05 & R > 0)
neg_net <- network_df %>% filter(FDR < 0.05 & R < 0)

write.csv(pos_net, "C4A_Positive_Network_post03.csv", row.names = FALSE)
write.csv(neg_net, "C4A_Negative_Network_post03.csv", row.names = FALSE)
write.table(expr_res, "postnatal_expression.residual.txt", col.names=NA,quote=F)
## volcano plot

library(tidyverse)
library(ggplot2)
library(patchwork)

volcano_data <- network_df %>%
  mutate(
    log10_fdr = -log10(FDR),
    group = case_when(
      FDR < 0.05 & R > 0 ~ "Positive",
      FDR < 0.05 & R < 0 ~ "Negative",
      TRUE ~ "NS"
    )
  )


volcano_data$group <- factor(volcano_data$group, levels = c("Positive", "Negative", "NS"))

p_hist <- ggplot(volcano_data, aes(x = R)) +
  geom_histogram(bins = 50, fill = "gray70", color = "black", alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
  theme_classic(base_size = 14) +
  labs(
    title = "C4A Correlation Distribution (postnatal)",
    x = "Pearson Correlation (R)",
    y = "Frequency"
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


p_volcano <- ggplot(volcano_data, aes(x = R, y = log10_fdr, color = group)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_color_manual(values = c("Positive" = "#d7191c", "Negative" = "#2c7bb6", "NS" = "gray80")) +
  theme_classic(base_size = 14) +
  labs(
    title = "C4A Seed Network Volcano Plot (postnatal)",
    x = "Pearson Correlation (R)",
    y = expression(-Log[10](FDR)),
    color = "Significance"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  )

final_plot <- p_hist | p_volcano
print(final_plot)
ggsave("C4A_correlation_distribution_and_volcano_post.png", final_plot, width = 11, height = 5, dpi = 400)

pdf("./C4A_correlation_distribution_and_post_volcano.pdf", width = 9, height = 5)
print(final_plot)
dev.off()


# Using clusterProfiler, translate Ensembl IDs into Entrez IDs and run the enrichments.

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(stringr)


run_enrichment <- function(gene_list, background, prefix) {

  genes_clean <- gsub("\\..*", "", gene_list)
  bg_clean <- gsub("\\..*", "", background)

  gene_entrez <- tryCatch(
    bitr(genes_clean, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID,
    error = function(e) return(NULL)
  )
  bg_entrez <- tryCatch(
    bitr(bg_clean, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID,
    error = function(e) return(NULL)
  )
  
  if(is.null(gene_entrez) || length(gene_entrez) < 5) return(NULL)
  ego <- enrichGO(gene = gene_entrez, universe = bg_entrez, OrgDb = org.Hs.eg.db, 
                  ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, 
                  qvalueCutoff = 0.2, readable = TRUE)
  if (!is.null(ego) && nrow(ego@result) > 0) {
    write.table(ego@result, paste0("C4A_GO_Result_postnatal", prefix, ".txt"), 
                quote = FALSE, sep = "\t", row.names = FALSE)
  }
  return(ego)
}

plot_dual_axis_label_on_bar <- function(enrich_obj, color_theme, title_text, top_n = 10, keywords = NULL) {
  
  if (is.null(enrich_obj) || nrow(enrich_obj@result) == 0) return(NULL)
  
  df <- enrich_obj@result %>% filter(p.adjust < 0.05)
  
  if (!is.null(keywords) && length(keywords) > 0) {
    pattern <- paste(keywords, collapse = "|")
    df <- df %>% filter(grepl(pattern, Description, ignore.case = TRUE))
  }
  
  if (nrow(df) == 0) {
    message("No significant pathway was found for: ", title_text)
    return(NULL)
  }

  df <- df %>%
    arrange(desc(Count)) %>% 
    slice_head(n = top_n) %>%
    mutate(
      GeneRatioNum = as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio)) * 100,
      MinusLog10P = -log10(p.adjust),
      Description = str_wrap(Description, width = 50) 
    ) %>%
    arrange(Count) 
  
  df$Description <- factor(df$Description, levels = df$Description)
  max_p <- max(df$MinusLog10P)
  max_gr <- max(df$GeneRatioNum)
  coeff <- if(max_gr == 0) 1 else max_p / max_gr
  
  p <- ggplot(df, aes(y = Description)) +
    geom_col(aes(x = GeneRatioNum * coeff), fill = color_theme, alpha = 0.3, width = 0.8) +
    geom_path(aes(x = MinusLog10P, group = 1), color = "grey50", linewidth = 0.8) +
    geom_point(aes(x = MinusLog10P, size = Count), color = color_theme) +
    geom_text(aes(x = 0.05, label = Description), 
              hjust = 0, size = 4.5, color = "black", lineheight = 0.9) + 
    scale_x_continuous(
      name = expression(-Log[10](p.adjust)),
      sec.axis = sec_axis(~ . / coeff, name = "GeneRatio (%)"),
      expand = expansion(mult = c(0, 0.1)) 
    ) +
    scale_size_continuous(range = c(4, 10)) + 
    labs(title = title_text, y = NULL) +
    theme_classic(base_size = 15) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(), 
      axis.text.x = element_text(color = "black", size = 12),
      axis.title.x.top = element_text(vjust = 2),
      legend.position = "bottom",
      plot.margin = margin(t = 20, r = 20, b = 20, l = 10)
    )
  return(p)
}


bg_genes <- network_df$Gene
pos_go <- run_enrichment(pos_net$Gene, bg_genes, "Positive")
cat("Running GO enrichment for Positive Network...\n")
cat("Running GO enrichment for Negative Network...\n")
neg_go <- run_enrichment(neg_net$Gene, bg_genes, "Negative")


if(!is.null(pos_go)) {
  p1 <- plot_dual_axis_label_on_bar(
    enrich_obj = pos_go, 
    color_theme = "#83bcd7",     
    title_text = "C4A Positive Seed Network (postnatal)", 
    top_n = 8, 
    keywords = NULL              
  )
  
  if(!is.null(p1)) {
    print(p1)
    ggsave("C4A_GO_postnatal_Positive.png", p1, width = 6, height = 5.5, dpi = 400)
    #ggsave("C4A_GO_postnatal_Positive.pdf", p1, width = 7, height = 6)
  }
}

if(!is.null(neg_go)) {
  p2 <- plot_dual_axis_label_on_bar(
    enrich_obj = neg_go, 
    color_theme = "#83bcd7",     
    title_text = "C4A Negative Seed Network (postnatal)", 
    top_n = 8, 
    keywords = NULL             
  )
 #  #4393C3 postantal  ; #D6604D #postnatal
 #  #83bcd7 postantal  ; #e89a92 #postnatal
  if(!is.null(p2)) {
    print(p2)
    ggsave("C4A_GO_postnatal_Negative.png", p2, width = 6, height =5.5, dpi = 400)
    #ggsave("C4A_GO_postnatal_Negative.pdf", p2, width = 7, height = 6)
  }
}


# ==============================================================================
# 2. postnatal: C4A co-expression network: downsampling & resampling 100 
# Goal: Standardize Postnatal power to match Prenatal (N=170)
# ==============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

library(data.table)
library(dplyr)
library(tidyverse)

n_iterations <- 100
target_n     <- 170 # Target size to match prenatal dataset
expr_file    <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/bulk_data/postnatal.logTMM.ComBat.txt"
cov_file     <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/bulk_data/postnatalcovariatesToUse.PCAforQTLfromINT.txt"
#expr_file    <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/bulk_data/newpostnatal.logTMM.ComBat.txt"
#cov_file     <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/bulk_data/newpostnatalcovariatesToUse.PCAforQTLfromINT.txt"
output_name  <- "C4A_newResampled_Median_Network_Postnatal_N170.csv"

message("Loading high-throughput data...")
# Load expression matrix
datExpr <- fread(expr_file, header = TRUE, sep = "\t", check.names = FALSE)
gene_ids <- datExpr[[1]]
datExpr  <- as.matrix(datExpr[, -1, with = FALSE])
rownames(datExpr) <- gene_ids

# Load and transpose covariates
datCov <- fread(cov_file, header = TRUE, check.names = FALSE)
cov_rownames <- datCov[[1]]
datCov <- as.matrix(datCov[, -1, with = FALSE])
rownames(datCov) <- cov_rownames
datCov <- as.data.frame(t(datCov))

# Synchronize sample IDs
common_samples <- intersect(colnames(datExpr), rownames(datCov))
datExpr <- datExpr[, common_samples]
datCov  <- datCov[common_samples, ]

# Identify C4A seed (Target Gene)
c4a_id <- grep("ENSG00000244731.9|C4A", rownames(datExpr), value = TRUE)[1]

# --- 3. Resampling Loop ---
message("Starting ", n_iterations, " iterations of bootstrapping...")

# Storage matrix for R values: [Genes x Iterations]
all_r_matrix <- matrix(NA, nrow = nrow(datExpr), ncol = n_iterations)
rownames(all_r_matrix) <- rownames(datExpr)

for (i in 1:n_iterations) {
  set.seed(i) # Ensure reproducibility
  
  # 3.1 Randomly downsample to N=170
  selected_samples <- sample(common_samples, target_n)
  sub_expr <- datExpr[, selected_samples]
  sub_cov  <- datCov[selected_samples, ]
  
  # 3.2 Prepare covariates for regression
  # Note: Postnatal uses 'Age' instead of 'GW' (Gestational Week)
  # Modify selection based on your postnatal covariate file headers
  covs_to_regress <- sub_cov %>% 
    dplyr::select(sex, RIN, Age, PC1, PC2, PC3) 
  
  # 3.3 Regress out technical noise in the subset
  # Why: Covariate structures vary across random subsets; local correction is safer.
  fit <- lm(t(sub_expr) ~ ., data = covs_to_regress)
  expr_res <- t(residuals(fit)) + rowMeans(sub_expr)
  
  # 3.4 Compute Pearson Correlation for C4A
  sub_c4a <- expr_res[c4a_id, ]
  r_vals  <- cor(t(expr_res), sub_c4a)[, 1]
  
  all_r_matrix[, i] <- r_vals
  
  if(i %% 10 == 0) message("Iteration ", i, " completed.")
}

# --- 4. Final Aggregation ---
message("Aggregating median results...")

# Use Median R for robustness against sampling outliers
#median_r <- apply(all_r_matrix, 1, median, na.rm = TRUE)

# Back-calculate P-values using N=170 (df = N-2)
# Formula: t = r * sqrt((n-2)/(1-r^2))
#n <- target_n
#t_stats <- median_r * sqrt((n - 2) / (1 - median_r^2))
#p_vals  <- 2 * pt(-abs(t_stats), df = n - 2)
#fdr     <- p.adjust(p_vals, method = "BH")


mean_z  <- apply(all_z_matrix, 1, mean, na.rm = TRUE)
final_r <- fisherz2r(mean_z) # 转换回重采样后的 R 值
se_z     <- 1 / sqrt(target_n - 3)
z_scores <- mean_z / se_z
p_vals   <- 2 * (1 - pnorm(abs(z_scores)))
fdr      <- p.adjust(p_vals, method = "BH")

# --- 5. Export Results ---
final_network <- data.frame(
  Gene = gsub("\\..*", "", rownames(all_r_matrix)),                    
  Median_R = median_r,
  P_value  = p_vals,
  FDR      = fdr
) %>%
  filter(Gene != c4a_id) %>% # Remove seed self-correlation
  arrange(desc(Median_R))

write.csv(final_network, output_name, row.names = FALSE)
message("Pipeline finished. Results saved to: ", output_name)

# --- 6. Quick Summary Stats ---
sig_pos <- sum(final_network$FDR < 0.05 & final_network$Median_R > 0)
sig_neg <- sum(final_network$FDR < 0.05 & final_network$Median_R < 0)
network_df<- final_network
pos_net <- network_df %>% filter(FDR < 0.05 & Median_R > 0)
neg_net <- network_df %>% filter(FDR < 0.05 & Median_R < 0)
cat("\nSummary (FDR < 0.05):\nPositive Connections:", sig_pos, "\nNegative Connections:", sig_neg, "\n")


# Using clusterProfiler, translate Ensembl IDs into Entrez IDs and run the enrichments.

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(stringr)


run_enrichment <- function(gene_list, background, prefix) {

  genes_clean <- gsub("\\..*", "", gene_list)
  bg_clean <- gsub("\\..*", "", background)

  gene_entrez <- tryCatch(
    bitr(genes_clean, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID,
    error = function(e) return(NULL)
  )
  bg_entrez <- tryCatch(
    bitr(bg_clean, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID,
    error = function(e) return(NULL)
  )
  
  if(is.null(gene_entrez) || length(gene_entrez) < 5) return(NULL)
  ego <- enrichGO(gene = gene_entrez, universe = bg_entrez, OrgDb = org.Hs.eg.db, 
                  ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, 
                  qvalueCutoff = 0.2, readable = TRUE)
  if (!is.null(ego) && nrow(ego@result) > 0) {
    write.table(ego@result, paste0("C4A_GO_Result_", "postnatal",prefix, ".txt"), 
                quote = FALSE, sep = "\t", row.names = FALSE)
  }
  
  return(ego)
}

plot_dual_axis_label_on_bar <- function(enrich_obj, color_theme, title_text, top_n = 10, keywords = NULL) {
  
  if (is.null(enrich_obj) || nrow(enrich_obj@result) == 0) return(NULL)
  
  df <- enrich_obj@result %>% filter(p.adjust < 0.05)
  
  if (!is.null(keywords) && length(keywords) > 0) {
    pattern <- paste(keywords, collapse = "|")
    df <- df %>% filter(grepl(pattern, Description, ignore.case = TRUE))
  }
  
  if (nrow(df) == 0) {
    message("No significant pathway was found for: ", title_text)
    return(NULL)
  }

  df <- df %>%
    arrange(desc(Count)) %>% 
    slice_head(n = top_n) %>%
    mutate(
      GeneRatioNum = as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio)) * 100,
      MinusLog10P = -log10(p.adjust),
      Description = str_wrap(Description, width = 50) 
    ) %>%
    arrange(Count) 
  
  df$Description <- factor(df$Description, levels = df$Description)
  max_p <- max(df$MinusLog10P)
  max_gr <- max(df$GeneRatioNum)
  coeff <- if(max_gr == 0) 1 else max_p / max_gr
  
  p <- ggplot(df, aes(y = Description)) +
    geom_col(aes(x = GeneRatioNum * coeff), fill = color_theme, alpha = 0.3, width = 0.8) +
    geom_path(aes(x = MinusLog10P, group = 1), color = "grey50", linewidth = 0.8) +
    geom_point(aes(x = MinusLog10P, size = Count), color = color_theme) +
    geom_text(aes(x = 0.05, label = Description), 
              hjust = 0, size = 4.5, color = "black", lineheight = 0.9) + 
    scale_x_continuous(
      name = expression(-Log[10](p.adjust)),
      sec.axis = sec_axis(~ . / coeff, name = "GeneRatio (%)"),
      expand = expansion(mult = c(0, 0.1)) 
    ) +
    scale_size_continuous(range = c(4, 10)) + 
    labs(title = title_text, y = NULL) +
    theme_classic(base_size = 15) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(), 
      axis.text.x = element_text(color = "black", size = 12),
      axis.title.x.top = element_text(vjust = 2),
      legend.position = "bottom",
      plot.margin = margin(t = 20, r = 20, b = 20, l = 10)
    )
  return(p)
}


bg_genes <- network_df$Gene
pos_go <- run_enrichment(pos_net$Gene, bg_genes, "Positive")
cat("Running GO enrichment for Positive Network...\n")
cat("Running GO enrichment for Negative Network...\n")
neg_go <- run_enrichment(neg_net$Gene, bg_genes, "Negative")


if(!is.null(pos_go)) {
  p1 <- plot_dual_axis_label_on_bar(
    enrich_obj = pos_go, 
    color_theme = "#83bcd7",     
    title_text = "C4A Positive Seed Network (Postnatal)", 
    top_n = 8, 
    keywords = NULL              
  )
  
  if(!is.null(p1)) {
    print(p1)
    ggsave("C4A_GO_postnatal_Positive.png", p1, width = 6, height = 5, dpi = 400)
    #ggsave("C4A_GO_prenatal_Positive.pdf", p1, width = 7, height = 6)
  }
}


if(!is.null(neg_go)) {
  p2 <- plot_dual_axis_label_on_bar(
    enrich_obj = neg_go, 
    color_theme = "#83bcd7",     
    title_text = "C4A Negative Seed Network (Postnatal)", 
    top_n = 8, 
    keywords = NULL             
  )
 #  #4393C3 postantal  ; #D6604D #prenatal
 #  #83bcd7 postantal  ; #e89a92 #prenatal
  if(!is.null(p2)) {
    print(p2)
    ggsave("C4A_GO_Postnatal_Negative.png", p2, width = 6, height = 6, dpi = 400)
    #ggsave("C4A_GO_prenatal_Negative.pdf", p2, width = 7, height = 6)
  }
}

# ==============================================================================
# 3. Goal: Build tiered seed networks based on effect size for prenatal (|R|)
# ==============================================================================


library(data.table)
library(dplyr)
library(tidyverse)

expr_file   <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/bulk_data/stage1.logTMM.ComBat.txt"
cov_file    <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/bulk_data/prenatalcovariatesToUse.PCAforQTLfromINT.txt"
output_pref <- "C4A_Network_tieredFullSet_Prenatal" 

# --- 2. Data Loading & Pre-processing ---
message("Loading and synchronizing data...")

datExpr <- fread(expr_file, header = TRUE, sep = "\t", check.names = FALSE)
gene_ids <- datExpr[[1]]
datExpr  <- as.matrix(datExpr[, -1, with = FALSE])
rownames(datExpr) <- gene_ids

datCov <- fread(cov_file, header = TRUE, check.names = FALSE)
cov_rownames <- datCov[[1]]
datCov <- as.matrix(datCov[, -1, with = FALSE])
rownames(datCov) <- cov_rownames
datCov <- as.data.frame(t(datCov))

# Match samples across expression and covariates
common_samples <- intersect(colnames(datExpr), rownames(datCov))
datExpr <- datExpr[, common_samples]
datCov  <- datCov[common_samples, ]

# Identify C4A Seed (ENSG00000244731)
c4a_id <- grep("ENSG00000244731.9|C4A", rownames(datExpr), value = TRUE)[1]

# --- 3. Covariate Regression (Noise Removal) ---
message("Regressing out technical covariates...")

# Ensure covariate selection matches your specific dataset (e.g., Age vs GW)
covs_to_regress <- datCov %>% 
  dplyr::select(sex, RIN, GW, PC1, PC2, PC3) 

# Perform global linear regression to get residuals
fit      <- lm(t(datExpr) ~ ., data = covs_to_regress)
# expr_res <- t(residuals(fit)) + rowMeans(datExpr)
expr_res <- t(residuals(fit)) 
# --- 4. Pearson Correlation Analysis ---
message("Calculating Pearson correlations for C4A...")

c4a_profile <- expr_res[c4a_id, ]
r_values    <- cor(t(expr_res), c4a_profile)[, 1]

# Statistical inference (N = total samples)
n       <- ncol(expr_res)
t_stats <- r_values * sqrt((n - 2) / (1 - r_values^2))
p_values <- 2 * pt(-abs(t_stats), df = n - 2)
fdr_vals <- p.adjust(p_values, method = "BH")

# Construct base dataframe
network_full <- data.frame(
  Gene = rownames(expr_res),
  R    = r_values,
  P    = p_values,
  FDR  = fdr_vals
) %>%
  filter(Gene != c4a_id) # Exclude self-correlation
write.csv(network_full, "C4A_Network_Prenatal_N170.csv", row.names = FALSE)
# --- 5. Dynamic Thresholding (|R| > 0, 0.3, 0.5) ---
message("Applying dynamic thresholds...")

# Function to extract and save network tiers
extract_tier <- function(df, r_thresh, prefix) {
  tier_pos <- df %>% filter(FDR < 0.05 & R > r_thresh)
  tier_neg <- df %>% filter(FDR < 0.05 & R < -r_thresh)
  
  write.csv(tier_pos, paste0(prefix, "_R", r_thresh, "_Positive.csv"), row.names = FALSE)
  write.csv(tier_neg, paste0(prefix, "_R", r_thresh, "_Negative.csv"), row.names = FALSE)
  
  return(data.frame(Threshold = r_thresh, Pos_Count = nrow(tier_pos), Neg_Count = nrow(tier_neg)))
}

# Run extraction for 0.0, 0.3, and 0.5
stats_00 <- extract_tier(network_full, 0.0, output_pref)
stats_03 <- extract_tier(network_full, 0.3, output_pref)
stats_05 <- extract_tier(network_full, 0.5, output_pref)

# --- 6. Summary Output ---
summary_table <- rbind(stats_00, stats_03, stats_05)
print(summary_table)

#  Threshold Pos_Count Neg_Count
#1       0.0      1469      1513
#2       0.3       404       313
#3       0.5         5         1

message("Analysis complete. CSV files generated for all tiers.")

# --- 1. Define thresholds and directions ---
thresholds <- c(0.0, 0.3, 0.5)
directions <- c("Positive", "Negative")
bg_genes   <- network_full$Gene # Background is all genes tested

# --- 2. Automated Loop for Enrichment ---
for (tr in thresholds) {
  for (dir in directions) {
    
    cat(sprintf("\nProcessing Tier: R > %s | Direction: %s\n", tr, dir))
    
    # Filter genes based on current criteria
    if (dir == "Positive") {
      current_list <- network_full %>% filter(FDR < 0.05 & R > tr) %>% pull(Gene)
    } else {
      current_list <- network_full %>% filter(FDR < 0.05 & R < -tr) %>% pull(Gene)
    }
    
    # Check if there are enough genes to run GO
    if (length(current_list) < 10) {
      message("Skipping: Too few genes in this category.")
      next
    }
    
    # Run Enrichment (using your function)
    label <- paste0("R", tr, "_", dir)
    ego_res <- run_enrichment(current_list, bg_genes, label)
    
    # Plot results if significant
    if (!is.null(ego_res)) {
      # Choose color theme based on direction
      color_theme <- ifelse(dir == "Positive", "#D6604D", "#D6604D")
      
      p <- plot_dual_axis_label_on_bar(
        enrich_obj  = ego_res, 
        color_theme = color_theme, 
        title_text  = paste("C4A", dir, "Seed Network (Prenatal)"),                                                
        top_n       = 8
      )
      
      if (!is.null(p)) {
        file_name <- sprintf("C4A_GO_Prenatal_%s_R%s.png", dir, tr)
        ggsave(file_name, p, width = 7, height = 6, dpi = 400)
        message("Saved: ", file_name)
      }
    }
  }
}

# ==============================================================================
# 4. Goal: Build tiered seed networks based on effect size for postnatal (After downsampling) (|R|)
# ==============================================================================

# Get final_network of postantal
output_post <- "C4A_Network_tieredFullSet_Postnatal" 

# Function to extract and save network tiers
extract_tier <- function(df, r_thresh, prefix) {
  tier_pos <- df %>% filter(FDR < 0.05 & Median_R > r_thresh)
  tier_neg <- df %>% filter(FDR < 0.05 & Median_R < -r_thresh)
  
  write.csv(tier_pos, paste0(prefix, "_R", r_thresh, "_Positive.csv"), row.names = FALSE)
  write.csv(tier_neg, paste0(prefix, "_R", r_thresh, "_Negative.csv"), row.names = FALSE)
  
  return(data.frame(Threshold = r_thresh, Pos_Count = nrow(tier_pos), Neg_Count = nrow(tier_neg)))
}

# Run extraction for 0.0, 0.3, and 0.5
stats_00 <- extract_tier(final_network, 0.0, output_post)
stats_03 <- extract_tier(final_network, 0.3, output_post)
stats_05 <- extract_tier(final_network, 0.5, output_post)

# --- 6. Summary Output ---
summary_table <- rbind(stats_00, stats_03, stats_05)
print(summary_table)

#  Threshold Pos_Count Neg_Count
#1       0.0      3839      2659
#2       0.3      2035       704
#3       0.5       275         0
message("Analysis complete. CSV files generated for all tiers.")

# --- 1. Define thresholds and directions ---
thresholds <- c(0.0, 0.3, 0.5)
directions <- c("Positive", "Negative")
bg_genes   <- final_network$Gene # Background is all genes tested

# --- 2. Automated Loop for Enrichment ---
for (tr in thresholds) {
  for (dir in directions) {
    
    cat(sprintf("\nProcessing Tier: R > %s | Direction: %s\n", tr, dir))
    
    # Filter genes based on current criteria
    if (dir == "Positive") {
      current_list <- final_network %>% filter(FDR < 0.05 & Median_R > tr) %>% pull(Gene)
    } else {
      current_list <- final_network %>% filter(FDR < 0.05 & Median_R < -tr) %>% pull(Gene)
    }
    
    # Check if there are enough genes to run GO
    if (length(current_list) < 10) {
      message("Skipping: Too few genes in this category.")
      next
    }
    
    # Run Enrichment (using your function)
    label <- paste0("R", tr, "_", dir)
    ego_res <- run_enrichment(current_list, bg_genes, label)
    
    # Plot results if significant
    if (!is.null(ego_res)) {
      # Choose color theme based on direction
      color_theme <- ifelse(dir == "Positive", "#4393C3", "#4393C3")
      #4393C3
      p <- plot_dual_axis_label_on_bar(
        enrich_obj  = ego_res, 
        color_theme = color_theme, 
        title_text  = paste("C4A", dir, "Seed Network (Postnatal)"), 
        top_n       = 8
      )
      
      if (!is.null(p)) {
        file_name <- sprintf("C4A_GO_Postnatal_%s_R%s.png", dir, tr)
        ggsave(file_name, p, width = 7, height = 6, dpi = 400)
        message("Saved: ", file_name)
      }
    }
  }
}

