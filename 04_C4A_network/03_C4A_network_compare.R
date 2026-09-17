suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(VennDiagram)
  library(grid)
})


load_genes <- function(filepath) {
  if(!file.exists(filepath)) {
    warning(paste("File not found:", filepath))
    return(character(0))
  }
  dt <- fread(filepath)
  gene_col <- intersect(c("Gene", "Gene_version", "GID"), colnames(dt))
  if(length(gene_col) > 0) {
    genes <- unique(gsub("\\..*", "", dt[[gene_col[1]]]))
  } else {
    genes <- unique(gsub("\\..*", "", dt[[1]]))
  }
  return(genes)
}

pos_pre <- load_genes("C4A_Positive_Network_pre03.csv")
neg_pre <- load_genes("C4A_Negative_Network_pre03.csv")
pos_post <- load_genes("C4A_Positive_Network_post03.csv")
neg_post <- load_genes("C4A_Negative_Network_post03.csv")


prenatal_fill  <- "#f4a582" # 柔和红
postnatal_fill <- "#92c5de" # 柔和蓝
prenatal_col   <- "#b2182b" # 深红边框
postnatal_col  <- "#2166ac" # 深蓝边框

# 3. Plot C4A-positive genes 

png("C4A_Positive_Network_Overlap_Venn.png", width = 2050, height = 2000, res = 350)
grid.newpage()
venn_pos <- venn.diagram(
  x = list(Prenatal = pos_pre, Postnatal = pos_post),
  category.names = c("Prenatal", "Postnatal"),
  filename = NULL,
  fill = c(prenatal_fill, postnatal_fill),
  alpha = 0.5,
  col = c(prenatal_col, postnatal_col),
  cex = 1.5,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27),
  cat.dist = c(0.05, 0.05)
)
grid.draw(venn_pos)
grid.text("C4A-positive genes", x = 0.5, y = 0.90, gp = gpar(fontsize = 16, fontface = "bold", col = "black"))
dev.off()


# 4. Plot C4A-negative genes venn 

png("C4A_Negative_Network_Overlap_Venn.png", width = 2050, height = 2000, res = 350)
grid.newpage()
venn_neg <- venn.diagram(
  x = list(Prenatal = neg_pre, Postnatal = neg_post),
  category.names = c("Prenatal", "Postnatal"),
  filename = NULL,
  fill = c(prenatal_fill, postnatal_fill),
  alpha = 0.5,
  col = c(prenatal_col, postnatal_col),
  cex = 1.5,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27),
  cat.dist = c(0.05, 0.05)
)
grid.draw(venn_neg)
grid.text("C4A-negative genes", x = 0.5, y = 0.90, gp = gpar(fontsize = 16, fontface = "bold", col = "black"))
dev.off()

cat("Venn diagrams updated and saved successfully:\n")
cat(" - C4A_Positive_Network_Overlap_Venn.png (Title: C4A-positive genes)\n")
cat(" - C4A_Negative_Network_Overlap_Venn.png (Title: C4A-negative genes)\n")



library(data.table)
library(dplyr)
library(ggplot2)
rm(list = ls())
options(stringsAsFactors = FALSE)
pre_df  <- fread("C4A_Negative_Network_pre03.csv")
post_df <- fread("C4A_Negative_Network_post03.csv")
pre_df$Gene_clean  <- gsub("\\..*", "", pre_df$Gene)
post_df$Gene_clean <- gsub("\\..*", "", post_df$Gene)

genes_pre  <- unique(pre_df$Gene_clean)
genes_post <- unique(post_df$Gene_clean)

# --- 3. Dimension 1: Gene Overlap and Fisher's Exact Test ---
overlap_genes <- intersect(genes_pre, genes_post)
union_genes   <- union(genes_pre, genes_post)

# Calculate the Jaccard similarity index
jaccard_index <- length(overlap_genes) / length(union_genes)

# Hypergeometric distribution test (assuming a background gene count of 28,146)
universe_size <- 28146 

a <- length(overlap_genes)                 # Shared by both
b <- length(genes_pre) - a                  # Only Prenatal
c <- length(genes_post) - a                 # Only Postnatal
d <- universe_size - (a + b + c)            # Neither

contingency_table <- matrix(c(a, b, c, d), nrow = 2)
fisher_res <- fisher.test(contingency_table, alternative = "greater")

# --- 4. Dimension 2: Analysis of R-value consistency for shared genes ---
shared_df <- inner_join(
  pre_df  %>% dplyr::select(Gene_clean, R_pre = R),
  post_df %>% dplyr::select(Gene_clean, R_post = R),
  by = "Gene_clean"
)

# Calculate the Spearman correlation coefficient (to evaluate the direction and relative rank consistency of the R-value)
r_cor <- cor(shared_df$R_pre, shared_df$R_post, method = "spearman")

cat("\n======================================================\n")
cat("          Co-expression Network Similarity             \n")
cat("======================================================\n")
cat("Prenatal Network Genes count  :", length(genes_pre), "\n")
cat("Postnatal Network Genes count :", length(genes_post), "\n")
cat("Overlapping Genes count        :", length(overlap_genes), "\n")
cat("Jaccard Index                  :", round(jaccard_index, 4), "\n")
cat("Fisher's Exact Test Odds Ratio :", round(fisher_res$estimate, 2), "\n")
cat("Fisher's Exact Test P-value    :", fisher_res$p.value, "\n")
cat("Shared Genes R-value Correlation:", round(r_cor, 4), "\n")
cat("======================================================\n")

# --- 6. Plot scatter plot of R-values ​​for shared genes (optional) ---
ggplot(shared_df, aes(x = R_pre, y = R_post)) +
  geom_point(alpha = 0.5, color = "steelblue") +
  geom_smooth(method = "lm", color = "darkred", se = TRUE) +
  theme_bw() +
  labs(
    title = paste0("R-value Preservation (Spearman R = ", round(r_cor, 3), ")"),
    x = "Prenatal R-value",
    y = "Postnatal R-value"
  )


# ==============================================================================
# Compare Prenatal vs Postnatal GO Enrichment Results Directly
# ==============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)


suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyverse)
  library(GOSemSim)        # GO Semantic Similarity
  library(org.Hs.eg.db)  
  library(ggplot2)         
})

# ==============================================================================
# Part 1: Reading and Preprocessing GO Result Files
# ==============================================================================
file_pre  <- "C4A_GO_Result_prenatal03_Positive.txt"
file_post <- "C4A_GO_Result_postnatal03_Positive.txt"

go_pre  <- fread(file_pre, header = TRUE, sep = "\t")
go_post <- fread(file_post, header = TRUE, sep = "\t")

# Significance filtering (default p.adjust < 0.05)
p_cutoff <- 0.05
go_pre_sig  <- go_pre  %>% filter(p.adjust < p_cutoff)
go_post_sig <- go_post %>% filter(p.adjust < p_cutoff)

# Part 2: GO Term Level Overlap

ids_pre  <- unique(go_pre_sig$ID)
ids_post <- unique(go_post_sig$ID)

common_go_ids <- intersect(ids_pre, ids_post)
union_go_ids  <- union(ids_pre, ids_post)

jaccard_go <- length(common_go_ids) / length(union_go_ids)

cat("==========================================\n")
cat("1. GO Term Overlap Statistics (p.adjust <", p_cutoff, ")\n")
cat("------------------------------------------\n")
cat("Prenatal Significant GO Count  :", length(ids_pre), "\n")
cat("Postnatal Significant GO Count :", length(ids_post), "\n")
cat("Shared (Overlap) GO Count      :", length(common_go_ids), "\n")
cat("Jaccard Index (GO Level)       :", round(jaccard_go, 4), "\n\n")

# Export detailed comparison of shared pathways
if (length(common_go_ids) > 0) {
  shared_df <- inner_join(
    go_pre_sig  %>% dplyr::select(ID, Description, GeneRatio_pre = GeneRatio, p.adjust_pre = p.adjust, Count_pre = Count),
    go_post_sig %>% dplyr::select(ID, GeneRatio_post = GeneRatio, p.adjust_post = p.adjust, Count_post = Count),
    by = "ID"
  ) %>% arrange(p.adjust_pre)
  
  write.csv(shared_df, "GO_Pos_Shared_Terms_Comparison.csv", row.names = FALSE)
  cat("Shared GO terms saved to 'GO_Shared_Terms_Comparison.csv'\n")
}

# ==============================================================================
# Part 3: GO Semantic Similarity
# ==============================================================================
cat("\n==========================================\n")
cat("2. GO Semantic Similarity Analysis (GOSemSim)\n")
cat("------------------------------------------\n")

if (length(ids_pre) > 0 && length(ids_post) > 0) {
 # Constructing a BP (Biological Process) database object
  hsGO <- godata('org.Hs.eg.db', ont = "BP")
  
  # Computing Semantic Similarity (Wang's Method + BMA Combination Strategy)
  sem_sim <- mgoSim(
    ids_pre, 
    ids_post, 
    semData = hsGO, 
    measure = "Wang", 
    combine = "BMA"
  )
  cat("GO Semantic Similarity Score (Wang/BMA):", round(sem_sim, 4), "\n")
} else {
  cat("Insufficient GO terms to calculate semantic similarity.\n")
}
#Positive" preparing gene to GO mapping data...
#preparing IC data...
#GO Semantic Similarity Score (Wang/BMA): 0.772
# ==============================================================================
# Part 4: Back-extracting enriched genes from GO results for gene-level comparison
# ==============================================================================
cat("\n==========================================\n")
cat("3. Enriched Genes Extraction & Overlap Analysis\n")
cat("------------------------------------------\n")

extract_genes <- function(df) {
  if ("geneID" %in% colnames(df)) {
    genes <- unlist(strsplit(df$geneID, "/"))
    return(unique(genes))
  } else {
    return(character(0))
  }
}

genes_in_pre_go  <- extract_genes(go_pre_sig)
genes_in_post_go <- extract_genes(go_post_sig)

common_genes_in_go <- intersect(genes_in_pre_go, genes_in_post_go)

cat("Genes driving Prenatal GOs  :", length(genes_in_pre_go), "\n")
cat("Genes driving Postnatal GOs :", length(genes_in_post_go), "\n")
cat("Shared Genes in GO terms    :", length(common_genes_in_go), "\n")
if (length(common_genes_in_go) > 0) {
  cat("Shared Gene List            :", paste(common_genes_in_go, collapse = ", "), "\n")
}

#Genes driving Prenatal GOs  : 77 
#Genes driving Postnatal GOs : 237 
#Shared Genes in GO terms    : 10 
#Shared Gene List            : CYFIP2, EPHB2, NTNG1, APBB1, KCNN1, KCNH3, KCNA2, DAGLA, TTBK1, BASP1 

cat("\n==========================================\n")
cat("4. Generating Comparison Visualization Plot...\n")
cat("------------------------------------------\n")

top_pre  <- go_pre_sig  %>% arrange(p.adjust) %>% head(10) %>% mutate(Group = "Prenatal")
top_post <- go_post_sig %>% arrange(p.adjust) %>% head(10) %>% mutate(Group = "Postnatal")

combined_top <- bind_rows(
  top_pre  %>% dplyr::select(ID, Description, p.adjust, Count, Group),
  top_post %>% dplyr::select(ID, Description, p.adjust, Count, Group)
)

# Ensure data from both sides is present in the comparison plot (fill in missing values ​​to enable the plotting of comparison bubbles)
all_selected_ids <- unique(combined_top$ID)

plot_df_pre <- go_pre %>% 
  filter(ID %in% all_selected_ids) %>% 
  mutate(Group = "Prenatal") %>% 
  dplyr::select(ID, Description, p.adjust, Count, Group)

plot_df_post <- go_post %>% 
  filter(ID %in% all_selected_ids) %>% 
  mutate(Group = "Postnatal") %>% 
  dplyr::select(ID, Description, p.adjust, Count, Group)

plot_df <- bind_rows(plot_df_pre, plot_df_post) %>%
  mutate(
    Log10P = -log10(p.adjust),
    Log10P = ifelse(Log10P > 20, 20, Log10P) 
  )

# Dotplot
p <- ggplot(plot_df, aes(x = Group, y = reorder(Description, Log10P))) +
  geom_point(aes(size = Count, color = Log10P)) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(p.adjust)") +
  scale_size_continuous(range = c(3, 8), name = "Gene Count") +
  theme_bw() +
  labs(
    x = "Developmental Stage",
    y = "GO Term (Biological Process)",
    title = "Functional Pathway Comparison(C4A-Negtivate): Prenatal vs Postnatal"
  ) +
  theme(
    axis.text.x  = element_text(size = 11, face = "bold"),
    axis.text.y  = element_text(size = 10),
    plot.title   = element_text(hjust = 0.5, face = "bold")
  )

ggsave("GO_Comparison_Dotplot.pdf", plot = p, width = 9, height = 7)
cat("Comparison plot saved as 'GO_Comparison_Dotplot.pdf'\n")
cat("==========================================\n")
#2. GO Semantic Similarity Analysis (GOSemSim)
------------------------------------------
#preparing gene to GO mapping data...
#preparing IC data...
#GO Semantic Similarity Score (Wang/BMA): 0.521 


# ==============================================================================
# GO Pathways Comparison: Shared vs Prenatal-Specific vs Postnatal-Specific
# ==============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyverse)
  library(ggplot2)
})

# 1. Read GO result file
file_pre  <- "C4A_GO_Result_prenatal03_Negative.txt"
file_post <- "C4A_GO_Result_postnatal03_Negative.txt"

go_pre  <- fread(file_pre, header = TRUE, sep = "\t")
go_post <- fread(file_post, header = TRUE, sep = "\t")

# 2. Set significance threshold and filter
p_cutoff <- 0.05

go_pre_sig  <- go_pre  %>% filter(p.adjust < p_cutoff)
go_post_sig <- go_post %>% filter(p.adjust < p_cutoff)

ids_pre  <- go_pre_sig$ID
ids_post <- go_post_sig$ID

# 3. Three-Class Set Extraction
shared_ids    <- intersect(ids_pre, ids_post)
pre_only_ids  <- setdiff(ids_pre, ids_post)
post_only_ids <- setdiff(ids_post, ids_pre)


cat("==========================================\n")
cat("Pathway Classification Summary (p.adjust <", p_cutoff, ")\n")
cat("------------------------------------------\n")
cat("Total Prenatal Significant GOs  :", length(ids_pre), "\n")
cat("Total Postnatal Significant GOs :", length(ids_post), "\n")
cat("1. Shared GO Pathways           :", length(shared_ids), "\n")
cat("2. Prenatal-Specific GOs        :", length(pre_only_ids), "\n")
cat("3. Postnatal-Specific GOs       :", length(post_only_ids), "\n")
cat("==========================================\n\n")
#==========================================
#Pathway Classification Summary (p.adjust < 0.05 )
#------------------------------------------
#Total Prenatal Significant GOs  : 50 
#Total Postnatal Significant GOs : 145 
#1. Shared GO Pathways           : 13 
#2. Prenatal-Specific GOs        : 37 
#3. Postnatal-Specific GOs       : 132 
#==========================================


# 4.1 Shared Pathway Table
df_shared <- inner_join(
  go_pre_sig  %>% dplyr::select(ID, Description, GeneRatio_pre = GeneRatio, p.adjust_pre = p.adjust, Count_pre = Count, geneID_pre = geneID),
  go_post_sig %>% dplyr::select(ID, GeneRatio_post = GeneRatio, p.adjust_post = p.adjust, Count_post = Count, geneID_post = geneID),
  by = c("ID")
) %>% arrange(p.adjust_pre)

# 4.2 Table of Prenatal-Specific Pathways
df_pre_specific <- go_pre_sig %>% 
  filter(ID %in% pre_only_ids) %>% 
  arrange(p.adjust)

# 4.3 Table of Postnatal-Specific Pathways
df_post_specific <- go_post_sig %>% 
  filter(ID %in% post_only_ids) %>% 
  arrange(p.adjust)


write.csv(df_shared,        "GO_Shared_Pathways.csv",        row.names = FALSE)
write.csv(df_pre_specific,  "GO_Prenatal_Specific_Pathways.csv", row.names = FALSE)
write.csv(df_post_specific, "GO_Postnatal_Specific_Pathways.csv", row.names = FALSE)

# 5. Visualization: Select the top 8 pathways for each category to generate a faceted dotplot for comparison.

top_shared <- df_shared %>% 
  head(8) %>% 
  mutate(Category = "Shared") %>% 
  dplyr::select(ID, Description, p.adjust = p.adjust_pre, Count = Count_pre, Category)

top_pre_spec <- df_pre_specific %>% 
  head(8) %>% 
  mutate(Category = "Prenatal Specific") %>% 
  dplyr::select(ID, Description, p.adjust, Count, Category)

top_post_spec <- df_post_specific %>% 
  head(8) %>% 
  mutate(Category = "Postnatal Specific") %>% 
  dplyr::select(ID, Description, p.adjust, Count, Category)

plot_df <- bind_rows(top_shared, top_pre_spec, top_post_spec) %>%
  mutate(
    Log10P = -log10(p.adjust),
    Category = factor(Category, levels = c("Prenatal Specific", "Shared", "Postnatal Specific"))
  )

p <- ggplot(plot_df, aes(x = Log10P, y = reorder(Description, Log10P))) +
  geom_point(aes(size = Count, color = Log10P)) +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(p.adjust)") +
  scale_size_continuous(range = c(3, 7), name = "Gene Count") +
  theme_bw() +
  labs(
    x = "-log10(p.adjust)",
    y = "GO Term (Biological Process)",
    title = "C4A-seeded negative network_GO Pathways Comparison"
  ) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold", size = 12),
    axis.text.y      = element_text(size = 12),
    plot.title       = element_text(hjust = 0.5, face = "bold")
  )

ggsave("GO_negtiveSpecific_vs_Shared_Dotplot.pdf", plot = p, width = 8, height = 6)
ggsave("GO_negtiveSpecific_vs_Shared_Dotplot.png", plot = p, width = 8, height = 6,dpi=500)
message("Saved visualization plot as 'GO_Specific_vs_Shared_Dotplot.pdf'")


