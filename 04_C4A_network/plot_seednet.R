suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(igraph)
  library(ggraph)      
  library(tidygraph)   
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

# Load the expression matrix (used to calculate the internal spiderweb connections between target genes)

expr_file <- "/gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/02_C4_network/networkv2/prenatal_expression.residual.txt"
expr_data <- fread(expr_file, header = TRUE)
colnames(expr_data)[1] <- "GID"

expr_mat <- expr_data %>%
  column_to_rownames("GID") %>%
  as.matrix() %>%
  t()

colnames(expr_mat) <- gsub("\\..*", "", colnames(expr_mat))
seed_gene <- "ENSG00000244731"  # C4A
if(!seed_gene %in% colnames(expr_mat)) stop("Seed gene not found in expression matrix!")

# Define plotting function based on CSV results (excluding unnamed ENSG genes)

plot_single_module_from_csv <- function(nodes_csv, seed_name, direction_type, color_theme, edge_color, expr_matrix) {
  nodes_df <- fread(nodes_csv, header = TRUE) %>% as.data.frame()
  if("R" %in% colnames(nodes_df)) {
    nodes_df$Cor <- nodes_df$R
  } else if("R" %in% colnames(nodes_df)) {
    nodes_df$Cor <- nodes_df$R
  } else {
    stop("Error: Cannot find 'R' or 'Median_R' column in CSV.")
  }
  nodes_df$Gene_clean <- gsub("\\..*", "", nodes_df$Gene)
  gene_symbols <- suppressMessages(mapIds(org.Hs.eg.db, keys = nodes_df$Gene_clean, 
                                          column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"))
  
  nodes_df$SYMBOL <- gene_symbols

  nodes_df <- nodes_df %>% 
    filter(!is.na(SYMBOL)) %>%                     
    filter(!grepl("^ENSG", SYMBOL))              
  
# From the filtered set of known genes, select the top 50 most significant ones.
  nodes_df <- nodes_df %>% 
    arrange(FDR) %>% 
    head(60)
  #desc(abs(Cor));arrange(FDR)
  if(nrow(nodes_df) == 0) {
    cat("Warning: No valid named genes left after filtering.\n")
    return(NULL)
  }
  
 # 6. Preparing Network Nodes
  target_genes <- c(seed_name, nodes_df$Gene_clean)
  seed_info <- data.frame(Gene_clean = seed_name, Cor = 1, SYMBOL = "C4A")
  
  nodes <- bind_rows(seed_info, nodes_df %>% dplyr::select(Gene_clean, Cor, SYMBOL)) %>%
    mutate(
      Type = ifelse(Gene_clean == seed_name, "Seed", direction_type),
      Abs_Cor = abs(Cor)
    )
  
  # 7. Extract internal interaction edges of target genes (Internal Edges)
  available_genes <- intersect(target_genes, colnames(expr_matrix))
  expr_sub <- expr_matrix[, available_genes, drop = FALSE] 
  
  adj_mat <- cor(expr_sub, method = "pearson")
  adj_mat[lower.tri(adj_mat, diag = TRUE)] <- 0 
  
  internal_edges <- as.data.frame(adj_mat) %>%
    rownames_to_column(var = "from") %>%
    pivot_longer(cols = -from, names_to = "to", values_to = "Cor") %>%
    filter(abs(Cor) > 0.1 & Cor != 0) %>% 
    mutate(from = as.character(from), to = as.character(to))
  
 # 8. Preparing the connection from the seed to the target gene
  seed_edges <- data.frame(from = seed_name, to = nodes_df$Gene_clean, Cor = nodes_df$Cor)
  
  final_edges <- bind_rows(internal_edges, seed_edges) %>%
    filter(from %in% nodes$Gene_clean & to %in% nodes$Gene_clean) %>%
    distinct(from, to, .keep_all = TRUE) %>%
    mutate(weight = abs(Cor)) 

  graph <- tbl_graph(nodes = nodes, edges = final_edges, directed = FALSE)
  
  p <- ggraph(graph, layout = 'fr') +
    geom_edge_link(aes(edge_alpha = weight, edge_width = weight), 
                   color = edge_color, show.legend = FALSE) +
    scale_edge_width_continuous(range = c(0.2, 1.5)) +
    scale_edge_alpha_continuous(range = c(0.3, 0.8)) +
    
    geom_node_point(aes(color = Type, size = Abs_Cor), stroke = 0.5) +
    scale_color_manual(values = setNames(c("black", color_theme), c("Seed", direction_type)), 
                       guide = "none") +
    scale_size_continuous(range = c(3, 12),breaks = function(x) { pretty(x, n = 3) }) +
    geom_node_text(aes(label = SYMBOL), repel = TRUE, size = 3.8, fontface = "plain",
                   bg.color = "white", bg.r = 0.15, max.overlaps = Inf) +
    
    theme_graph(base_family = "sans") +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, size = 16, face = "plain")) +
    labs(title = paste("C4A-Seeded", direction_type, "Network (Prenatal)"),
         size = "Absolute Correlation")
  
  return(p)
}

# ==============================================================================
# 3. 运行并生成图片 (以 Prenatal 为例)
# ==============================================================================
cat("Plotting network modules from CSV results...\n")

# --- 1. 画产前正相关网络 ---
p_pos <- plot_single_module_from_csv(
  nodes_csv = "/gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/02_C4_network/networkv2/C4A_Positive_Network_pre03.csv", 
  seed_name = seed_gene, 
  direction_type = "Positive", 
  color_theme = "#e89a92",       
  edge_color = "#f4cbc7",        
  expr_matrix = expr_mat
)

if(!is.null(p_pos)){
  ggsave("C4A_Seed_Network_Positive_Prenatal_New.png", p_pos, width =6, height =6, dpi = 300)
  cat("Saved Positive Network Plot.\n")
}

# --- 2. 画产前负相关网络 ---
p_neg <- plot_single_module_from_csv(
  nodes_csv = "C4A_Negative_Network_pre.csv", 
  seed_name = seed_gene, 
  direction_type = "Negative", 
  color_theme = "#e89a92",       
  edge_color = "#f4cbc7",      
  expr_matrix = expr_mat
)

if(!is.null(p_neg)){
  ggsave("C4A_Seed_Network_Negative_Prenatal_New.png", p_neg, width = 6, height = 6, dpi = 300)
  cat("Saved Negative Network Plot.\n")
}


###############################################Postnatal##################################################################



# "#83bcd7","#c8e2f0"
# "#e89a92","#f4cbc7"
library(tidyverse)
library(data.table)
library(igraph)
library(ggraph)      
library(tidygraph)   
library(clusterProfiler)
library(org.Hs.eg.db)

expr_file <- "/gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/02_C4_network/networkv2/postnatal_expression.residual.txt"
expr_data <- fread(expr_file, header = TRUE)
colnames(expr_data)[1] <- "GID"


expr_mat <- expr_data %>%
  column_to_rownames("GID") %>%
  as.matrix() %>%
  t()

colnames(expr_mat) <- gsub("\\..*", "", colnames(expr_mat))
seed_gene <- "ENSG00000244731"  # C4A
if(!seed_gene %in% colnames(expr_mat)) stop("Seed gene not found in expression matrix!")

# Define plotting function based on CSV results (excluding unnamed ENSG genes)

plot_single_module_from_csv <- function(nodes_csv, seed_name, direction_type, color_theme, edge_color, expr_matrix) {
  nodes_df <- fread(nodes_csv, header = TRUE) %>% as.data.frame()
  if("R" %in% colnames(nodes_df)) {
    nodes_df$Cor <- nodes_df$R
  } else if("R" %in% colnames(nodes_df)) {
    nodes_df$Cor <- nodes_df$R
  } else {
    stop("Error: Cannot find 'R' or 'Median_R' column in CSV.")
  }
  nodes_df$Gene_clean <- gsub("\\..*", "", nodes_df$Gene)
  gene_symbols <- suppressMessages(mapIds(org.Hs.eg.db, keys = nodes_df$Gene_clean, 
                                          column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"))
  
  nodes_df$SYMBOL <- gene_symbols

  nodes_df <- nodes_df %>% 
    filter(!is.na(SYMBOL)) %>%                     
    filter(!grepl("^ENSG", SYMBOL))              
  
# From the filtered set of known genes, select the top 50 most significant ones.
  nodes_df <- nodes_df %>% 
    arrange(FDR) %>% 
    head(60)
  #desc(abs(Cor));arrange(FDR)
  if(nrow(nodes_df) == 0) {
    cat("Warning: No valid named genes left after filtering.\n")
    return(NULL)
  }
  
 # 6. Preparing Network Nodes
  target_genes <- c(seed_name, nodes_df$Gene_clean)
  seed_info <- data.frame(Gene_clean = seed_name, Cor = 1, SYMBOL = "C4A")
  
  nodes <- bind_rows(seed_info, nodes_df %>% dplyr::select(Gene_clean, Cor, SYMBOL)) %>%
    mutate(
      Type = ifelse(Gene_clean == seed_name, "Seed", direction_type),
      Abs_Cor = abs(Cor)
    )
  
  # 7. Extract internal interaction edges of target genes (Internal Edges)
  available_genes <- intersect(target_genes, colnames(expr_matrix))
  expr_sub <- expr_matrix[, available_genes, drop = FALSE] 
  
  adj_mat <- cor(expr_sub, method = "pearson")
  adj_mat[lower.tri(adj_mat, diag = TRUE)] <- 0 
  
  internal_edges <- as.data.frame(adj_mat) %>%
    rownames_to_column(var = "from") %>%
    pivot_longer(cols = -from, names_to = "to", values_to = "Cor") %>%
    filter(abs(Cor) > 0.1 & Cor != 0) %>% 
    mutate(from = as.character(from), to = as.character(to))
  
 # 8. Preparing the connection from the seed to the target gene
  seed_edges <- data.frame(from = seed_name, to = nodes_df$Gene_clean, Cor = nodes_df$Cor)
  
  final_edges <- bind_rows(internal_edges, seed_edges) %>%
    filter(from %in% nodes$Gene_clean & to %in% nodes$Gene_clean) %>%
    distinct(from, to, .keep_all = TRUE) %>%
    mutate(weight = abs(Cor)) 

  graph <- tbl_graph(nodes = nodes, edges = final_edges, directed = FALSE)
  
  p <- ggraph(graph, layout = 'fr') +
    geom_edge_link(aes(edge_alpha = weight, edge_width = weight), 
                   color = edge_color, show.legend = FALSE) +
    scale_edge_width_continuous(range = c(0.2, 1.5)) +
    scale_edge_alpha_continuous(range = c(0.3, 0.8)) +
    
    geom_node_point(aes(color = Type, size = Abs_Cor), stroke = 0.5) +
    scale_color_manual(values = setNames(c("black", color_theme), c("Seed", direction_type)), 
                       guide = "none") +
    scale_size_continuous(range = c(3, 12),breaks = function(x) { pretty(x, n = 3) }) +
    geom_node_text(aes(label = SYMBOL), repel = TRUE, size = 3.8, fontface = "plain",
                   bg.color = "white", bg.r = 0.15, max.overlaps = Inf) +
    
    theme_graph(base_family = "sans") +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, size = 16, face = "plain")) +
    labs(title = paste("C4A-Seeded", direction_type, "Network (Postnatal)"),
         size = "Absolute Correlation")
  
  return(p)
}

# ==============================================================================
# 3. 运行并生成图片 (以 Prenatal 为例)
# ==============================================================================
cat("Plotting network modules from CSV results...\n")

# --- 1. 画产前正相关网络 ---
p_pos <- plot_single_module_from_csv(
  nodes_csv = "/gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/02_C4_network/networkv2/C4A_Positive_Network_post03.csv", 
  seed_name = seed_gene, 
  direction_type = "Positive", 
  color_theme = "#83bcd7",       
  edge_color = "#c8e2f0",        
  expr_matrix = expr_mat
)

if(!is.null(p_pos)){
  ggsave("C4A_Seed_Network_Positive_Postnatal_New.png", p_pos, width =6, height =6, dpi = 300)
  cat("Saved Positive Network Plot.\n")
}

# --- 2. 画产前负相关网络 ---
p_neg <- plot_single_module_from_csv(
  nodes_csv = "/gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/02_C4_network/networkv2/C4A_Negative_Network_post03.csv", 
  seed_name = seed_gene, 
  direction_type = "Negative", 
  color_theme = "#83bcd7",       
  edge_color = "#c8e2f0",      
  expr_matrix = expr_mat
)

if(!is.null(p_neg)){
  ggsave("C4A_Seed_Network_Negative_Postnatal_New.png", p_neg, width = 6, height = 6, dpi = 300)
  cat("Saved Negative Network Plot.\n")
}




