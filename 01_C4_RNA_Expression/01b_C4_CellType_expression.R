# C4_CellType_Analysis.R: Paralog Cell-Type Specific Expression 
# This script performs deconvolution-based cell-type correlation, 
#'proportion-corrected expression estimation, sex-bias analysis for C4A and C4B in human brain.

### Search by
# 1a. correlation of C4 exp and prenatal each cell-type correlation: scatter plot
# 1b. correlation of C4 exp and postnatal each cell-type correlation
# 2a: Barplot: correltaion of C4 gene expression and total cell proportions of prenatal and postnatal.
# 2b: Barplot: correltaion of C4 gene expression and core cell (4 types) proportions of prenatal and postnatal.
# 3. c4 gene cell type proportion-corrected expression estimation (each cell type)
# 

# 1a. correlation of C4 exp and prenatal each cell-type correlation: scatter plot
setwd("/gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/figures/fig1")
prop_file <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/pre_devbulk_prop.csv" 
expr_file <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/bulk_data/stage1.logTMM.ComBat.txt"

target_genes <- c("C4A" = "ENSG00000244731.9", 
                  "C4B" = "ENSG00000224389.9")
props <- read.csv(prop_file, row.names = 1, check.names = FALSE)
props$sample <- rownames(props)
print("Cell Types found:")
print(colnames(props))

expr_data <- fread(expr_file, header = TRUE)
colnames(expr_data)[1] <- "GID"
expr_c4 <- expr_data %>% 
  filter(GID %in% target_genes) %>%
  column_to_rownames("GID") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample")

if(target_genes["C4A"] %in% colnames(expr_c4)) colnames(expr_c4)[colnames(expr_c4) == target_genes["C4A"]] <- "C4A"
if(target_genes["C4B"] %in% colnames(expr_c4)) colnames(expr_c4)[colnames(expr_c4) == target_genes["C4B"]] <- "C4B"

df_merged <- merge(expr_c4, props, by = "sample")

cell_types <- colnames(props)[colnames(props) != "sample"]
df_long <- df_merged %>%
  pivot_longer(cols = all_of(cell_types), 
               names_to = "CellType", 
               values_to = "Proportion") %>%
  pivot_longer(cols = c("C4A", "C4B"),
               names_to = "Gene",
               values_to = "Expression")


# Scatter plot of correlations among key cell types (focusing on MG, AST, and ExNeu)

#target_cells <- c("GLIALPROG","prenatal_MG", "prenatal_AST", "prenatal_ExNeu", "prenatal_IN", "prenatal_OPC")
target_cells <- c("GLIALPROG","prenatal_MG", "prenatal_AST", "prenatal_ExNeu", "prenatal_IN", "prenatal_OPC")
plot_subset <- df_long %>% filter(CellType %in% target_cells)
my_palette <- c("#CB6463", "#889ABF", "#A9D69C", "#BB9BCA", "#E59194", "#9BCCEA")
df_long_prenatal<- plot_subset
p_scatter <- ggplot(plot_subset, aes(x = Proportion, y = Expression)) +
  geom_point(aes(color = CellType), alpha = 0.6, size = 1.5) +
  geom_smooth(method = "lm", color = "grey30", size = 0.5, se = TRUE) + 
  
# Faceting: Row = C4 gene, Column = Cell type
  facet_grid(Gene ~ CellType, scales = "free") +
  stat_cor(method = "pearson", 
           label.x.npc = "left", 
           label.y.npc = "top", 
           size = 4, 
           color = "black") +

  theme_bw(base_size = 19) + 
  scale_color_manual(values = my_palette) +
  labs(title = "Correlation: C4 vs Cell Type (Prenatal Brain)",
       x = "Cell Type Proportion",
       y = "Bulk Expression (LogTMM)") +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"), 
        plot.title = element_text(hjust = 0.5, size = 19)) 

print(p_scatter)
ggsave("C4_CellType_Scatter_Prenatal_LargeFont.png", p_scatter, width = 12, height = 6)


# 1b. postnatal cell-type correlation
prop_file <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/post_devbulk_prop.csv" 
expr_file <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/bulk_data/postnatal.logTMM.ComBat.txt"

target_genes <- c("C4A" = "ENSG00000244731.9", 
                  "C4B" = "ENSG00000224389.9")

props <- read.csv(prop_file, row.names = 1, check.names = FALSE)
props$sample <- rownames(props)
print("Cell Types found:")
print(colnames(props))

expr_data <- fread(expr_file, header = TRUE)
colnames(expr_data)[1] <- "GID"

expr_c4 <- expr_data %>% 
  filter(GID %in% target_genes) %>%
  column_to_rownames("GID") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample")

if(target_genes["C4A"] %in% colnames(expr_c4)) colnames(expr_c4)[colnames(expr_c4) == target_genes["C4A"]] <- "C4A"
if(target_genes["C4B"] %in% colnames(expr_c4)) colnames(expr_c4)[colnames(expr_c4) == target_genes["C4B"]] <- "C4B"

df_merged <- merge(expr_c4, props, by = "sample")
cell_types <- colnames(props)[colnames(props) != "sample"] 
df_long <- df_merged %>%
  pivot_longer(cols = all_of(cell_types), 
               names_to = "CellType", 
               values_to = "Proportion") %>%
  pivot_longer(cols = c("C4A", "C4B"),
               names_to = "Gene",
               values_to = "Expression")


target_cells <- c("OL","postnatal_MG", "postnatal_AST", "postnatal_ExNeu", "postnatal_IN", "postnatal_OPC")
plot_subset <- df_long %>% filter(CellType %in% target_cells)
df_long_postnatal<- plot_subset
my_palette <- c("#CB6463", "#889ABF", "#A9D69C", "#BB9BCA", "#E59194", "#9BCCEA")

p_scatter <- ggplot(plot_subset, aes(x = Proportion, y = Expression)) +
  geom_point(aes(color = CellType), alpha = 0.6, size = 1.5) +
  geom_smooth(method = "lm", color = "black", size = 0.5, se = TRUE) + 
  facet_grid(Gene ~ CellType, scales = "free") +
 
  stat_cor(method = "pearson", 
           label.x.npc = "left", 
           label.y.npc = "top", 
           size = 4, #Increase the font size to 7 (approximately 20pt).
           color = "black") +
  
  theme_bw(base_size = 19) + 
  scale_color_manual(values = my_palette) +
  labs(title = "Correlation: C4 vs Cell Type (Postnatal Brain)",
       x = "Cell Type Proportion",
       y = "Bulk Expression (LogTMM)") +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"), 
        plot.title = element_text(hjust = 0.5, size = 19)) 

print(p_scatter)
ggsave("C4_CellType_Scatter_Postnatal_LargeFont.png", p_scatter, width = 13, height = 6)

# 2a: Barplot: correltaion of C4 gene expression and total cell proportions of prenatal and postnatal.

library(tidyverse)
library(ggpubr)

my_palette <- c(
  "prenatal_AST" = "#889ABF", "postnatal_AST" = "#889ABF", "AST" = "#889ABF",
  "prenatal_ExNeu" = "#A9D69C", "postnatal_ExNeu" = "#A9D69C", "ExNeu" = "#A9D69C",
  "prenatal_IN" = "#BB9BCA", "postnatal_IN" = "#BB9BCA", "IN" = "#BB9BCA",
  "prenatal_MG" = "#E59194", "postnatal_MG" = "#E59194", "MG" = "#E59194",
  "prenatal_OPC" = "#9BCCEA", "postnatal_OPC" = "#9BCCEA", "OPC" = "#9BCCEA",
  "prenatal_GLIALPROG" = "#CB6463", "postnatal_OL" = "#CB6463", "OL" = "#CB6463" 
)

get_cor_stats <- function(df_long, stage_label) {
  df_long %>%
    group_by(Gene, CellType) %>%
    summarise(
      R = cor(Proportion, Expression, method = "pearson"),
      P = cor.test(Proportion, Expression)$p.value,
      .groups = "drop"
    ) %>%
    mutate(Stage = stage_label)
}

stats_pre <- get_cor_stats(df_long_prenatal, "Prenatal")
stats_post <- get_cor_stats(df_long_postnatal, "Postnatal")
all_stats <- bind_rows(stats_pre, stats_post)

plot_bar_data <- all_stats %>%
  mutate(
    Cell_Label = gsub("prenatal_|postnatal_", "", CellType),
    Stars = case_when(
      P < 0.001 ~ "***",
      P < 0.01  ~ "**",
      P < 0.05  ~ "*",
      TRUE      ~ ""
    ),
    
# 3. Determine the Y-axis position of the asterisk (to avoid obstructing the bar chart)
# If it's a positive correlation, the position is at R + 0.02
# If it's a negative correlation, the position is at R - 0.02

    Label_Y = ifelse(R > 0, R + 0.02, R - 0.02),

    Stage = factor(Stage, levels = c("Prenatal", "Postnatal")),
    Gene = factor(Gene, levels = c("C4A", "C4B"))
  )

plot_bar_sorted <- all_stats %>%
  mutate(
    Cell_Label = gsub("prenatal_|postnatal_", "", CellType),

    Stars = case_when(
      P < 0.001 ~ "***",
      P < 0.01  ~ "**",
      P < 0.05  ~ "*",
      TRUE      ~ ""
    ),

    Label_Y = ifelse(R > 0, R + 0.02, R - 0.02),

    Stage = factor(Stage, levels = c("Prenatal", "Postnatal")),
    Gene = factor(Gene, levels = c("C4A", "C4B"))
  )
# We sorted them from highest to lowest by "average correlation" for each cell type.
cell_order <- plot_bar_sorted %>%
  group_by(Cell_Label) %>%
  summarise(Mean_R = mean(R, na.rm = TRUE)) %>%
  arrange(desc(Mean_R)) %>% # From most positive -> 0 -> most negative
  pull(Cell_Label)

plot_bar_sorted$Cell_Label <- factor(plot_bar_sorted$Cell_Label, levels = cell_order)


p_bar_sorted <- ggplot(plot_bar_sorted, aes(x = Cell_Label, y = R, fill = Gene, alpha = Stage)) +

  geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(aes(y = Label_Y, label = Stars), 
            position = position_dodge(width = 0.8), 
            vjust = ifelse(plot_bar_sorted$R > 0, -0.5, 1.2), 
            size = 4, fontface = "bold", color = "black") +
  scale_fill_manual(
    values = c("C4A" = "#d25756", "C4B" = "#7eb4db"),
    name = "Gene"
  ) +

  scale_alpha_manual(
    values = c("Prenatal" = 0.5, "Postnatal" = 1),
    name = "Development Stage"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 15),
    plot.title = element_text(hjust = 0.5, size =18),
    panel.grid.major.x = element_blank(), 
    panel.grid.major.y = element_line(color = "grey90", linetype = "dashed"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11)
  ) +
  labs(title = "C4 Expression Correlation with Cell type Proportion",
       x = NULL,
       y = "Pearson Correlation (R)") +
  scale_y_continuous(expand = expansion(mult = 0.15))

print(p_bar_sorted)

ggsave("C4_Correlation_BarPlot_Sorted.png", p_bar_sorted, width = 11, height = 5, dpi = 400)


# 2b: Barplot: correltaion of C4 gene expression and core cell (4 types) proportions of prenatal and postnatal.

target_cells_raw <- c(
  "prenatal_MG", "postnatal_MG",
  "prenatal_AST", "postnatal_AST",
  "prenatal_IN", "postnatal_IN",
  "prenatal_ExNeu", "postnatal_ExNeu"
)

plot_bar_final <- all_stats %>%
  filter(CellType %in% target_cells_raw) %>%
  mutate(
    Cell_Label = gsub("prenatal_|postnatal_", "", CellType),

    Stars = case_when(
      P < 0.001 ~ "***",
      P < 0.01  ~ "**",
      P < 0.05  ~ "*",
      TRUE      ~ ""
    ),
    Label_Y = ifelse(R > 0, R + 0.03, R - 0.03),

    Stage = factor(Stage, levels = c("Prenatal", "Postnatal")),
    Gene = factor(Gene, levels = c("C4A", "C4B"))
  )


cell_order <- plot_bar_final %>%
  group_by(Cell_Label) %>%
  summarise(Mean_R = mean(R, na.rm = TRUE)) %>%
  arrange(desc(Mean_R)) %>%
  pull(Cell_Label)


plot_bar_final$Cell_Label <- factor(plot_bar_final$Cell_Label, levels = cell_order)


p_bar_final <- ggplot(plot_bar_final, aes(x = Cell_Label, y = R, fill = Gene, alpha = Stage)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
  #geom_hline(yintercept = c(-0.4, 0.4), color = "grey85", linetype = "dashed") +
  
  geom_col(position = position_dodge(width = 0.8), width = 0.75, color = "white", linewidth = 0.2) +
  

  geom_text(aes(y = Label_Y, label = Stars), 
            position = position_dodge(width = 0.8), 
            vjust = 0.5, 
            size = 4, fontface = "bold", color = "black", show.legend = FALSE) +
  
  scale_fill_manual(
    values = c("C4A" = "#d25756", "C4B" = "#7eb4db"),
    name = "Gene"
  ) +
  
# Transparency Mapping (Developmental Stage) - Prenatally more transparent, postnatally more solid, for visual differentiation.
  scale_alpha_manual(
    values = c("Prenatal" = 0.6, "Postnatal" = 1),
    name = "Stage"
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(color = "black", size = 14, face = "bold"),
    axis.text.y = element_text(color = "black", size = 12),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.grid.major.x = element_blank(), 
    panel.grid.major.y = element_blank(),
    axis.line.y = element_line(color = "grey50"),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11)
  ) +
  labs(title = "C4 Expression Correlation with Core Cell Types",
       x = NULL,
       y = "Pearson Correlation (R)") 

print(p_bar_final)

ggsave("C4_Correlation_BarPlot_Final_CoreCells.png", p_bar_final, width = 10, height = 4, dpi = 400)


# 3. proportion-corrected expression estimation
library(tidyverse)
library(data.table)
library(ggpubr)
library(reshape2)
library(ppcor)
library(clusterProfiler)
library(org.Hs.eg.db)

theme_set(theme_bw(base_size = 17))
setwd("/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/removeversionIDmapping/cts_expr") 
files <- list.files(pattern = "^stage1_.*\\.txt$")
target_genes <- c("C4A", "C4B")
my_palette <- c(
  "prenatal_AST" = "#889ABF", "postnatal_AST" = "#889ABF", "AST" = "#889ABF",
  "prenatal_ExNeu" = "#A9D69C", "postnatal_ExNeu" = "#A9D69C", "ExNeu" = "#A9D69C",
  "prenatal_IN" = "#BB9BCA", "postnatal_IN" = "#BB9BCA", "IN" = "#BB9BCA",
  "prenatal_MG" = "#E59194", "postnatal_MG" = "#E59194", "MG" = "#E59194",
  "prenatal_OPC" = "#9BCCEA", "postnatal_OPC" = "#9BCCEA", "OPC" = "#9BCCEA",
  "prenatal_GLIALPROG" = "#CB6463", "postnatal_OL" = "#CB6463", "OL" = "#CB6463",
  "GLIALPROG" = "#CB6463" 
)

prop_file <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/pre_devbulk_prop.csv" 
prop_raw <- fread(prop_file, header = TRUE)
colnames(prop_raw)[1] <- "Sample"
prop_long <- prop_raw %>%
  pivot_longer(cols = -Sample, names_to = "CellType", values_to = "Proportion") %>%
  mutate(CellType = as.character(CellType)) 

data_list <- list()

for (f in files) {
  cell_type <- gsub("stage1_|\\.txt", "", f)
  dt <- fread(f, header = TRUE)
  colnames(dt)[1] <- "Gene"
  dt_subset <- dt %>% filter(Gene %in% target_genes)
  if (nrow(dt_subset) > 0) {
    dt_long <- dt_subset %>%
      pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Expression") %>%
      mutate(CellType = cell_type)
    
    data_list[[cell_type]] <- dt_long
  }
}
final_df <- bind_rows(data_list)


# Data merging, mathematical ratio correction

plot_df <- final_df %>% 
  filter(!grepl("VASC", CellType)) %>%
  inner_join(prop_long, by = c("Sample", "CellType")) %>%
  filter(Proportion > 0.05) %>%
  mutate(
    # 1. Convert back to linear space
    Linear_Expr = 2^Expression - 1,
    # 2. Perform scaling correction
    Corrected_Linear = Linear_Expr / Proportion,
    #3. Reset negative values ​​to 0 to prevent NaN
    Corrected_Linear = pmax(Corrected_Linear, 0),
    # 4. Perform another log transformation for visualization
    Corrected_LogExpr = log2(Corrected_Linear + 1),
    CellType = gsub("prenatal_|postnatal_", "", CellType)
  )

# Convert CellType to factors and sort by median expression level from highest to lowest (.desc = TRUE)
# Add na.rm = TRUE to prevent missing values ​​from disrupting the sorting logic.

plot_df <- plot_df %>%
  mutate(CellType = fct_reorder(CellType, Corrected_Linear, .fun = median, .desc = TRUE))
#na.rm = TRUE,
print("Check formatted and ORDERED Cell Types for plotting:")
print(levels(plot_df$CellType)) 

p <- ggplot(plot_df, aes(x = CellType, y = Corrected_LogExpr, fill = CellType)) +
  geom_boxplot(alpha = 0.8, color = "black") + 
  facet_wrap(~Gene, scales = "free_y") + 
  scale_fill_manual(values = my_palette) +
  theme_bw(base_size = 15) +
  labs(title = "Expression of C4A/C4B across cell types (Prenatal Brain)",
       x = NULL, 
       y = "Proportion-Corrected expression") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black", size = 15),
    axis.text.y = element_text(color = "black"),
    legend.position = "none", 
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold", size = 14)
  )
print(p)
setwd("/gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/figures/fig1")
ggsave("C4A_C4B_CellType_Expression_Corrected_prenatal.png", width = 7.2, height = 4)


########### 
library(tidyverse)
library(data.table)
library(ggpubr) 

# First, remove covariates (technical noise such as PC and PEER factors) at the gene expression level.
# Then, remove proportional interference at the cell abundance level (restoring the true expression level of single cells).


setwd("/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/removeversionIDmapping/cts_expr") 
cov_dir <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/cts-qtl-output-rmgeneversionID/stage1/"
files <- list.files(pattern = "^stage1_.*\\.txt$")
target_genes <- c("C4A", "C4B")

my_palette <- c(
  "prenatal_AST" = "#889ABF", "postnatal_AST" = "#889ABF", "AST" = "#889ABF",
  "prenatal_ExNeu" = "#A9D69C", "postnatal_ExNeu" = "#A9D69C", "ExNeu" = "#A9D69C",
  "prenatal_IN" = "#BB9BCA", "postnatal_IN" = "#BB9BCA", "IN" = "#BB9BCA",
  "prenatal_MG" = "#E59194", "postnatal_MG" = "#E59194", "MG" = "#E59194",
  "prenatal_OPC" = "#9BCCEA", "postnatal_OPC" = "#9BCCEA", "OPC" = "#9BCCEA",
  "prenatal_GLIALPROG" = "#CB6463", "postnatal_OL" = "#CB6463", "OL" = "#CB6463",
  "GLIALPROG" = "#CB6463" 
)

prop_file <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/pre_devbulk_prop.csv" 

prop_raw <- fread(prop_file, header = TRUE)
colnames(prop_raw)[1] <- "Sample"

prop_long <- prop_raw %>%
  pivot_longer(cols = -Sample, names_to = "CellType", values_to = "Proportion") %>%
  mutate(CellType = as.character(CellType)) 

data_list <- list()

for (f in files) {
  cell_type <- gsub("stage1_|\\.txt", "", f)
  cat("Processing:", cell_type, "\n")
  dt <- fread(f, header = TRUE)
  colnames(dt)[1] <- "Gene"
  dt_subset <- dt %>% filter(Gene %in% target_genes)
  
  if (nrow(dt_subset) > 0) {
    dt_long <- dt_subset %>%
      pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Expression") %>%
      mutate(CellType = cell_type)
    cov_file <- paste0(cov_dir, cell_type, "/covariates.cov.txt")
    
    if (file.exists(cov_file)) {
      cov_raw <- fread(cov_file, header = TRUE)
      colnames(cov_raw)[1] <- "Covariate"
      cov_df <- cov_raw %>%
        column_to_rownames("Covariate") %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column("Sample") %>%
        mutate(across(-Sample, as.numeric))
      
      dt_merged <- inner_join(dt_long, cov_df, by = "Sample")
      
      if (nrow(dt_merged) > 0) {
        cov_names <- setdiff(colnames(cov_df), "Sample")
        formula_str <- paste("Expression ~", paste(cov_names, collapse = " + "))
        
        corrected_list <- list()
        for (g in unique(dt_merged$Gene)) {
          tmp <- dt_merged[dt_merged$Gene == g, ]
          fit <- lm(as.formula(formula_str), data = tmp)
          tmp$Cov_Corrected_Expr <- residuals(fit) + mean(tmp$Expression, na.rm = TRUE)
          corrected_list[[g]] <- tmp
        }
      

        dt_res <- bind_rows(corrected_list) %>%
        dplyr::select(Gene, Sample, Expression, Cov_Corrected_Expr, CellType)
        
        data_list[[cell_type]] <- dt_res
      }
    } else {
     cat(" -> Warning: Covariate file not found", cov_file, "Using raw expression levels\n")
# If some cell types do not have covariates, use the raw expression levels directly to avoid errors.
      dt_long$Cov_Corrected_Expr <- dt_long$Expression
      data_list[[cell_type]] <- dt_long
    }
  }
}

final_df <- bind_rows(data_list)

plot_df <- final_df %>% 
  filter(!grepl("VASC", CellType)) %>%
  # Compare cell ratios
  inner_join(prop_long, by = c("Sample", "CellType")) %>%
  filter(Proportion > 0.05) %>%
  mutate(
   # Here, the starting input value has been changed to Cov_Corrected_Expr with covariates removed.
    Linear_Expr = 2^Cov_Corrected_Expr - 1,
    Corrected_Linear = Linear_Expr / Proportion,
    Corrected_Linear = pmax(Corrected_Linear, 0),
    Corrected_LogExpr = log2(Corrected_Linear + 1),
    CellType = gsub("prenatal_|postnatal_", "", CellType)
  )

print("Check formatted Cell Types for plotting:")
print(unique(plot_df$CellType))


p <- ggplot(plot_df, aes(x = reorder(CellType, -Corrected_LogExpr, median), y = Corrected_LogExpr, fill = CellType)) +
  geom_boxplot(alpha = 0.8, color = "black") + 
  facet_wrap(~Gene, scales = "free_y") + 
  scale_fill_manual(values = my_palette) +
  theme_bw(base_size = 16) +
  labs(title = "Expression of C4A/C4B across cell types (Prenatal Brain)",
       x = NULL, 
       y = "Corrected Expression (Log Scale)") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black", size = 15),
    axis.text.y = element_text(color = "black"),
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 15)
  )
print(p)

setwd("/gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/figures/fig1")
ggsave("C4A_C4B_CellType_Expression_FullycovCorrected2.png", width = 7.3, height = 4)

