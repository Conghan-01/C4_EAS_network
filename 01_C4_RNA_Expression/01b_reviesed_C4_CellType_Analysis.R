# C4_CellType_Analysis.R: Paralog Cell-Type Specific Expression 
# This script performs deconvolution-based cell-type correlation, 
# cell-type-specific expression analysis (bMIND-inferred), and statistical comparisons for C4A and C4B in human brain.
# sex-bias analysis for C4A and C4B in human brain

### Sections
### Search by
# 1a. Correlation of C4 exp and prenatal each cell-type correlation: scatter plot
# 1b. Correlation of C4 exp and postnatal each cell-type correlation: scatter plot
# 2a. Barplot: Correlation of C4 gene expression and total cell proportions of prenatal and postnatal
# 2b. Barplot: Correlation of C4 gene expression and core cell (4 types) proportions of prenatal and postnatal
# 3. Cell-type-specific expression analysis (direct bMIND inferred expression without double proportion correction)

# ==============================================================================
# 1a. Correlation of C4 exp and prenatal each cell-type correlation: scatter plot
# ==============================================================================

library(data.table)
library(tidyverse)
library(ggpubr)

setwd("/gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/figures/fig1")
prop_file <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/pre_devbulk_prop.csv" 
expr_file <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/bulk_data/stage1.logTMM.ComBat.txt"

target_genes <- c("C4A" = "ENSG00000244731.9", 
                  "C4B" = "ENSG00000224389.9")
props <- read.csv(prop_file, row.names = 1, check.names = FALSE)
props$sample <- rownames(props)
print("Cell Types found (Prenatal):")
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

target_cells <- c("GLIALPROG","prenatal_MG", "prenatal_AST", "prenatal_ExNeu", "prenatal_IN", "prenatal_OPC")
plot_subset <- df_long %>% filter(CellType %in% target_cells)
my_palette <- c("#CB6463", "#889ABF", "#A9D69C", "#BB9BCA", "#E59194", "#9BCCEA")
df_long_prenatal <- plot_subset

p_scatter <- ggplot(plot_subset, aes(x = Proportion, y = Expression)) +
  geom_point(aes(color = CellType), alpha = 0.6, size = 1.5) +
  geom_smooth(method = "lm", color = "grey30", size = 0.5, se = TRUE) + 
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


# ==============================================================================
# 1b. Postnatal cell-type correlation: scatter plot
# ==============================================================================

prop_file <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/post_devbulk_prop.csv" 
expr_file <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/bulk_data/postnatal.logTMM.ComBat.txt"
target_genes <- c("C4A" = "ENSG00000244731.9", 
                  "C4B" = "ENSG00000224389.9")

props <- read.csv(prop_file, row.names = 1, check.names = FALSE)
props$sample <- rownames(props)
print("Cell Types found (Postnatal):")
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
df_long_postnatal <- plot_subset
my_palette <- c("#CB6463", "#889ABF", "#A9D69C", "#BB9BCA", "#E59194", "#9BCCEA")

p_scatter <- ggplot(plot_subset, aes(x = Proportion, y = Expression)) +
  geom_point(aes(color = CellType), alpha = 0.6, size = 1.5) +
  geom_smooth(method = "lm", color = "black", size = 0.5, se = TRUE) + 
  facet_grid(Gene ~ CellType, scales = "free") +
  stat_cor(method = "pearson", 
           label.x.npc = "left", 
           label.y.npc = "top", 
           size = 4, 
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


# ==============================================================================
# 2a. Barplot: Correlation of C4 gene expression and total cell proportions of prenatal and postnatal
# ==============================================================================

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

cell_order <- plot_bar_sorted %>%
  group_by(Cell_Label) %>%
  summarise(Mean_R = mean(R, na.rm = TRUE)) %>%
  arrange(desc(Mean_R)) %>%
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
    plot.title = element_text(hjust = 0.5, size = 18),
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


# Generate Supplementary Table 3: Correlation & Cell Proportion Stats
cat("\nGenerating Supplementary Table 3...\n")
get_comprehensive_stats <- function(df_long, stage_label) {
  df_long %>%
    group_by(Gene, CellType) %>%
    summarise(
      Median_Prop = median(Proportion, na.rm = TRUE),
      Q1_Prop = quantile(Proportion, 0.25, na.rm = TRUE),
      Q3_Prop = quantile(Proportion, 0.75, na.rm = TRUE),
      Pearson_r = cor(Proportion, Expression, method = "pearson", use = "complete.obs"),
      P_value = cor.test(Proportion, Expression)$p.value,
      .groups = "drop"
    ) %>%
    mutate(Developmental_Stage = stage_label)
}

table_pre <- get_comprehensive_stats(df_long_prenatal, "Prenatal")
table_post <- get_comprehensive_stats(df_long_postnatal, "Postnatal")
supp_table3 <- bind_rows(table_pre, table_post)
supp_table3 <- supp_table3 %>%
  group_by(Developmental_Stage) %>%
  mutate(FDR = p.adjust(P_value, method = "BH")) %>%
  ungroup()

supp_table3_final <- supp_table3 %>%
  mutate(
    Cell_Type = gsub("prenatal_|postnatal_", "", CellType),
    `Median_Proportion_[IQR]_(%)` = sprintf("%.1f%% [%.1f%% - %.1f%%]", 
                                            Median_Prop * 100, 
                                            Q1_Prop * 100, 
                                            Q3_Prop * 100),
    Pearson_r = round(Pearson_r, 3),
    P_value = signif(P_value, 3),
    FDR = signif(FDR, 3)
  ) %>%
  dplyr::select(Developmental_Stage, Gene, Cell_Type, `Median_Proportion_[IQR]_(%)`, Pearson_r, P_value, FDR) %>%
  arrange(Developmental_Stage, Gene, desc(Pearson_r))

print(head(supp_table3_final))
output_csv <- "Supplementary_Table_3_CellType_Correlation.csv"
write.csv(supp_table3_final, output_csv, row.names = FALSE, quote = FALSE)


# ==============================================================================
# 2b. Barplot: Correlation of C4 gene expression and core cell (4 types) proportions
# ==============================================================================

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
  geom_col(position = position_dodge(width = 0.8), width = 0.75, color = "white", linewidth = 0.2) +
  geom_text(aes(y = Label_Y, label = Stars), 
            position = position_dodge(width = 0.8), 
            vjust = 0.5, 
            size = 4, fontface = "bold", color = "black", show.legend = FALSE) +
  scale_fill_manual(
    values = c("C4A" = "#d25756", "C4B" = "#7eb4db"),
    name = "Gene"
  ) +
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


# ==============================================================================
# 3. Cell-Type-Specific Expression Estimation (Postnatal & Prenatal)
# Direct use of bMIND-inferred expression (no double proportion correction)
# ==============================================================================

theme_set(theme_bw(base_size = 17))

# ------------------------------------------------------------------------------
# 3a. Postnatal Brain: Cell-Type Expression Processing
# ------------------------------------------------------------------------------

setwd("/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/removeversionIDmapping/cts_expr") 
files <- list.files(pattern = "^adult_.*\\.txt$")
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

prop_file <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/post_devbulk_prop.csv" 
prop_raw <- fread(prop_file, header = TRUE)
colnames(prop_raw)[1] <- "Sample"
prop_long <- prop_raw %>%
  pivot_longer(cols = -Sample, names_to = "CellType", values_to = "Proportion") %>%
  mutate(CellType = as.character(CellType)) 

data_list <- list()

for (f in files) {
  cell_type <- gsub("adult_|\\.txt", "", f)
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

# Data merging and filtering (Directly using bMIND expression, skipping double proportion division)
plot_df <- final_df %>% 
  filter(!grepl("VASC|OL|OPC", CellType)) %>%
  inner_join(prop_long, by = c("Sample", "CellType")) %>%
  filter(Proportion > 0.05) %>%
  mutate(
    # bMIND inferred expression is used directly without dividing by Proportion
    Corrected_LogExpr = Expression,
    CellType = gsub("prenatal_|postnatal_", "", CellType)
  )

plot_df <- plot_df %>%
  mutate(CellType = fct_reorder(CellType, Corrected_LogExpr, .fun = median, .desc = TRUE))

print("Check formatted and ORDERED Cell Types for plotting (Postnatal):")
print(levels(plot_df$CellType)) 

p <- ggplot(plot_df, aes(x = CellType, y = Corrected_LogExpr, fill = CellType)) +
  geom_boxplot(alpha = 0.8, color = "black") + 
  facet_wrap(~Gene, scales = "free_y") + 
  scale_fill_manual(values = my_palette) +
  theme_bw(base_size = 15) +
  labs(title = "Expression of C4A/C4B across cell types (Postnatal Brain)",
       x = NULL, 
       y = "bMIND Inferred Expression") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black", size = 15),
    axis.text.y = element_text(color = "black"),
    legend.position = "none", 
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold", size = 14)
  )
print(p)

setwd("/gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/figures/fig1")
ggsave("C4A_C4B_CellType_Expression_Corrected_postnatal_reviesed.png", width = 7.2, height = 4)
plot_df_postnatal <- plot_df %>% mutate(Stage = "Postnatal")


# ------------------------------------------------------------------------------
# 3b. Prenatal Brain: Cell-Type Expression Processing (with Covariate Adjustment)
# ------------------------------------------------------------------------------

setwd("/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/removeversionIDmapping/cts_expr") 
cov_dir <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/cts-qtl-output-rmgeneversionID/stage1/"
files <- list.files(pattern = "^stage1_.*\\.txt$")
target_genes <- c("C4A", "C4B")

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
      dt_long$Cov_Corrected_Expr <- dt_long$Expression
      data_list[[cell_type]] <- dt_long
    }
  }
}

final_df <- bind_rows(data_list)

plot_df <- final_df %>% 
  filter(!grepl("VASC", CellType)) %>%
  inner_join(prop_long, by = c("Sample", "CellType")) %>%
  filter(Proportion > 0.05) %>%
  mutate(
    # Covariate-adjusted bMIND expression is used directly without dividing by Proportion
    Corrected_LogExpr = Cov_Corrected_Expr,
    CellType = gsub("prenatal_|postnatal_", "", CellType)
  )

print("Check formatted Cell Types for plotting (Prenatal):")
print(unique(plot_df$CellType))

p <- ggplot(plot_df, aes(x = reorder(CellType, -Corrected_LogExpr, median), y = Corrected_LogExpr, fill = CellType)) +
  geom_boxplot(alpha = 0.8, color = "black") + 
  facet_wrap(~Gene, scales = "free_y") + 
  scale_fill_manual(values = my_palette) +
  theme_bw(base_size = 16) +
  labs(title = "Expression of C4A/C4B across cell types (Prenatal Brain)",
       x = NULL, 
       y = "bMIND Inferred Expression") +
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
ggsave("C4A_C4B_CellType_Expression_Corrected_prenatal_reviesed.png", width = 7.3, height = 4)
plot_df_prenatal <- plot_df %>% mutate(Stage = "Prenatal")


# ------------------------------------------------------------------------------
# 3c. Generate Supplementary Table 4: AST Comparison Statistical Tests
# ------------------------------------------------------------------------------

cat("\nGenerating Supplementary Table 4 (AST Comparison Stats)...\n")
all_corrected_df <- bind_rows(plot_df_prenatal, plot_df_postnatal)

generate_ast_stats <- function(data) {
  results <- list()
  for(stg in unique(data$Stage)) {
    for(g in unique(data$Gene)) {
      tmp <- data %>% filter(Stage == stg, Gene == g)
      if(!("AST" %in% tmp$CellType)) next 
      
      ast_data <- tmp %>% filter(CellType == "AST") %>% pull(Corrected_LogExpr)
      ast_median <- median(ast_data, na.rm = TRUE)
      other_cells <- setdiff(unique(tmp$CellType), "AST")
      for(cell in other_cells) {
        comp_data <- tmp %>% filter(CellType == cell) %>% pull(Corrected_LogExpr)
        comp_median <- median(comp_data, na.rm = TRUE)
        
        test_res <- wilcox.test(ast_data, comp_data, alternative = "two.sided")
        
        results[[length(results) + 1]] <- data.frame(
          Developmental_Stage = stg,
          Gene = g,
          Reference_Cell = "AST",
          Comparison_Cell = cell,
          Median_Expr_AST = round(ast_median, 3),
          Median_Expr_Comparison = round(comp_median, 3),
          Log2_Fold_Difference = round(ast_median - comp_median, 3), 
          P_value = test_res$p.value
        )
      }
    }
  }
  res_df <- bind_rows(results)
  res_df <- res_df %>%
    group_by(Developmental_Stage, Gene) %>%
    mutate(FDR = p.adjust(P_value, method = "BH")) %>% 
    ungroup() %>%
    mutate(
      P_value = signif(P_value, 3),
      FDR = signif(FDR, 3)
    ) %>%
    arrange(Developmental_Stage, Gene, desc(Log2_Fold_Difference))
  
  return(res_df)
}

supp_table4_final <- generate_ast_stats(all_corrected_df)
print(head(supp_table4_final))
output_csv4 <- "Supplementary_Table_4_AST_Predominant_Expression_Stats_revised.csv"
setwd("/gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/figures/fig1")
write.csv(supp_table4_final, output_csv4, row.names = FALSE, quote = FALSE)

