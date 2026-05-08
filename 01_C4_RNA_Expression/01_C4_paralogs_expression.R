# This script contains two parts: 1. The changes in C4 expression in Lifespan
#                                 2. The sex differences of C4 in brain development

######### 1. C4 Lifespan Trajectory Analysis

# ==============================================================================
# Description: Analyzes and visualizes the lifespan expression trajectories of 
#              C4 paralogs (C4A and C4B) compared to background genes using 
#              human prefrontal cortex RNA-seq data.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
})

setwd("/gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/figures/fig1")
meta_path <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/01mixbatches_group/meta.txt2"
tpm_path  <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/01mixbatches_group/subset_tpm.txt"

# Read metadata and TPM matrix
meta <- read.table(meta_path, header = TRUE, stringsAsFactors = FALSE)
tpm <- read.table(tpm_path, header = TRUE, row.names = 1, check.names = FALSE)

# Target Genes Preparation (C4A / C4B)
target_genes <- c("C4A" = "ENSG00000244731.9", 
                  "C4B" = "ENSG00000224389.9")

valid_genes <- target_genes[target_genes %in% rownames(tpm)]
if(length(valid_genes) == 0) stop("Error: Target gene IDs not found in TPM matrix.")

# Extract target gene expression and merge with metadata
tpm_t <- as.data.frame(t(tpm[valid_genes, , drop = FALSE]))
colnames(tpm_t) <- names(valid_genes)
tpm_t$sample <- rownames(tpm_t)
df_target <- merge(meta, tpm_t, by = "sample")

# Background Genes Preparation (Random sampling)

set.seed(123) # Ensure reproducibility
all_genes <- rownames(tpm)
bg_gene_ids <- sample(setdiff(all_genes, valid_genes), 1000) 

# Extract background gene expression and merge with metadata
bg_tpm_t <- as.data.frame(t(tpm[bg_gene_ids, , drop = FALSE]))
bg_tpm_t$sample <- rownames(bg_tpm_t)
df_bg <- merge(meta, bg_tpm_t, by = "sample")

# Data Transformation and Z-score Normalization

prep_time_cols <- function(df_input) {
  df_input %>%
    mutate(
      PCD = case_when(
        Stage == "stage1" ~ age * 7, # Gestational weeks to days
        Stage %in% c("stage2", "stage3") ~ age * 365 + 280 # Years to post-conception days
      ),
      Log_PCD = log10(PCD),
      Phase = case_when(
        Stage == "stage1" ~ "Prenatal",
        Stage %in% c("stage2", "stage3") ~ "Postnatal"
      ),
      Age_Visual = age
    ) %>%
    mutate(Phase = factor(Phase, levels = c("Prenatal", "Postnatal")))
}


long_target <- df_target %>% 
  pivot_longer(cols = all_of(names(valid_genes)), names_to = "Gene", values_to = "TPM") %>%
  prep_time_cols() %>%
  mutate(Type = "Target") 

long_bg <- df_bg %>%
  pivot_longer(cols = all_of(bg_gene_ids), names_to = "Gene", values_to = "TPM") %>%
  prep_time_cols() %>%
  mutate(Type = "Background") 

# Combine datasets and calculate Z-scores per gene
combined_data_z <- bind_rows(long_target, long_bg) %>%
  group_by(Gene) %>%
  mutate(Z_Score = scale(TPM)[,1]) %>%
  ungroup()


# Visualization: Trajectory Plot

p_zscore <- ggplot() +
  # Layer 1: Background genes (Grey)
  geom_smooth(data = combined_data_z %>% filter(Type == "Background"), 
              aes(x = Age_Visual, y = Z_Score, group = Gene), 
              method = "loess", se = FALSE, color = "grey90", alpha = 0.4, linewidth = 0.5) +
  
  # Layer 2: Target genes (C4A and C4B)
  geom_smooth(data = combined_data_z %>% filter(Type == "Target"), 
              aes(x = Age_Visual, y = Z_Score, color = Gene, fill = Gene),
              method = "loess", alpha = 0.2, linewidth = 1.5, span = 0.8) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.6) +
  facet_wrap(~Phase, scales = "free_x", strip.position = "top") +
  scale_color_manual(values = c("C4A" = "#E69191", "C4B" = "#92B5CA")) +
  scale_fill_manual(values = c("C4A" = "#E69191", "C4B" = "#92B5CA")) +
  coord_cartesian(ylim = c(-2.5, 4)) +
  theme_bw() + 
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90", color = "black", size = 1),
    strip.text = element_text(size = 16),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16),
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(size = 20, hjust = 0.5)
  ) +
  labs(
    title = "Lifespan Expression Trajectory",
    x = "Age (Left: Gestational Weeks / Right: Postnatal Years)", 
    y = "Standardized Expression (Z-score)"
  )


ggsave("C4_expression_trajectory.png", plot = p_zscore, width =6, height =4.3, dpi = 500)
cat("Plot successfully saved to C4_expression_trajectory.png\n")


# Statistical Analysis
stat_data <- df_target %>%
  pivot_longer(cols = c("C4A", "C4B"), names_to = "Gene", values_to = "TPM") %>%
  mutate(
    Phase = case_when(
      Stage == "stage1" ~ "Prenatal",
      Stage %in% c("stage2", "stage3") ~ "Postnatal"
    ),
    Age_Group = case_when(
      Phase == "Prenatal" ~ "Prenatal",
      Phase == "Postnatal" & age <= 65 ~ "0-65y",
      Phase == "Postnatal" & age > 65 ~ ">65y"
    )
  ) %>%
  mutate(Age_Group = factor(Age_Group, levels = c("Prenatal", "0-65y", ">65y")))


# Prenatal Stage (Stability testing)
cat("\n[1] Prenatal Stage (Pearson Correlation - Age vs TPM)\n")
for (g in c("C4A", "C4B")) {
  tmp <- stat_data %>% filter(Gene == g, Phase == "Prenatal")
  res <- cor.test(tmp$age, tmp$TPM, method = "pearson")
  cat(sprintf("%s : r = %.3f, P = %.2e\n", g, res$estimate, res$p.value))
}

#[1] Prenatal Stage (Pearson Correlation - Age vs TPM)
#C4A : r = 0.140, P = 2.14e-01
#C4B : r = 0.107, P = 3.43e-01
#[2] Postnatal Stage (Pearson Correlation - Age vs TPM)
#C4A : r = 0.225, P = 4.50e-02
#C4B : r = 0.316, P = 4.35e-03

# Postnatal Stage (Increase testing)
cat("\n[2] Postnatal Stage (Pearson Correlation - Age vs TPM)\n")
for (g in c("C4A", "C4B")) {
  tmp <- stat_data %>% filter(Gene == g, Phase == "Postnatal")
  res <- cor.test(tmp$age, tmp$TPM, method = "pearson")
  cat(sprintf("%s : r = %.3f, P = %.2e\n", g, res$estimate, res$p.value))
}

# Developmental Stages Comparison (ANOVA)
cat("\n[3] Stage Comparison (ANOVA: Prenatal vs 0-65y vs >65y)\n")
for (g in c("C4A", "C4B")) {
  tmp <- stat_data %>% filter(Gene == g)
  res <- aov(TPM ~ Age_Group, data = tmp)
  f_val <- summary(res)[[1]][["F value"]][1]
  p_val <- summary(res)[[1]][["Pr(>F)"]][1]
  cat(sprintf("%s : F = %.2f, P = %.2e\n", g, f_val, p_val))
}

#[3] Stage Comparison (Prenatal vs 0-65y vs >65y)
#C4A : F = 13.50, P = 3.87e-06
#C4B : F = 9.13, P = 1.78e-04


# supplemental fig1 expression
plot_data <- df_target %>%
  pivot_longer(cols = c("C4A", "C4B"), names_to = "Gene", values_to = "TPM") %>%
  mutate(
    Phase = case_when(
      Stage == "stage1" ~ "Prenatal",
      Stage %in% c("stage2", "stage3") ~ "Postnatal"
    ),

    Age_Group = case_when(
      Phase == "Prenatal" ~ "Prenatal",
      Phase == "Postnatal" & age <= 65 ~ "0-65y",
      Phase == "Postnatal" & age > 65 ~ ">65y"
    ),
    Log2TPM = log2(TPM + 1) 
  ) %>%
  mutate(
    Age_Group = factor(Age_Group, levels = c("Prenatal", "0-65y", ">65y")),
    Gene = factor(Gene, levels = c("C4A", "C4B"))
  )

my_comparisons <- list(
  c("Prenatal", "0-65y"), 
  c("0-65y", ">65y"), 
  c("Prenatal", ">65y")
)

p_violin <- ggplot(plot_data, aes(x = Age_Group, y = Log2TPM, fill = Age_Group)) +
  geom_violin(trim = FALSE, alpha = 0.5, color = NA) +
  geom_boxplot(width = 0.15, fill = "white", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 1.2, color = "black", alpha = 0.4) +
  facet_wrap(~ Gene) +
  
  scale_fill_manual(values = c("Prenatal" = "#AECBEB", 
                               "0-65y"    = "#C1DDB4", 
                               ">65y"     = "#F4ECA1")) +
  
stat_compare_means(
    comparisons = my_comparisons, 
    method = "wilcox.test",     
    label = "p.signif",    
    step.increase = 0.12,   
    tip.length = 0.02,
    hide.ns = FALSE         
  ) +
  
  labs(
    title = "C4A & C4B Expression Across Stages",
    x = NULL, 
    y = "Expression (log2(TPM + 1))"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5, color = "black"),
    axis.text.x = element_text(size = 15, color = "black"),
    axis.text.y = element_text(size = 15, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    strip.text = element_text(size = 16, color = "black"),
    strip.background = element_rect(fill = "white", color = "black"), 
    legend.position = "none", 
    panel.grid.major.x = element_blank(), 
    panel.border = element_rect(color = "black", linewidth = 1)
  )

print(p_violin)
ggsave("C4_expression_violin_with_p.png", plot = p_violin, width =6, height =4, dpi = 500)

#p value for c("Prenatal", "0-65y"): C4A 0.0056 ; C4B 0.0078
#p value for  c("0-65y", ">65y"): C4A  0.0042 ; C4B 0.00073
#p value for c("Prenatal", ">65y"): C4A  1.8e-6 ; C4B:3.9e-7



######### 2. Sex difference of C4 expression ########
# ==============================================================================
# C4 Sex Difference Analysis (Corrected for Age and RIN)
# ==============================================================================
prep_time_cols <- function(df_input) {
  df_input %>%
    mutate(
      PCD = case_when(
        Stage == "stage1" ~ age * 7, # Gestational weeks to days
        Stage %in% c("stage2", "stage3") ~ age * 365 + 280 # Years to post-conception days
      ),
      Log_PCD = log10(PCD),
      Phase = case_when(
        Stage == "stage1" ~ "Prenatal",
        Stage %in% c("stage2", "stage3") ~ "Postnatal"
      ),
      Age_Visual = age
    ) %>%
    mutate(Phase = factor(Phase, levels = c("Prenatal", "Postnatal")))
}

df_target$sex<- as.factor(df_target$sex)
df_sex_analysis <- df_target %>%
  prep_time_cols() %>%
  mutate(across(c(age, RIN, sex), as.numeric))
get_corrected_sex_expr <- function(data_subset) {
  res_list <- list()
  for (g in c("C4A", "C4B")) {
    tmp <- data_subset %>% 
      select(sample, Gene = all_of(g), age, RIN, sex, Phase) %>%
      filter(!is.na(Gene))
    fit <- lm(Gene ~ age + RIN, data = tmp)
    tmp$Corrected_Value <- residuals(fit) + mean(tmp$Gene, na.rm = TRUE)
    tmp$Target_Gene <- g
    res_list[[g]] <- tmp
  }
  return(bind_rows(res_list))
}


sex_plot_df <- df_sex_analysis %>%
  group_split(Phase) %>%
  map_dfr(~get_corrected_sex_expr(.x)) %>%
  mutate(
    Sex_Label = ifelse(sex == 2, "Female", "Male"),
    Sex_Label = factor(Sex_Label, levels = c("Female", "Male")),
    Phase = factor(Phase, levels = c("Prenatal", "Postnatal")),
    Target_Gene = factor(Target_Gene, levels = c("C4A", "C4B"))
  )

sex_colors <- c("Female" = "#CFA9C5", "Male" = "#96D1CD")

p_sex <- ggplot(sex_plot_df, aes(x = Phase, y = Corrected_Value, fill = Sex_Label)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.8, color = "black", 
               position = position_dodge(width = 0.7)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7),
             size = 1, alpha = 0.3, color = "black") +
  facet_wrap(~Target_Gene, scales = "free_y") +
  stat_compare_means(aes(group = Sex_Label), 
                     method = "wilcox.test", 
                     label = "p.signif", 
                     label.y.npc = "top",
                     vjust = 0.5, size = 5) +
  scale_fill_manual(values = sex_colors) +
  theme_bw(base_size = 17) +
  labs(
    title = "Sex Differences in C4 Expression",
    x = "Developmental Stages",
    y = "Expression (Covariate-Corrected TPM)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey30"),
    legend.position = "top",
    legend.title = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )
print(p_sex)
ggsave("C4_Sex_Difference_Age_RIN_Corrected.png", p_sex, width = 6.7, height = 5, dpi = 500)


# ==============================================================================
# C4 Sex Difference Analysis with Hidden Covariate Correction
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(data.table)
})


prep_data <- function(df) {
  df %>% mutate(
    Phase = ifelse(Stage == "stage1", "Prenatal", "Postnatal"),
    Phase = factor(Phase, levels = c("Prenatal", "Postnatal")),
    sex = as.numeric(sex)
  )
}

get_hidden_covs <- function(expr_mat, n_pcs = 5) {
  dat <- t(expr_mat[rowMeans(expr_mat) > 1, ])
  pca <- prcomp(dat, scale. = TRUE)
  pcs <- as.data.frame(pca$x[, 1:n_pcs])
  colnames(pcs) <- paste0("Hidden_PC", 1:n_pcs)
  pcs$sample <- rownames(pcs)
  return(pcs)
}

run_sex_analysis <- function(target_df, expr_mat) {
  hidden_pcs <- get_hidden_covs(expr_mat, n_pcs = 3)
  combined <- target_df %>% prep_data() %>% inner_join(hidden_pcs, by = "sample")
  
  results_list <- list()
  
  for (ph in c("Prenatal", "Postnatal")) {
    subset_df <- combined %>% filter(Phase == ph)
    clean_pcs <- c()
    for (i in 1:5) {
      pc_col <- paste0("Hidden_PC", i)
      p_sex <- t.test(subset_df[[pc_col]] ~ subset_df$sex)$p.value
      if (p_sex > 0.05) { clean_pcs <- c(clean_pcs, pc_col) }
    }
    covs <- c("age", "RIN", clean_pcs)
    formula_str <- paste("Value ~", paste(covs, collapse = " + "))
    
    for (gene in c("C4A", "C4B")) {
      tmp <- subset_df %>% select(sample, Value = all_of(gene), sex, Phase, all_of(covs)) %>% filter(!is.na(Value))
      
      fit <- lm(as.formula(formula_str), data = tmp)
      tmp$Corrected_Expr <- residuals(fit) + mean(tmp$Value, na.rm = TRUE)
      tmp$Gene_Name <- gene
      results_list[[paste(ph, gene)]] <- tmp
    }
  }
  return(bind_rows(results_list))
}

final_plot_df <- run_sex_analysis(df_target, tpm)

final_plot_df <- final_plot_df %>%
  mutate(Sex_Label = ifelse(sex == 2, "Female", "Male"),
         Sex_Label = factor(Sex_Label, levels = c("Female", "Male")))

sex_colors <- c("Female" = "#CFA9C5", "Male" = "#96D1CD")

p_sex_final <- ggplot(final_plot_df, aes(x = Phase, y = Corrected_Expr, fill = Sex_Label)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.8, color = "black", 
               position = position_dodge(width = 0.7)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7),
             size = 0.8, alpha = 0.3, color = "black") +
  facet_wrap(~Gene_Name, scales = "free_y") +
  stat_compare_means(aes(group = Sex_Label), 
                     method = "t.test", 
                     label = "p.signif", 
                     label.y.npc = "top", 
                     size = 5) +
  scale_fill_manual(values = sex_colors) +
  theme_bw(base_size = 16) +
  labs(
    title = "Sex Differences in C4 Expression (Hidden-PC Corrected)",
    subtitle = "Corrected for Age, RIN, and Sex-independent Hidden Factors",
    x = NULL, y = "Corrected Expression (TPM)"
  ) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p_sex_final)
ggsave("C4_Sex_Difference_Hidden_Corrected.png", width = 6.7, height = 5, dpi = 400)


# Before comparing genders, we should filter outliers from the samples and then correct for RIN and age.

