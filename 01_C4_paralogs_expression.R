# C4 Lifespan Trajectory Analysis
#!/usr/bin/env Rscript
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

# 2. Target Genes Preparation (C4A / C4B)
target_genes <- c("C4A" = "ENSG00000244731.9", 
                  "C4B" = "ENSG00000224389.9")

valid_genes <- target_genes[target_genes %in% rownames(tpm)]
if(length(valid_genes) == 0) stop("Error: Target gene IDs not found in TPM matrix.")

# Extract target gene expression and merge with metadata
tpm_t <- as.data.frame(t(tpm[valid_genes, , drop = FALSE]))
colnames(tpm_t) <- names(valid_genes)
tpm_t$sample <- rownames(tpm_t)
df_target <- merge(meta, tpm_t, by = "sample")

# 3. Background Genes Preparation (Random sampling)

set.seed(123) # Ensure reproducibility
all_genes <- rownames(tpm)
bg_gene_ids <- sample(setdiff(all_genes, valid_genes), 1000) 

# Extract background gene expression and merge with metadata
bg_tpm_t <- as.data.frame(t(tpm[bg_gene_ids, , drop = FALSE]))
bg_tpm_t$sample <- rownames(bg_tpm_t)
df_bg <- merge(meta, bg_tpm_t, by = "sample")

# 4. Data Transformation and Z-score Normalization

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


# 5. Visualization: Trajectory Plot

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

# Save plot
ggsave("C4_expression_trajectory.png", plot = p_zscore, width =6, height =4.3, dpi = 500)
cat("Plot successfully saved to C4_expression_trajectory.png\n")


# 6. Statistical Analysis
# Prepare data for statistics
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

cat("\n================ STATISTICAL RESULTS ================\n")

# 6.1 Prenatal Stage (Stability testing)
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

# 6.2 Postnatal Stage (Increase testing)
cat("\n[2] Postnatal Stage (Pearson Correlation - Age vs TPM)\n")
for (g in c("C4A", "C4B")) {
  tmp <- stat_data %>% filter(Gene == g, Phase == "Postnatal")
  res <- cor.test(tmp$age, tmp$TPM, method = "pearson")
  cat(sprintf("%s : r = %.3f, P = %.2e\n", g, res$estimate, res$p.value))
}

# 6.3 Developmental Stages Comparison (ANOVA)
cat("\n[3] Stage Comparison (ANOVA: Prenatal vs 0-65y vs >65y)\n")
for (g in c("C4A", "C4B")) {
  tmp <- stat_data %>% filter(Gene == g)
  res <- aov(TPM ~ Age_Group, data = tmp)
  f_val <- summary(res)[[1]][["F value"]][1]
  p_val <- summary(res)[[1]][["Pr(>F)"]][1]
  cat(sprintf("%s : F = %.2f, P = %.2e\n", g, f_val, p_val))
}

#[3] Stage Comparison (ANOVA: Prenatal vs 0-65y vs >65y)
#C4A : F = 13.50, P = 3.87e-06
#C4B : F = 9.13, P = 1.78e-04