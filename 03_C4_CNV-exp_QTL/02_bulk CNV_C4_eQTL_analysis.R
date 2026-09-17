#############  Prenatal ######################################

library(tidyverse)
library(data.table)
library(patchwork)

# 1. Path Configuration

hap_file <- "/gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/topmed_analysis/final_check/fetal_onlychb_imputed_haps_refined.R5.txt"
expr_file <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/bulk_data/stage1.logTMM.ComBat.txt" 
cov_file <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/bulk_data/prenatalcovariatesToUse.PCAforQTLfromINT.txt"

# ==============================================================================
# 2. Parse Haplotype Data (Copy Number Calculation)
# ==============================================================================

haps <- fread(hap_file, header = TRUE)

get_cn <- function(hap_vec, pos) {
  as.numeric(sapply(strsplit(hap_vec, "_"), `[`, pos))
}

cn_df <- data.frame(
  Sample = haps$SAMPLE,
  Total_C4 = get_cn(haps$H1, 2) + get_cn(haps$H2, 2),
  C4A_CN   = get_cn(haps$H1, 3) + get_cn(haps$H2, 3),
  C4B_CN   = get_cn(haps$H1, 4) + get_cn(haps$H2, 4)
)

# ==============================================================================
# 3. Process Expression Data
# ==============================================================================

expr_data <- fread(expr_file, header = TRUE)
colnames(expr_data)[1] <- "Gene"
# Filter for C4A (including Ensembl ID) and transpose
c4a_expr <- expr_data %>%
  #filter(Gene == "C4B" | grepl("ENSG00000224389", Gene)) %>% 
  filter(Gene == "C4A" | grepl("ENSG00000244731", Gene)) %>%
  head(1) %>% 
  column_to_rownames("Gene") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Sample")

colnames(c4a_expr)[2] <- "Raw_Expr"

# ==============================================================================
# 4.Process Covariates
# ==============================================================================
cov_raw <- fread(cov_file, header = TRUE)
colnames(cov_raw)[1] <- "Covariate"

cov_df <- cov_raw %>%
  column_to_rownames("Covariate") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  mutate(across(-Sample, as.numeric)) 

# ==============================================================================
# 5. Data Merging & Covariate Correction
# ==============================================================================
df_merged <- inner_join(cn_df, c4a_expr, by = "Sample") %>%
  inner_join(cov_df, by = "Sample")


cov_names <- setdiff(colnames(cov_df), "Sample")
formula_str <- paste("Raw_Expr ~", paste(cov_names, collapse = " + "))

fit_cov <- lm(as.formula(formula_str), data = df_merged)
df_merged$Corrected_Expr <- residuals(fit_cov) + mean(df_merged$Raw_Expr, na.rm = TRUE)
df_prenatal <- df_merged
# ==============================================================================
# 6. Visualization
# ==============================================================================
plot_c4_eqtl <- function(data, x_col, title_text, y_label = "Normalized Expression (C4A)") {
  formula_eqtl <- as.formula(paste("Corrected_Expr ~", x_col))
  fit_eqtl <- lm(formula_eqtl, data = data)
  
  slope <- coef(fit_eqtl)[2]
  p_val <- summary(fit_eqtl)$coefficients[2, 4]
  
  title_color <- ifelse(p_val < 0.05, "#d7191c", "black")
  subtitle_text <- sprintf("Slope = %.3f, P = %.2e", slope, p_val)
  
  p <- ggplot(data, aes(x = .data[[x_col]], y = Corrected_Expr)) +
    geom_violin(aes(group = .data[[x_col]]), fill = "#a6cee3", color = NA, alpha = 0.6) +
    geom_boxplot(aes(group = .data[[x_col]]), width = 0.15, fill = "white", color = "black", outlier.shape = NA) +
    geom_jitter(width = 0.1, color = "grey40", size = 1.2, alpha = 0.8) +
    geom_smooth(method = "lm", color = "#b2182b", linetype = "dashed", se = TRUE, fill = "grey80", linewidth = 1) +
    scale_x_continuous(breaks = seq(min(data[[x_col]], na.rm=TRUE), max(data[[x_col]], na.rm=TRUE), by = 1)) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 13, color = title_color),
      axis.text = element_text(color = "black")
    ) +
    labs(title = title_text, subtitle = subtitle_text, x = x_col, y = y_label)  
  return(p)
}

# Combine plot
p1 <- plot_c4_eqtl(df_merged, "Total_C4", "Total C4 CN vs C4A Expr")
p2 <- plot_c4_eqtl(df_merged, "C4A_CN", "C4A CN vs C4A Expr")
p3 <- plot_c4_eqtl(df_merged, "C4B_CN", "C4B CN vs C4A Expr")

final_plot <- p1 | p2 | p3

final_plot <- final_plot + plot_annotation(
  title = "Prenatal Brain: C4A Expression vs Copy Numbers",
  theme = theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))
)

print(final_plot) 

ggsave("Fetal_Brain_C4A_eQTL_Bulk_Corrected.png", final_plot, width = 10, height = 4, dpi = 400)


#############  Postnatal ######################################

library(tidyverse)
library(data.table)
library(patchwork)


hap_file <- "/gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/topmed_analysis/final_check/adult/adult_only_imputed_haps_refinedR5revisedID.txt"
expr_file <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/bulk_data/postnatal.logTMM.ComBat.txt" 
cov_file <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/bulk_data/postnatalcovariatesToUse.PCAforQTLfromINT.txt"


haps <- fread(hap_file, header = TRUE)

get_cn <- function(hap_vec, pos) {
  as.numeric(sapply(strsplit(hap_vec, "_"), `[`, pos))
}

cn_df <- data.frame(
  Sample = haps$SAMPLE,
  Total_C4 = get_cn(haps$H1, 2) + get_cn(haps$H2, 2),
  C4A_CN   = get_cn(haps$H1, 3) + get_cn(haps$H2, 3),
  C4B_CN   = get_cn(haps$H1, 4) + get_cn(haps$H2, 4)
)

expr_data <- fread(expr_file, header = TRUE)
colnames(expr_data)[1] <- "Gene"

c4a_expr <- expr_data %>%
  #filter(Gene == "C4B" | grepl("ENSG00000224389", Gene)) %>% 
  filter(Gene == "C4A" | grepl("ENSG00000244731", Gene)) %>%
  head(1) %>% 
  column_to_rownames("Gene") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Sample")

colnames(c4a_expr)[2] <- "Raw_Expr"

cov_raw <- fread(cov_file, header = TRUE)
colnames(cov_raw)[1] <- "Covariate"

cov_df <- cov_raw %>%
  column_to_rownames("Covariate") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  mutate(across(-Sample, as.numeric)) 

df_merged <- inner_join(cn_df, c4a_expr, by = "Sample") %>%
  inner_join(cov_df, by = "Sample")

cov_names <- setdiff(colnames(cov_df), "Sample")
formula_str <- paste("Raw_Expr ~", paste(cov_names, collapse = " + "))

fit_cov <- lm(as.formula(formula_str), data = df_merged)
df_merged$Corrected_Expr <- residuals(fit_cov) + mean(df_merged$Raw_Expr, na.rm = TRUE)
df_postnatal <- df_merged
# plot
plot_c4_eqtl <- function(data, x_col, title_text, y_label = "Normalized Expression (C4A)") {
  formula_eqtl <- as.formula(paste("Corrected_Expr ~", x_col))
  fit_eqtl <- lm(formula_eqtl, data = data)
  
  slope <- coef(fit_eqtl)[2]
  p_val <- summary(fit_eqtl)$coefficients[2, 4]
  
  title_color <- ifelse(p_val < 0.05, "#d7191c", "black")
  subtitle_text <- sprintf("Slope = %.3f, P = %.2e", slope, p_val)
  
  p <- ggplot(data, aes(x = .data[[x_col]], y = Corrected_Expr)) +
    geom_violin(aes(group = .data[[x_col]]), fill = "#a6cee3", color = NA, alpha = 0.6) +
    geom_boxplot(aes(group = .data[[x_col]]), width = 0.15, fill = "white", color = "black", outlier.shape = NA) +
    geom_jitter(width = 0.1, color = "grey40", size = 1.2, alpha = 0.8) +
    geom_smooth(method = "lm", color = "#b2182b", linetype = "dashed", se = TRUE, fill = "grey80", linewidth = 1) +
    scale_x_continuous(breaks = seq(min(data[[x_col]], na.rm=TRUE), max(data[[x_col]], na.rm=TRUE), by = 1)) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 13, color = title_color),
      axis.text = element_text(color = "black")
    ) +
    labs(title = title_text, subtitle = subtitle_text, x = x_col, y = y_label)
  
  return(p)
}


p1 <- plot_c4_eqtl(df_merged, "Total_C4", "Total C4 CN vs C4A Expr")
p2 <- plot_c4_eqtl(df_merged, "C4A_CN", "C4A CN vs C4A Expr")
p3 <- plot_c4_eqtl(df_merged, "C4B_CN", "C4B CN vs C4A Expr")

final_plot <- p1 | p2 | p3

final_plot <- final_plot + plot_annotation(
  title = "Postnatal Brain: C4A Expression vs Copy Numbers",
  theme = theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))
)

print(final_plot) 

ggsave("Postnata_Brain_C4A_eQTL_Bulk_Corrected.png", final_plot, width = 10, height = 4, dpi = 400)

# Generate Supplementary Table 5: eQTL Linear Regression Stats
get_eqtl_stats <- function(data, stage_name, target_gene) {
  predictors <- c("Total_C4", "C4A_CN", "C4B_CN")
  res_list <- list()
  
  for (pred in predictors) {
    formula_eqtl <- as.formula(paste("Corrected_Expr ~", pred))
    fit_eqtl <- lm(formula_eqtl, data = data)
    
    summary_fit <- summary(fit_eqtl)
    slope <- summary_fit$coefficients[2, "Estimate"]
    std_error <- summary_fit$coefficients[2, "Std. Error"]
    p_val <- summary_fit$coefficients[2, "Pr(>|t|)"]
    
    res_list[[pred]] <- data.frame(
      Developmental_Stage = stage_name,
      Target_Gene = target_gene,
      CNV_Predictor = pred,
      Slope = slope,
      Std_Error = std_error,
      P_value = p_val
    )
  }
  
  return(bind_rows(res_list))
}

stats_prenatal <- get_eqtl_stats(df_prenatal, "Prenatal", "C4A")
stats_postnatal <- get_eqtl_stats(df_postnatal, "Postnatal", "C4A")

supp_table5 <- bind_rows(stats_prenatal, stats_postnatal)
supp_table5 <- supp_table5 %>%
  group_by(Developmental_Stage) %>%
  mutate(FDR = p.adjust(P_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    Slope = round(Slope, 3),
    Std_Error = round(Std_Error, 4),
    P_value = signif(P_value, 3),
    FDR = signif(FDR, 3)
  )
print(supp_table5)

output_csv5 <- "Supplementary_Table_5_eQTL_Stats.csv"
write.csv(supp_table5, output_csv5, row.names = FALSE, quote = FALSE)

# Calculate 95% CI and P-value for Slope Comparison (Prenatal vs Postnatal)

cat("\nCalculating 95% CIs and Slope Comparison P-values...\n")
calculate_slope_comparison <- function(df_pre, df_post, x_col) {
  # 1
  formula_str <- paste("Corrected_Expr ~", x_col)
  fit_pre <- lm(as.formula(formula_str), data = df_pre)
  fit_post <- lm(as.formula(formula_str), data = df_post)
  
  # 2. Slope, SE, 95% CI
  sum_pre <- summary(fit_pre)
  b_pre <- sum_pre$coefficients[2, "Estimate"]
  se_pre <- sum_pre$coefficients[2, "Std. Error"]
  ci_pre <- confint(fit_pre)[2, ]
  
  # 3.Slope, SE, 95% CI
  sum_post <- summary(fit_post)
  b_post <- sum_post$coefficients[2, "Estimate"]
  se_post <- sum_post$coefficients[2, "Std. Error"]
  ci_post <- confint(fit_post)[2, ]
  
  # Z = (Beta_post - Beta_pre) / sqrt(SE_post^2 + SE_pre^2)
  z_stat <- (b_post - b_pre) / sqrt(se_post^2 + se_pre^2)
  p_val_diff <- 2 * (1 - pnorm(abs(z_stat)))
  
  cat(sprintf("Prenatal Slope  = %.3f, 95%% CI = [%.3f, %.3f]\n", b_pre, ci_pre[1], ci_pre[2]))
  cat(sprintf("Postnatal Slope = %.3f, 95%% CI = [%.3f, %.3f]\n", b_post, ci_post[1], ci_post[2]))
  cat(sprintf("Comparison P-value = %.2e\n", p_val_diff))
  
  return(list(ci_pre = ci_pre, ci_post = ci_post, p_diff = p_val_diff))
}

# Prenatal Slope  = 0.153, 95% CI = [0.083, 0.224]
# Postnatal Slope = 0.370, 95% CI = [0.293, 0.447]
# Comparison P-value = 4.49e-05
# Prenatal Slope  = 0.141, 95% CI = [0.070, 0.212]
# Postnatal Slope = 0.325, 95% CI = [0.256, 0.394]
# Comparison P-value = 2.47e-04