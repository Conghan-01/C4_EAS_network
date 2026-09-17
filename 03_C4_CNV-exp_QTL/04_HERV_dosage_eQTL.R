library(tidyverse)
library(data.table)
library(patchwork)


parse_herv_dosage <- function(hap_file) {
  haps <- fread(hap_file, header = TRUE)

  get_part <- function(hap_vec, pos) {
    as.numeric(sapply(strsplit(as.character(hap_vec), "_"), `[`, pos))
  }
  
  df_res <- data.frame(
    Sample = haps$SAMPLE,
    Total_C4    = get_part(haps$H1, 2) + get_part(haps$H2, 2),
    C4A_CN      = get_part(haps$H1, 3) + get_part(haps$H2, 3),
    C4B_CN      = get_part(haps$H1, 4) + get_part(haps$H2, 4),
    HERV_Dosage = get_part(haps$H1, 5) + get_part(haps$H2, 5)
  )
  
  return(df_res)
}

process_data_and_correct_filtered <- function(hap_file, expr_file, cov_file, target_gene_symbol, target_ensembl_id, max_herv = 6) {
  cn_df <- parse_herv_dosage(hap_file)
  expr_data <- fread(expr_file, header = TRUE)
  colnames(expr_data)[1] <- "Gene"
  
  gene_expr <- expr_data %>%
    filter(Gene == target_gene_symbol | grepl(target_ensembl_id, Gene)) %>%
    head(1) %>% 
    column_to_rownames("Gene") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Sample")
  
  colnames(gene_expr)[2] <- "Raw_Expr"
  cov_raw <- fread(cov_file, header = TRUE)
  colnames(cov_raw)[1] <- "Covariate"
  
  cov_df <- cov_raw %>%
    column_to_rownames("Covariate") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Sample") %>%
    mutate(across(-Sample, as.numeric))
  df_merged <- inner_join(cn_df, gene_expr, by = "Sample") %>%
    inner_join(cov_df, by = "Sample")
  cov_names <- setdiff(colnames(cov_df), "Sample")
  formula_str <- paste("Raw_Expr ~", paste(cov_names, collapse = " + "))
  fit_cov <- lm(as.formula(formula_str), data = df_merged)
  df_merged$Corrected_Expr <- residuals(fit_cov) + mean(df_merged$Raw_Expr, na.rm = TRUE)
  df_filtered <- df_merged %>%
    filter(HERV_Dosage <= max_herv) %>% 
    group_by(HERV_Dosage) %>%
    filter(n() >= 2) %>%               
    ungroup()
  return(df_filtered)
}

plot_herv_eqtl <- function(data, title_text, y_label = "Corrected Expression") {
  fit_eqtl <- lm(Corrected_Expr ~ HERV_Dosage, data = data)
  slope <- coef(fit_eqtl)[2]
  p_val <- summary(fit_eqtl)$coefficients[2, 4]
  
  title_color <- ifelse(p_val < 0.05, "#d7191c", "black")
  subtitle_text <- sprintf("Slope = %.3f, P = %.2e", slope, p_val)
  
  ggplot(data, aes(x = HERV_Dosage, y = Corrected_Expr)) +
    geom_violin(aes(group = HERV_Dosage), fill = "#b2df8a", color = NA, alpha = 0.6) +
    geom_boxplot(aes(group = HERV_Dosage), width = 0.15, fill = "white", color = "black", outlier.shape = NA) +
    geom_jitter(width = 0.1, color = "grey40", size = 1.2, alpha = 0.8) +
    geom_smooth(method = "lm", color = "#33a02c", linetype = "dashed", se = TRUE, fill = "grey80", linewidth = 1) +
    scale_x_continuous(breaks = seq(min(data$HERV_Dosage, na.rm=TRUE), max(data$HERV_Dosage, na.rm=TRUE), by = 1)) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 12, color = title_color),
      axis.text = element_text(color = "black")
    ) +
    labs(title = title_text, subtitle = subtitle_text, x = "C4-HERV Copy Number", y = y_label)
}

cat("Processing Prenatal Cohort for HERV Dosage eQTL...\n")

pre_hap <- "/gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/topmed_analysis/final_check/fetal_onlychb_imputed_haps_refined.R5.txt"
pre_expr <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/bulk_data/stage1.logTMM.ComBat.txt"
pre_cov <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/bulk_data/prenatalcovariatesToUse.PCAforQTLfromINT.txt"

df_pre_c4a <- process_data_and_correct_filtered(pre_hap, pre_expr, pre_cov, "C4A", "ENSG00000244731")
df_pre_c4b <- process_data_and_correct_filtered(pre_hap, pre_expr, pre_cov, "C4B", "ENSG00000224389")

p_pre_c4a <- plot_herv_eqtl(df_pre_c4a, "Prenatal: HERV Dosage vs C4A Expr", "C4A Expression")
p_pre_c4b <- plot_herv_eqtl(df_pre_c4b, "Prenatal: HERV Dosage vs C4B Expr", "C4B Expression")


# 2.  Postnatal

cat("Processing Postnatal Cohort for HERV Dosage eQTL...\n")

post_hap <- "/gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/topmed_analysis/final_check/adult/adult_only_imputed_haps_refinedR5revisedID.txt"
post_expr <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/bulk_data/postnatal.logTMM.ComBat.txt"
post_cov <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/bulk_data/postnatalcovariatesToUse.PCAforQTLfromINT.txt"

df_post_c4a <- process_data_and_correct_filtered(post_hap, post_expr, post_cov, "C4A", "ENSG00000244731")
df_post_c4b <- process_data_and_correct_filtered(post_hap, post_expr, post_cov, "C4B", "ENSG00000224389")

p_post_c4a <- plot_herv_eqtl(df_post_c4a, "Postnatal: HERV Dosage vs C4A Expr", "C4A Expression")
p_post_c4b <- plot_herv_eqtl(df_post_c4b, "Postnatal: HERV Dosage vs C4B Expr", "C4B Expression")


final_herv_plot <- (p_pre_c4a | p_pre_c4b) / (p_post_c4a | p_post_c4b) +
  plot_annotation(
    title = "C4-HERV Dosage Effects on C4A and C4B Gene Expression Across Development",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )

print(final_herv_plot)
setwd("/gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/figures/fig1")
ggsave("C4_HERV_Dosage_eQTL_Developmental1.png", final_herv_plot, width = 10, height = 8, dpi = 400)


get_herv_stats <- function(df, stage, gene_name) {
  fit_marginal <- lm(Corrected_Expr ~ HERV_Dosage, data = df)
  sum_m <- summary(fit_marginal)
  fit_joint <- lm(Corrected_Expr ~ HERV_Dosage + C4A_CN, data = df)
  sum_j <- summary(fit_joint)
  
  data.frame(
    Stage = stage,
    Gene = gene_name,
    Marginal_HERV_Slope = sum_m$coefficients["HERV_Dosage", "Estimate"],
    Marginal_HERV_SE    = sum_m$coefficients["HERV_Dosage", "Std. Error"],
    Marginal_HERV_P     = sum_m$coefficients["HERV_Dosage", "Pr(>|t|)"],
    Joint_HERV_Slope    = sum_j$coefficients["HERV_Dosage", "Estimate"],
    Joint_HERV_P        = sum_j$coefficients["HERV_Dosage", "Pr(>|t|)"]
  )
}

herv_stats <- bind_rows(
  get_herv_stats(df_pre_c4a, "Prenatal", "C4A"),
  get_herv_stats(df_pre_c4b, "Prenatal", "C4B"),
  get_herv_stats(df_post_c4a, "Postnatal", "C4A"),
  get_herv_stats(df_post_c4b, "Postnatal", "C4B")
)

output_csv_herv <- "Supplementary_Table_HERV_eQTL_Stats.csv"
write.csv(herv_stats, output_csv_herv, row.names = FALSE, quote = FALSE)

cat("\nHERV Dosage eQTL Analysis Complete! Output saved to Supplementary_Table_HERV_eQTL_Stats.csv\n")