    #B2182B #CFA9C5
     #2166AC #96D1CD
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(patchwork)
})

plot_sex_interaction_eqtl <- function(data, x_col, title_text, y_label = "Corrected Expression (C4A)") {
# Likelihood Ratio Test / ANOVA
  formula_null <- as.formula(paste("Corrected_Expr ~", x_col, "+ sex_factor"))
  formula_int  <- as.formula(paste("Corrected_Expr ~", x_col, "* sex_factor"))
  
  fit_null <- lm(formula_null, data = data)
  fit_int  <- lm(formula_int, data = data)
  int_p_val <- anova(fit_null, fit_int)[2, "Pr(>F)"]
  
 # 2. Calculate gender-stratified independent slopes and p-values
  data_m <- data %>% filter(sex_factor == "Male")
  data_f <- data %>% filter(sex_factor == "Female")
  
  fit_m <- lm(as.formula(paste("Corrected_Expr ~", x_col)), data = data_m)
  fit_f <- lm(as.formula(paste("Corrected_Expr ~", x_col)), data = data_f)
  
  slope_m <- coef(fit_m)[2]; p_m <- summary(fit_m)$coefficients[2, 4]
  slope_f <- coef(fit_f)[2]; p_f <- summary(fit_f)$coefficients[2, 4]

  int_color <- ifelse(int_p_val < 0.05, "black", "black")
  sub_text <- sprintf("Interaction P = %.2e\nMale: Slope=%.3f (P=%.2e)\nFemale: Slope=%.3f (P=%.2e)", 
                      int_p_val, slope_m, p_m, slope_f, p_f)
  
  #plot
  p <- ggplot(data, aes(x = .data[[x_col]], y = Corrected_Expr, color = sex_factor, fill = sex_factor)) +
    geom_boxplot(aes(group = interaction(.data[[x_col]], sex_factor)), 
                 width = 0.4, position = position_dodge(width = 0.7), 
                 alpha = 0.3, outlier.shape = NA, color = "grey50") +
    geom_point(position = position_jitterdodge(jitter.width = 0.7, dodge.width = 0.7), 
               size = 1.2, alpha = 0.95) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 1) +
    scale_color_manual(values = c("Male" = "#96D1CD", "Female" = "#CFA9C5")) +
    scale_fill_manual(values = c("Male" = "#96D1CD", "Female" = "#CFA9C5")) +
    scale_x_continuous(breaks = seq(min(data[[x_col]], na.rm=TRUE), max(data[[x_col]], na.rm=TRUE), by = 1)) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      plot.subtitle = element_text(hjust = 0.5, size = 11, lineheight = 1.2, color = int_color),
      axis.text = element_text(color = "black"),
      legend.title = element_blank()
    ) +
    labs(title = title_text, subtitle = sub_text, x = x_col, y = y_label)
  
  return(p)
}


run_analysis <- function(hap_file, expr_file, cov_file, stage_title, out_file) {
  haps <- fread(hap_file, header = TRUE)
  get_cn <- function(hap_vec, pos) as.numeric(sapply(strsplit(hap_vec, "_"), `[`, pos))
  cn_df <- data.frame(
    Sample = haps$SAMPLE,
    Total_C4 = get_cn(haps$H1, 2) + get_cn(haps$H2, 2),
    C4A_CN   = get_cn(haps$H1, 3) + get_cn(haps$H2, 3),
    C4B_CN   = get_cn(haps$H1, 4) + get_cn(haps$H2, 4)
  )
  
  expr_data <- fread(expr_file, header = TRUE)
  colnames(expr_data)[1] <- "Gene"
  c4a_expr <- expr_data %>%
    filter(Gene == "C4A" | grepl("ENSG00000244731", Gene)) %>%
    head(1) %>% column_to_rownames("Gene") %>% t() %>%
    as.data.frame() %>% rownames_to_column("Sample")
  colnames(c4a_expr)[2] <- "Raw_Expr"
  
  cov_raw <- fread(cov_file, header = TRUE)
  colnames(cov_raw)[1] <- "Covariate"
  cov_df <- cov_raw %>%
    column_to_rownames("Covariate") %>% t() %>%
    as.data.frame() %>% rownames_to_column("Sample") %>%
    mutate(across(-Sample, as.numeric))
  
  df_merged <- inner_join(cn_df, c4a_expr, by = "Sample") %>%
    inner_join(cov_df, by = "Sample") %>%
    mutate(sex_factor = factor(sex, levels = c(1, 2), labels = c("Male", "Female")))
  
# Regress out covariates other than sex
  cov_names <- setdiff(colnames(cov_df), c("Sample", "sex"))
  formula_str <- paste("Raw_Expr ~", paste(cov_names, collapse = " + "))
  fit_cov <- lm(as.formula(formula_str), data = df_merged)
  df_merged$Corrected_Expr <- residuals(fit_cov) + mean(df_merged$Raw_Expr, na.rm = TRUE)
  
  p1 <- plot_sex_interaction_eqtl(df_merged, "Total_C4", "Total C4 CN vs C4A Expr")
  p2 <- plot_sex_interaction_eqtl(df_merged, "C4A_CN", "C4A CN vs C4A Expr")
  p3 <- plot_sex_interaction_eqtl(df_merged, "C4B_CN", "C4B CN vs C4A Expr")

  final_plot <- (p1 | p2 | p3) + 
    plot_layout(guides = "collect") +
    plot_annotation(
      title = stage_title,
      theme = theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
                    legend.position = "bottom")
    ) & theme(legend.position = "bottom")
  
  print(final_plot)
  ggsave(out_file, final_plot, width = 11, height = 5.5, dpi = 400)
}

# Prenatal

run_analysis(
  hap_file  = "fetal_onlychb_imputed_haps_refined.R5.txt",
  expr_file = "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/bulk_data/stage1.logTMM.ComBat.txt",
  cov_file  = "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/bulk_data/prenatalcovariatesToUse.PCAforQTLfromINT.txt",
  stage_title = "Prenatal Brain: Sex Interaction eQTL (C4A Expression)",
  out_file  = "Fetal_Brain_C4A_Sex_Interaction_eQTL.png"
)

# Postnatal

run_analysis(
  hap_file  = "adult_only_imputed_haps_refinedR5revisedID.txt",
  expr_file = "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/bulk_data/postnatal.logTMM.ComBat.txt",
  cov_file  = "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/bulk_data/postnatalcovariatesToUse.PCAforQTLfromINT.txt",
  stage_title = "Postnatal Brain: Sex Interaction eQTL (C4A Expression)",
  out_file  = "Postnatal_Brain_C4A_Sex_Interaction_eQTL.png"
)