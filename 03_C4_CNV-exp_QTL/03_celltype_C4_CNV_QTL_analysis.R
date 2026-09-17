################## Prenatal ##################################
library(tidyverse)
library(data.table)
library(patchwork)


hap_file <- "/gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/topmed_analysis/final_check/fetal_onlychb_imputed_haps_refined.R5.txt"
expr_dir <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/removeversionIDmapping/cts_expr/"
cov_dir  <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/cts-qtl-output-rmgeneversionID/stage1/"


cell_types <- c("prenatal_AST", "prenatal_MG", "prenatal_IN", "prenatal_ExNeu","GLIALPROG","prenatal_OPC")

# ================================================================================
# 2. Read and parse Haplotype (copy num)
# ===================================================================================
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
# 3. Define a generic plotting function (for a single eQTL)
# ==============================================================================
plot_c4_eqtl <- function(data, x_col, title_text, y_label = "Normalized Expr (C4A)") {
  
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
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 12, color = title_color),
      axis.text = element_text(color = "black")
    ) +
    labs(title = title_text, subtitle = subtitle_text, x = x_col, y = y_label)
  
  return(p)
}

# ==============================================================================
# 4. Iterate through each Cell Type: Adjust covariates + Plot data.
# ==============================================================================
stats_list_prenatal <- list()
plot_list <- list()

for (ct in cell_types) {
  cat(sprintf("Processing %s...\n", ct))
  ct_label <- gsub("prenatal_", "", ct)
  cov_file <- paste0(cov_dir, ct, "/covariates.cov.txt")
  cov_raw <- fread(cov_file, header = TRUE)
  

  colnames(cov_raw)[1] <- "Covariate"
  cov_df <- cov_raw %>%
    column_to_rownames("Covariate") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Sample")
  cov_df[, -1] <- lapply(cov_df[, -1], as.numeric)
  expr_file <- paste0(expr_dir, "stage1_", ct, ".txt")
  expr_raw <- fread(expr_file, header = TRUE)
  colnames(expr_raw)[1] <- "Gene"
  c4a_expr <- expr_raw %>%
    filter(Gene == "C4A" | grepl("ENSG00000244731", Gene)) %>%
    column_to_rownames("Gene") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Sample")
  colnames(c4a_expr)[2] <- "Raw_Expr"
  
# --- 4.3 Merging Data and Adjusting Covariates ---

  df_merged <- inner_join(cn_df, c4a_expr, by = "Sample") %>%
    inner_join(cov_df, by = "Sample")
  
  if(nrow(df_merged) == 0) {
    cat("Warning: No matching samples found for", ct, "after inner_join! Skipping...\n")
    next
  }
  
  cov_names <- setdiff(colnames(cov_df), "Sample")
  formula_str <- paste("Raw_Expr ~", paste(cov_names, collapse = " + "))
  
  fit_cov <- lm(as.formula(formula_str), data = df_merged)
  df_merged$Corrected_Expr <- residuals(fit_cov) + mean(df_merged$Raw_Expr, na.rm = TRUE)
  

  p1 <- plot_c4_eqtl(df_merged, "Total_C4", paste0(ct_label, ": Total C4 vs C4A Expr"), y_label = paste(ct_label, "C4A Expr"))
  p2 <- plot_c4_eqtl(df_merged, "C4A_CN",   paste0(ct_label, ": C4A CN vs C4A Expr"), y_label = NULL)
  p3 <- plot_c4_eqtl(df_merged, "C4B_CN",   paste0(ct_label, ": C4B CN vs C4A Expr"), y_label = NULL)
  predictors <- c("Total_C4", "C4A_CN", "C4B_CN")
  ct_stats <- list()
  for(pred in predictors) {
    fit <- lm(as.formula(paste("Corrected_Expr ~", pred)), data = df_merged)
    sum_fit <- summary(fit)
    ct_stats[[pred]] <- data.frame(
      Developmental_Stage = "Prenatal",
      Cell_Type = ct_label,
      Target_Gene = "C4A",
      CNV_Predictor = pred,
      Slope = sum_fit$coefficients[2, "Estimate"],
      Std_Error = sum_fit$coefficients[2, "Std. Error"],
      P_value = sum_fit$coefficients[2, "Pr(>|t|)"]
    )
  }
  stats_list_prenatal[[ct]] <- bind_rows(ct_stats)
  row_plot <- p1 | p2 | p3
  plot_list[[ct]] <- row_plot
}

# ==============================================================================
# 5. Assemble the panoramic image and output
# ==============================================================================
cat("Assembling final combined plot...\n")

final_grid <- wrap_plots(plot_list, ncol = 1) +
  plot_annotation(
    title = "Cell-Type Specific eQTL: C4 Copy Numbers vs C4A Expression (Prenatal Brain)",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))
  )

print(final_grid)

ggsave("prenatal_Brain_C4A_eQTL_CellType_Corrected.png", final_grid, width = 14, height = 18, dpi = 400)
ggsave("Fetal_Brain_C4A_eQTL_CellType_Corrected.pdf", final_grid, width = 14, height = 18)
cat("Done! Saved to Fetal_Brain_C4A_eQTL_CellType_Corrected.png\n")

################## Postnatal ##################################

library(tidyverse)
library(data.table)
library(patchwork)


hap_file <- "/gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/topmed_analysis/final_check/adult/adult_only_imputed_haps_refinedR5revisedID.txt"
expr_dir <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/removeversionIDmapping/cts_expr/"
cov_dir  <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/cts-qtl-output-rmgeneversionID/adult/"


cell_types <- c("postnatal_AST", "postnatal_MG", "postnatal_IN", "postnatal_ExNeu")

# ================================================================================
# 2. Read and parse Haplotype (copy num)
# ===================================================================================
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
# 3. Define a generic plotting function (for a single eQTL)
# ==============================================================================
plot_c4_eqtl <- function(data, x_col, title_text, y_label = "Normalized Expr (C4A)") {
  
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
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 12, color = title_color),
      axis.text = element_text(color = "black")
    ) +
    labs(title = title_text, subtitle = subtitle_text, x = x_col, y = y_label)
  
  return(p)
}

# ==============================================================================
# 4. Iterate through each Cell Type: Adjust covariates + Plot data.
# ==============================================================================
stats_list_postnatal <- list()
plot_list <- list()

for (ct in cell_types) {
  cat(sprintf("Processing %s...\n", ct))
  ct_label <- gsub("postnatal_", "", ct)
  cov_file <- paste0(cov_dir, ct, "/covariates.cov.txt")
  cov_raw <- fread(cov_file, header = TRUE)
  

  colnames(cov_raw)[1] <- "Covariate"
  cov_df <- cov_raw %>%
    column_to_rownames("Covariate") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Sample")
  cov_df[, -1] <- lapply(cov_df[, -1], as.numeric)
  expr_file <- paste0(expr_dir, "adult_", ct, ".txt")
  expr_raw <- fread(expr_file, header = TRUE)
  colnames(expr_raw)[1] <- "Gene"
  c4a_expr <- expr_raw %>%
    filter(Gene == "C4A" | grepl("ENSG00000244731", Gene)) %>%
    column_to_rownames("Gene") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Sample")
  colnames(c4a_expr)[2] <- "Raw_Expr"
  
# --- 4.3 Merging Data and Adjusting Covariates ---

  df_merged <- inner_join(cn_df, c4a_expr, by = "Sample") %>%
    inner_join(cov_df, by = "Sample")
  
  if(nrow(df_merged) == 0) {
    cat("Warning: No matching samples found for", ct, "after inner_join! Skipping...\n")
    next
  }
  
  cov_names <- setdiff(colnames(cov_df), "Sample")
  formula_str <- paste("Raw_Expr ~", paste(cov_names, collapse = " + "))
  
  fit_cov <- lm(as.formula(formula_str), data = df_merged)
  df_merged$Corrected_Expr <- residuals(fit_cov) + mean(df_merged$Raw_Expr, na.rm = TRUE)
  

  p1 <- plot_c4_eqtl(df_merged, "Total_C4", paste0(ct_label, ": Total C4 vs C4A Expr"), y_label = paste(ct_label, "C4A Expr"))
  p2 <- plot_c4_eqtl(df_merged, "C4A_CN",   paste0(ct_label, ": C4A CN vs C4A Expr"), y_label = NULL)
  p3 <- plot_c4_eqtl(df_merged, "C4B_CN",   paste0(ct_label, ": C4B CN vs C4A Expr"), y_label = NULL)
  
predictors <- c("Total_C4", "C4A_CN", "C4B_CN")
  ct_stats <- list()
  for(pred in predictors) {
    fit <- lm(as.formula(paste("Corrected_Expr ~", pred)), data = df_merged)
    sum_fit <- summary(fit)
    ct_stats[[pred]] <- data.frame(
      Developmental_Stage = "Postnatal",
      Cell_Type = ct_label,
      Target_Gene = "C4A",
      CNV_Predictor = pred,
      Slope = sum_fit$coefficients[2, "Estimate"],
      Std_Error = sum_fit$coefficients[2, "Std. Error"],
      P_value = sum_fit$coefficients[2, "Pr(>|t|)"]
    )
  }
  stats_list_postnatal[[ct]] <- bind_rows(ct_stats)

  row_plot <- p1 | p2 | p3
  plot_list[[ct]] <- row_plot
}

# ==============================================================================
# 5. Assemble the panoramic image and output
# ==============================================================================
cat("Assembling final combined plot...\n")

final_grid <- wrap_plots(plot_list, ncol = 1) +
  plot_annotation(
    title = "Cell-Type Specific eQTL: C4 Copy Numbers vs C4A Expression (Postnatal Brain)",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))
  )

print(final_grid)

ggsave("postnatal_Brain_C4A_eQTL_CellType_Corrected2.png", final_grid, width = 14, height = 18, dpi = 400)
ggsave("Adult_Brain_C4A_eQTL_CellType_Corrected2.pdf", final_grid, width = 14, height = 18)

# 6. Generate Supplementary Table 6: Cell-Type Specific eQTL Stats
supp_table6 <- bind_rows(
  bind_rows(stats_list_prenatal),
  bind_rows(stats_list_postnatal)
)

# Calculate FDR (multiple test correction by developmental stage grouping)
supp_table6 <- supp_table6 %>%
  group_by(Developmental_Stage) %>%
  mutate(FDR = p.adjust(P_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    Slope = round(Slope, 3),
    Std_Error = round(Std_Error, 4),
    P_value = signif(P_value, 3),
    FDR = signif(FDR, 3)
  )

print(as.data.frame(supp_table6))
output_csv6 <- "Supplementary_Table_6_CellType_eQTL_Stats.csv"
write.csv(supp_table6, output_csv6, row.names = FALSE, quote = FALSE)