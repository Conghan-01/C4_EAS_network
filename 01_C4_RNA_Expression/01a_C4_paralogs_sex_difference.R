# ==============================================================================
# Title: Comparative Sex Difference Analysis of C4A and C4B
# Methods: 1. Known Covariate (Age/RIN) Correction 
#          2. SVA (Surrogate Variable Analysis) Latent Correction

library(tidyverse)
library(data.table)
library(ggpubr)
library(sva)

# Define file paths for expression matrices and covariate files
pre_expr_file  <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/bulk_data/stage1.logTMM.ComBat.txt"
post_expr_file <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/bulk_data/postnatal.expr.qn.comb"
pre_cov_file   <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/bulk_data/prenatalcovariatesToUse.PCAforQTLfromINT.txt"
post_cov_file  <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/bulk_data/postnatalcovariatesToUse.PCAforQTLfromINT.txt"

# Map Ensembl IDs to Gene Symbols
target_genes <- c("C4A" = "ENSG00000244731", "C4B" = "ENSG00000224389")
sex_colors <- c("Female" = "#CFA9C5", "Male" = "#96D1CD")

# Process Prenatal Data
expr_pre <- fread(pre_expr_file) %>% 
            mutate(V1 = gsub("\\..*$", "", V1)) %>% 
            filter(V1 %in% target_genes) %>% 
            column_to_rownames("V1") %>% t() %>% as.data.frame() %>% 
            rownames_to_column("Sample")

cov_pre  <- fread(pre_cov_file) %>% 
            column_to_rownames("cov") %>% t() %>% as.data.frame() %>% 
            rownames_to_column("Sample")

df_pre   <- inner_join(expr_pre, cov_pre, by = "Sample") %>% mutate(Phase = "Prenatal")
df_pre$Age<- df_pre$GW
df_pre<- df_pre[,c("Sample","ENSG00000224389","ENSG00000244731","RIN","sex","Age","Phase")]
# Process Postnatal Data
expr_post <- fread(post_expr_file) %>% 
             mutate(V1 = gsub("\\..*$", "", V1)) %>% 
             filter(V1 %in% target_genes) %>% 
             column_to_rownames("V1") %>% t() %>% as.data.frame() %>% 
             rownames_to_column("Sample")

cov_post  <- fread(post_cov_file) %>% 
             column_to_rownames("cov") %>% t() %>% as.data.frame() %>% 
             rownames_to_column("Sample")

df_post   <- inner_join(expr_post, cov_post, by = "Sample") %>% mutate(Phase = "Postnatal")
df_post<- df_post[,c("Sample","ENSG00000224389","ENSG00000244731","RIN","sex","Age","Phase")]
# Combine into a single raw dataframe for Method 1
full_raw_df <- bind_rows(df_pre, df_post) %>% 
  rename(C4A = ENSG00000244731, C4B = ENSG00000224389) %>%
  pivot_longer(cols = c("C4A", "C4B"), names_to = "Gene", values_to = "Raw")

#  Method 1: Known Covariate Correction (Age + RIN) ---

# This step uses linear regression to remove Age/RIN effects while preserving Sex variance
results_known <- full_raw_df %>%
  group_by(Phase, Gene) %>%
  group_modify(~ {
    fit <- lm(Raw ~ Age + RIN, data = .x)
    # Add residuals back to the mean to maintain original scale (logTMM)
    .x$Corrected_Value <- residuals(fit) + mean(.x$Raw, na.rm = TRUE)
    return(.x)
  }) %>% ungroup()

# Method 2: SVA (Surrogate Variable Analysis) ---

# Prenatal SVA Processing
mat_pre <- fread(pre_expr_file) %>% column_to_rownames("V1") %>% as.matrix()
mat_pre <- mat_pre[rowMeans(mat_pre) > 0.1, ] # Filter low-expression genes for SVA stability
meta_pre <- cov_pre; meta_pre$Age <- meta_pre$GW

# Design matrices: mod includes sex, mod0 (null) does not
mod_pre  <- model.matrix(~ sex + Age + RIN, data = meta_pre)
mod0_pre <- model.matrix(~ Age + RIN, data = meta_pre)
n_sv_pre <- num.sv(mat_pre, mod_pre, method = "be")
sv_pre   <- sva(mat_pre, mod_pre, mod0_pre, n.sv = n_sv_pre)

# Incorporate SVs into the Prenatal dataframe
sv_pre_df <- as.data.frame(sv_pre$sv); colnames(sv_pre_df) <- paste0("SV", 1:n_sv_pre)
df_pre_sva <- df_pre %>% rename(C4A = ENSG00000244731, C4B = ENSG00000224389) %>% 
              cbind(sv_pre_df) %>% pivot_longer(cols = c("C4A", "C4B"), names_to = "Gene", values_to = "Raw")

# Postnatal SVA Processing
mat_post <- fread(post_expr_file) %>% column_to_rownames("V1") %>% as.matrix()
mat_post <- mat_post[rowMeans(mat_post) > 0.1, ]
meta_post <- cov_post; meta_post$sex <- as.factor(meta_post$sex)

mod_post  <- model.matrix(~ sex + age + rin, data = meta_post)
mod0_post <- model.matrix(~ age + rin, data = meta_post)
n_sv_post <- num.sv(mat_post, mod_post, method = "be")
sv_post   <- sva(mat_post, mod_post, mod0_post, n.sv = n_sv_post)

# Incorporate SVs into the Postnatal dataframe
sv_post_df <- as.data.frame(sv_post$sv); colnames(sv_post_df) <- paste0("SV", 1:n_sv_post)
df_post_sva <- df_post %>% rename(C4A = ENSG00000244731, C4B = ENSG00000224389) %>% 
               cbind(sv_post_df) %>% pivot_longer(cols = c("C4A", "C4B"), names_to = "Gene", values_to = "Raw")

#Residual Regression including SVs
# Prenatal Regression
formula_pre <- as.formula(paste("Raw ~ age + rin +", paste(colnames(sv_pre_df), collapse = "+")))
results_sva_pre <- df_pre_sva %>% group_by(Gene) %>% group_modify(~ {
  fit <- lm(formula_pre, data = .x)
  .x$Corrected_Value <- residuals(fit) + mean(.x$Raw, na.rm = TRUE)
  return(.x)
})

# Postnatal Regression
formula_post <- as.formula(paste("Raw ~ age + rin +", paste(colnames(sv_post_df), collapse = "+")))
results_sva_post <- df_post_sva %>% group_by(Gene) %>% group_modify(~ {
  fit <- lm(formula_post, data = .x)
  .x$Corrected_Value <- residuals(fit) + mean(.x$Raw, na.rm = TRUE)
  return(.x)
})

results_sva <- bind_rows(results_sva_pre, results_sva_post) %>% ungroup()


# Function to map sex labels and set factor levels
prepare_plot <- function(dat) {
  dat %>% mutate(
    Sex_Label = ifelse(sex == 2, "Female", "Male"),
    Sex_Label = factor(Sex_Label, levels = c("Female", "Male")),
    Phase = factor(Phase, levels = c("Prenatal", "Postnatal"))
  )
}

plot_known_data <- prepare_plot(results_known)
plot_sva_data   <- prepare_plot(results_sva)

# Plot 1: Age & RIN Corrected
p1 <- ggplot(plot_known_data, aes(x = Phase, y = Corrected_Value, fill = Sex_Label)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.8, color = "black", position = position_dodge(width = 0.7)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.7), size = 0.8, alpha = 0.3, color = "black") +
  facet_wrap(~Gene, scales = "free_y") +
  stat_compare_means(aes(group = Sex_Label), method = "wilcox.test", label = "p.signif", size =5) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  scale_fill_manual(values = sex_colors) +
  theme_bw(base_size = 17) +
  labs(title = "Sex Differences (Age & RIN Corrected)", x = "Developmental Stages", y = "Corrected Expression (logTMM)") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "top", legend.title = element_blank())

# Plot 2: SVA Corrected
p2 <- ggplot(plot_sva_data, aes(x = Phase, y = Corrected_Value, fill = Sex_Label)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.8, color = "black", position = position_dodge(width = 0.7)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.7), size = 0.8, alpha = 0.3, color = "black") +
  facet_wrap(~Gene, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  stat_compare_means(aes(group = Sex_Label), method = "wilcox.test", label = "p.signif", size = 5) +
  scale_fill_manual(values = sex_colors) +
  theme_bw(base_size = 17) +
  labs(title = "Sex Differences (SVA Corrected)", x = "Developmental Stages", y = "Corrected Expression (logTMM)") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "top", legend.title = element_blank())


ggsave("C4_Sex_Diff_Known.png", p1, width = 6.7, height = 5, dpi = 500)
ggsave("C4_Sex_Diff_SVA.png", p2, width = 6.7, height = 5, dpi = 500)


cat("\n--- Known Correction Statistics ---\n")
print(plot_known_data %>% group_by(Phase, Gene) %>% summarise(p_val = wilcox.test(Corrected_Value ~ sex)$p.value))

cat("\n--- SVA Correction Statistics ---\n")
print(plot_sva_data %>% group_by(Phase, Gene) %>% summarise(p_val = wilcox.test(Corrected_Value ~ sex)$p.value))