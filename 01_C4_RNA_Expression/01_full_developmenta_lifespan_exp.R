##################     EAS        #####################################
#######################################################################
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
})

setwd("/gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/figures/fig1")
meta_path <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/01mixbatches_group/meta.txt2"
tpm_path  <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/01mixbatches_group/subset_tpm.txt"
meta <- read.table(meta_path, header = TRUE, stringsAsFactors = FALSE)
tpm <- read.table(tpm_path, header = TRUE, row.names = 1, check.names = FALSE)


target_genes <- c("C4A" = "ENSG00000244731.9", 
                  "C4B" = "ENSG00000224389.9")

valid_genes <- target_genes[target_genes %in% rownames(tpm)]
if(length(valid_genes) == 0) stop("Error: Target gene IDs not found in TPM matrix.")

tpm_t <- as.data.frame(t(tpm[valid_genes, , drop = FALSE]))
colnames(tpm_t) <- names(valid_genes)
tpm_t$sample <- rownames(tpm_t)

# 3. Merge and perform [developmental stage classification] (Infancy excluded)
df_plot <- merge(meta, tpm_t, by = "sample") %>%
  pivot_longer(cols = all_of(names(valid_genes)), names_to = "Gene", values_to = "TPM") %>%
  mutate(
    # Adjust prenatal staging thresholds based on the conversion between PCW and GW (GW = PCW + 2)
    Period = case_when(
      Stage == "stage1" & age < 15 ~ "Early, prenatal",                   
      Stage == "stage1" & age >= 15 & age < 26 ~ "Middle, prenatal",     
      Stage == "stage1" & age >= 26 ~ "Late, prenatal",                   
      # Omit Infancy (< 1)
      Stage %in% c("stage2", "stage3") & age >= 1 & age < 12 ~ "Childhood",
      Stage %in% c("stage2", "stage3") & age >= 12 & age < 20 ~ "Adolescence",
      Stage %in% c("stage2", "stage3") & age >= 20 & age < 40 ~ "Young adulthood",
      Stage %in% c("stage2", "stage3") & age >= 40 & age < 60 ~ "Middle adulthood",
      Stage %in% c("stage2", "stage3") & age >= 60 ~ "Late adulthood"
    ),
    Period = factor(Period, levels = c(
      "Early, prenatal", "Middle, prenatal", "Late, prenatal",
      "Childhood", "Adolescence", "Young adulthood", "Middle adulthood", "Late adulthood"
    )),
    Period_Num = as.numeric(Period)
  ) %>%
  filter(!is.na(Period)) 

# 4. Calculate the actual sample size (N) for each period
sample_counts <- df_plot %>%
  filter(Gene == names(valid_genes)[1]) %>% 
  group_by(Period_Num, Period) %>%
  summarise(N = n(), .groups = "drop") %>%
  mutate(Label_with_N = paste0(Period, "\n(n=", N, ")")) 


# ---------------------------------------------------------
# 5.plot for C4A and C4B
# ---------------------------------------------------------
plot_gene_trajectory <- function(gene_name, box_color, y_label = "TPM") {
  df_sub <- df_plot %>% filter(Gene == gene_name)
  p <- ggplot(df_sub, aes(x = Period_Num, y = TPM)) + 
    geom_smooth(method = "loess", se = TRUE, 
                color = "#4169E1", fill = "grey85", alpha = 0.4, span = 0.8) +
    geom_boxplot(aes(group = Period_Num), 
                 width = 0.35, fill = box_color, color = "black", 
                 outlier.shape = NA) +
    scale_x_continuous(breaks = sample_counts$Period_Num, 
                       labels = sample_counts$Label_with_N) +
                       coord_cartesian(ylim = c(0, 50)) +
    theme_classic(base_size = 22) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 18, lineheight = 0.8),
      axis.text.y = element_text(color = "black", size = 20),
      axis.title = element_text(face = "bold", size = 22),
      axis.line = element_line(linewidth = 1),
      axis.ticks = element_line(linewidth = 1),
      plot.title = element_text(face = "bold", size = 24, hjust = 0.5),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 40)
    ) +
    labs(title = gene_name, x = NULL, y = y_label)
  return(p)
}
# Save C4A
p_C4A <- plot_gene_trajectory("C4A", "#E69191")
ggsave("C4A_expression_boxplot_lifespantrajectory_new.png", plot = p_C4A, width = 8, height = 7, dpi = 500)
# Save C4B
p_C4B <- plot_gene_trajectory("C4B", "#92B5CA")
ggsave("C4B_expression_boxplot_lifespantrajectory_new.png", plot = p_C4B, width = 8, height = 7, dpi = 500)

cat("C4A and C4B separate plots successfully saved!\n")

##################     EUR        #####################################
#######################################################################

setwd("/gpfs/hpc/home/chenchao/hanc/project/04psychENCODE/addregion")
load("combat.RData")


suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
})


#meta_BrainSpan <- metaData[metaData$study == "BrainSpan" & metaData$tissue == "frontal cortex", ]

meta_BrainSpan <- droplevels(meta_BrainSpan)
meta_BrainSpan$sample <- rownames(meta_BrainSpan) 
common_samples <- intersect(rownames(meta_BrainSpan), colnames(tpm))
tpm_BrainSpan <- tpm[, common_samples]

# 2. Extract target genes (remove version numbers to match BrainSpan's Ensembl ID format)
target_genes <- c("C4A" = "ENSG00000244731", 
                  "C4B" = "ENSG00000224389")

rownames(tpm_BrainSpan) <- gsub("\\..*$", "", rownames(tpm_BrainSpan))
valid_genes <- target_genes[target_genes %in% rownames(tpm_BrainSpan)]
if(length(valid_genes) == 0) stop("Error: Target gene IDs not found in TPM matrix.")

tpm_t <- as.data.frame(t(tpm_BrainSpan[valid_genes, , drop = FALSE]))
colnames(tpm_t) <- names(valid_genes)
tpm_t$sample <- rownames(tpm_t)

df_merged <- as.data.frame(merge(meta_BrainSpan, tpm_t, by = "sample"))

df_plot <- df_merged %>%
  pivot_longer(cols = all_of(names(valid_genes)), names_to = "Gene", values_to = "TPM") %>%
  mutate(
    age_str = as.character(ageDeath),
    is_pcw = grepl("PCW", age_str, ignore.case = TRUE),
    num_val = suppressWarnings(as.numeric(gsub("[^0-9\\.\\-]", "", age_str))),
    
    # Standardize to gestational weeks (GW) and years (Years)
    GW = case_when(
      is_pcw ~ num_val + 2,                                      
      !is_pcw & !is.na(num_val) & num_val < 0 ~ 40 + (num_val * 52.1775), 
      TRUE ~ NA_real_
    ),
    Years = case_when(
      !is_pcw & !is.na(num_val) & num_val >= 0 ~ num_val,        
      TRUE ~ NA_real_
    )
  ) %>%
  mutate(
    Period = case_when(
      !is.na(GW) & GW < 15 ~ "Early, prenatal",                   # < 13 PCW
      !is.na(GW) & GW >= 15 & GW < 21 ~ "Early, mid-prenatal",    # 13-19 PCW
      !is.na(GW) & GW >= 21 & GW < 26 ~ "Late, mid-prenatal",     # 19-24 PCW
      !is.na(GW) & GW >= 26 ~ "Late, prenatal",                   # >= 24 PCW
      !is.na(Years) & Years >= 1 & Years < 12 ~ "Childhood",      # 1-12 Y
      !is.na(Years) & Years >= 12 & Years < 20 ~ "Adolescence",   # 12-20 Y
      !is.na(Years) & Years >= 20 ~ "Adulthood"                  
    ),
    Period = factor(Period, levels = c(
      "Early, prenatal", "Early, mid-prenatal", "Late, mid-prenatal", "Late, prenatal",
      "Childhood", "Adolescence", "Adulthood"
    )),
    Period_Num = as.numeric(Period)
  ) %>%
  filter(!is.na(Period)) 
sample_counts <- df_plot %>%
  filter(Gene == names(valid_genes)[1]) %>% 
  group_by(Period_Num, Period) %>%
  summarise(N = n(), .groups = "drop") %>%
  mutate(Label_with_N = paste0(Period, "\n(n=", N, ")")) 


plot_gene_trajectory <- function(gene_name, box_color, y_label = "TPM") {
  
  df_sub <- df_plot %>% filter(Gene == gene_name)
  p <- ggplot(df_sub, aes(x = Period_Num, y = TPM)) +  
    geom_smooth(method = "loess", se = TRUE, 
                color = "#4169E1", fill = "grey85", alpha = 0.4, span = 1.0) +              
    geom_boxplot(aes(group = Period_Num), 
                 width = 0.35, fill = box_color, color = "black", 
                 outlier.shape = NA) +               
    scale_x_continuous(breaks = sample_counts$Period_Num, 
                       labels = sample_counts$Label_with_N) +
                       
    coord_cartesian(ylim = c(0, 50)) +                     
    theme_classic(base_size = 22) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 18, lineheight = 0.8),
      axis.text.y = element_text(color = "black", size = 20),
      axis.title = element_text(face = "bold", size = 22),
      axis.line = element_line(linewidth = 1),
      axis.ticks = element_line(linewidth = 1),
      plot.title = element_text(face = "bold", size = 24, hjust = 0.5),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 40)
    ) +
    labs(title = paste0("BrainSpan: ", gene_name), x = NULL, y = y_label)
  
  return(p)
}

p_C4A <- plot_gene_trajectory("C4A", "#E69191")
ggsave("BrainSpan_Frontal_C4A_trajectory.pdf", plot = p_C4A, width = 8, height = 7)
ggsave("/gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/figures/fig1/BrainSpan_Frontal_C4A_trajectory.png", plot = p_C4A, width = 8, height = 7, dpi = 500)

p_C4B <- plot_gene_trajectory("C4B", "#92B5CA")
ggsave("BrainSpan_Frontal_C4B_trajectory.pdf", plot = p_C4B, width = 8, height = 7)
ggsave("/gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/figures/fig1/BrainSpan_Frontal_C4B_trajectory.png", plot = p_C4B, width = 8, height = 7, dpi = 500)


library(ggplot2)
library(patchwork)
p1 <- ggplot(mtcars, aes(wt, mpg)) + geom_point() + ggtitle("Plot 1")
p2 <- ggplot(mtcars, aes(gear, fill = as.factor(gear))) + geom_bar() + ggtitle("Plot 2")
p1 + p2
p1 / p2
plot_grid(p1, p2, labels = c("A", "B"), ncol = 1, align = "v")



##################      EAS      #####################################
#######################################################################
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
})

# Set paths (keep your actual paths unchanged)
setwd("/gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/figures/fig1")
meta_path <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/01mixbatches_group/meta.txt2"
tpm_path  <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/01mixbatches_group/subset_tpm.txt"

# 1. Read data
meta <- read.table(meta_path, header = TRUE, stringsAsFactors = FALSE)
tpm <- read.table(tpm_path, header = TRUE, row.names = 1, check.names = FALSE)

# 2. Extract target genes (C4A / C4B)
target_genes <- c("C4A" = "ENSG00000244731.9", 
                  "C4B" = "ENSG00000224389.9")

valid_genes <- target_genes[target_genes %in% rownames(tpm)]
if(length(valid_genes) == 0) stop("Error: Target gene IDs not found in TPM matrix.")

tpm_t <- as.data.frame(t(tpm[valid_genes, , drop = FALSE]))
colnames(tpm_t) <- names(valid_genes)
tpm_t$sample <- rownames(tpm_t)

# 3. Merge and classify by developmental stages (Infancy excluded)
df_plot <- merge(meta, tpm_t, by = "sample") %>%
  pivot_longer(cols = all_of(names(valid_genes)), names_to = "Gene", values_to = "TPM") %>%
  mutate(
    # Adjust prenatal stage thresholds based on PCW to GW conversion (GW = PCW + 2)
    Period = case_when(
      Stage == "stage1" & age < 15 ~ "Early, prenatal",                   
      Stage == "stage1" & age >= 15 & age < 26 ~ "Middle, prenatal",     
      # Fix gap: changed 28 to 26 to prevent losing samples from weeks 26-27
      Stage == "stage1" & age >= 26 ~ "Late, prenatal",                   
      # Exclude Infancy (< 1)
      Stage %in% c("stage2", "stage3") & age >= 1 & age < 12 ~ "Childhood",
      Stage %in% c("stage2", "stage3") & age >= 12 & age < 20 ~ "Adolescence",
      Stage %in% c("stage2", "stage3") & age >= 20 & age < 40 ~ "Young adulthood",
      Stage %in% c("stage2", "stage3") & age >= 40 & age < 60 ~ "Middle adulthood",
      Stage %in% c("stage2", "stage3") & age >= 60 ~ "Late adulthood"
    ),
    # Core modification: Update factor levels to ensure correct X-axis order
    Period = factor(Period, levels = c(
      "Early, prenatal", "Middle, prenatal", "Late, prenatal",
      "Childhood", "Adolescence", "Young adulthood", "Middle adulthood", "Late adulthood"
    )),
    # Generate a numeric X-axis variable for geom_smooth and continuous scaling
    Period_Num = as.numeric(Period)
  ) %>%
  filter(!is.na(Period)) # Automatically exclude Infancy and unmapped samples

# 4. Extract unique periods and their numeric mapping for the X-axis
axis_labels <- df_plot %>%
  select(Period_Num, Period) %>%
  distinct() %>%
  arrange(Period_Num)


# ---------------------------------------------------------
# 5. Create plotting function for C4A and C4B
# ---------------------------------------------------------
plot_gene_trajectory <- function(gene_name, box_color, y_label = "TPM") {
  
  # Filter data for a single gene
  df_sub <- df_plot %>% filter(Gene == gene_name)
  
  p <- ggplot(df_sub, aes(x = Period_Num, y = TPM)) + 
  
    geom_smooth(method = "loess", se = TRUE, 
                color = "#4169E1", fill = "grey85", alpha = 0.4, span = 0.8) +
    geom_boxplot(aes(group = Period_Num), 
                 width = 0.35, fill = box_color, color = "black", 
                 outlier.shape = NA) +
    # Map numeric breaks back to period names without sample sizes
    scale_x_continuous(breaks = axis_labels$Period_Num, 
                       labels = axis_labels$Period) +
    coord_cartesian(ylim = c(0, 50)) +
    theme_classic(base_size = 22) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 18, lineheight = 0.8),
      axis.text.y = element_text(color = "black", size = 20),
      axis.title = element_text(face = "bold", size = 22),
      axis.line = element_line(linewidth = 1),
      axis.ticks = element_line(linewidth = 1),
      plot.title = element_text(face = "bold", size = 24, hjust = 0.5),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 40)
    ) +
    labs(title = gene_name, x = NULL, y = y_label)
  
  return(p)
}

# 6. Generate and save C4A plot (using light red)
p_C4A <- plot_gene_trajectory("C4A", "#E69191")
ggsave("C4A_expression_boxplot_reNum_lifespantrajectory_new.png", plot = p_C4A, width = 8, height = 7, dpi = 500)
ggsave("C4A_expression_boxplot_reNum_lifespantrajectory_new.pdf", plot = p_C4A, width = 8, height = 7, dpi = 500)

# 7. Generate and save C4B plot (using light blue)
p_C4B <- plot_gene_trajectory("C4B", "#92B5CA")
ggsave("C4B_expression_boxplot_reNum_lifespantrajectory_new.png", plot = p_C4B, width = 8, height = 7, dpi = 500)
ggsave("C4B_expression_boxplot_reNum_lifespantrajectory_new.pdf", plot = p_C4B, width = 8, height = 7, dpi = 500)
