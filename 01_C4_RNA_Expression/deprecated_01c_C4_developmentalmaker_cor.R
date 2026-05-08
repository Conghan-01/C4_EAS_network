

library(tidyverse)
library(data.table)
library(ggpubr)


tpm_file <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/01_data/RNA-seq/stage1/raw_tpm.txt"
meta_file <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/01_data/RNA-seq/stage1/meta.filteredsamplessex.txt"

genes_to_extract <- c("C4A" = "ENSG00000244731", "C4B" = "ENSG00000224389", "SYP" = "ENSG00000102003")

message("Extracting target genes...")
tpm_raw <- fread(tpm_file, header = TRUE)

df_genes <- tpm_raw %>%
  filter(apply(., 1, function(row) any(sapply(genes_to_extract, function(g) grepl(g, row[1]))))) %>%
  column_to_rownames("exp") %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("ID")

for(i in seq_along(genes_to_extract)) {
  colnames(df_genes) <- gsub(paste0(".*", genes_to_extract[i], ".*"), names(genes_to_extract)[i], colnames(df_genes))
}


meta <- fread(meta_file, header = FALSE)
colnames(meta) <- c("ID", "sex", "RIN", "batch", "Age")

df_merged <- inner_join(df_genes, meta, by = "ID") %>%
  mutate(across(c(C4A, C4B, SYP, Age, RIN), as.numeric),
         batch = as.factor(batch),
         sex = as.factor(sex))

message("Correcting for batch effects and RIN...")

get_corrected_val <- function(gene_name, data) {
  data$y <- log2(data[[gene_name]] + 1)
  fit <- lm(y ~ batch + RIN, data = data)
  return(residuals(fit) + mean(data$y, na.rm = TRUE))
}

df_corrected <- df_merged %>%
  mutate(
    C4A_adj = get_corrected_val("C4A", .),
    C4B_adj = get_corrected_val("C4B", .),
    SYP_adj = get_corrected_val("SYP", .)
  )

message("Plotting batch-corrected correlations...")

p_corr <- ggplot(df_corrected, aes(x = SYP_adj, y = C4A_adj, color = sex)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = TRUE, size = 1.2) +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 5) +
  facet_wrap(~sex) +
  scale_color_manual(values = c("female" = "#AD5C8A", "male" = "#008B8B")) +
  theme_bw(base_size = 21) +
  labs(
    title = "C4A vs SYP expression Correlation",
    x = "Residuals SYP (log2 scale)",
    y = "Residuals C4A (log2 scale)"
  ) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

print(p_corr)
ggsave("C4A_SYP_Batch_Corrected_Correlation.png", width = 6.6, height = 5, dpi = 400)



################## Other genes##########

library(tidyverse)
library(data.table)
library(ggpubr)
library(pheatmap)

# --- 1. Configurations ---
tpm_file <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/01_data/RNA-seq/stage1/raw_tpm.txt"
meta_file <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/01_data/RNA-seq/stage1/meta.filteredsamplessex.txt"

extended_markers <- c(
  "C4A" = "ENSG00000244731", "C4B" = "ENSG00000224389",
  "SOX2" = "ENSG00000181449", "PAX6" = "ENSG00000007372", # Progenitors
  "DCX"  = "ENSG00000077279", "STMN2" = "ENSG00000104435", # Immature/Migrating
  "RBFOX3" = "ENSG00000167281", "SYP" = "ENSG00000102003", # Mature/Synapse
  "GFAP" = "ENSG00000131095"                             # Astrocytes
)

# --- 2. Data Extraction ---
message("Extracting extended marker set...")
tpm_raw <- fread(tpm_file, header = TRUE)

# Filter for IDs in the list (ignoring versions)
df_genes <- tpm_raw %>%
  filter(apply(., 1, function(row) any(sapply(extended_markers, function(g) grepl(g, row[1]))))) %>%
  column_to_rownames("exp") %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("ID")

# Standardize column names to Symbols
for(i in seq_along(extended_markers)) {
  colnames(df_genes) <- gsub(paste0(".*", extended_markers[i], ".*"), names(extended_markers)[i], colnames(df_genes))
}

# --- 3. Merging and Pre-processing ---
meta <- fread(meta_file, header = FALSE)
colnames(meta) <- c("ID", "sex", "RIN", "batch", "Age")

# Ensure all marker columns are numeric
df_merged <- inner_join(df_genes, meta, by = "ID") %>%
  mutate(across(any_of(names(extended_markers)), as.numeric),
         across(c(Age, RIN), as.numeric),
         batch = as.factor(batch),
         sex = as.factor(sex))

# --- 4. Batch & Age Correction ---
message("Correcting expression for covariates...")
get_corrected_val <- function(gene_name, data) {
  # Logic: log2 transformation -> linear model -> residuals + mean
  val_vector <- data[[gene_name]]
  y <- log2(val_vector + 1)
  fit <- lm(y ~ batch + RIN + Age, data = data)
  return(residuals(fit) + mean(y, na.rm = TRUE))
}

# Apply correction to all markers in the list
df_corrected <- df_merged
for(gene in names(extended_markers)) {
  df_corrected[[paste0(gene, "_adj")]] <- get_corrected_val(gene, df_merged)
}

marker_names <- c("SOX2", "PAX6", "DCX", "STMN2", "RBFOX3", "SYP", "GFAP","C4B")
adj_cols <- paste0(marker_names, "_adj")

df_long_total <- df_corrected %>%
  select(ID, C4A_adj, all_of(adj_cols)) %>%
  pivot_longer(cols = all_of(adj_cols), 
               names_to = "Marker", 
               values_to = "Marker_Value") %>%
  mutate(Marker = gsub("_adj", "", Marker)) 

df_long_total$Marker <- factor(df_long_total$Marker, levels = marker_names)

p_total_corr <- ggplot(df_long_total, aes(x = Marker_Value, y = C4A_adj)) +
  geom_point(alpha = 0.3, size = 1.5, color = "#2F4F4F") +
  geom_smooth(method = "lm", se = TRUE, color = "#004488", fill = "#CCE5FF", size = 1.2) +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", 
           size = 5, fontface = "bold.italic", color = "black") +
  facet_wrap(~ Marker, scales = "free_x", ncol = 4) +
  theme_bw(base_size = 18) +
  labs(
    title = "Global Correlation: C4A vs. Developmental Markers",
    subtitle = "Corrected for Batch, RIN, and Age (Residuals)",
    x = "Corrected Marker Expression (log2 TPM)",
    y = "Corrected C4A Expression (log2 TPM)"
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = "black"),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 22),
    plot.subtitle = element_text(hjust = 0.5, size = 14, color = "grey30"),
    axis.text = element_text(color = "black")
  )

print(p_total_corr)

ggsave("C4A_Global_Markers_Correlation.png", p_total_corr, width = 14, height = 8, dpi = 500)




#meta_path <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/01mixbatches_group/meta.txt2"
#tpm_path  <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/01mixbatches_group/subset_tpm.txt"

library(tidyverse)
library(data.table)
library(ggpubr)


tpm_path  <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/01mixbatches_group/subset_tpm.txt"
meta_path <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/01mixbatches_group/meta.txt2"

extended_markers <- c(
  "C4A"    = "ENSG00000244731", "C4B" = "ENSG00000224389",
  "SOX2"   = "ENSG00000181449", "PAX6" = "ENSG00000007372", 
  "DCX"    = "ENSG00000077279", "STMN2" = "ENSG00000104435", 
  "RBFOX3" = "ENSG00000167281", "SYP" = "ENSG00000102003", 
  "GFAP"   = "ENSG00000131095"
)


message("Loading and filtering for Stage 1 (Prenatal) samples...")
meta <- fread(meta_path, header = TRUE) 
meta_stage1 <- meta %>% filter(Stage == "stage1")

tpm_raw <- read.table(tpm_path,head=T,row.names=1,check.names=F)
stage1_ids <- meta_stage1$sample

tpm_raw <- as.data.frame(tpm_raw)
target_rows <- which(sapply(rownames(tpm_raw), function(x) {
  any(sapply(extended_markers, function(g) grepl(g, x)))
}))

df_genes <- tpm_raw[target_rows, ] %>%
  select(any_of(stage1_ids)) %>%
  t() %>% 
  as.data.frame() %>%
  rownames_to_column("sample")

for(i in seq_along(extended_markers)) {
  idx <- grep(extended_markers[i], colnames(df_genes))
  if(length(idx) > 0) {
    colnames(df_genes)[idx] <- names(extended_markers)[i]
  }
}
if(any(duplicated(colnames(df_genes)))) {
  colnames(df_genes) <- make.unique(colnames(df_genes))
}

df_final <- df_genes %>%
  mutate(across(-sample, ~ log2(as.numeric(.x) + 1)))

for(i in seq_along(extended_markers)) {
  colnames(df_genes) <- gsub(paste0(".*", extended_markers[i], ".*"), names(extended_markers)[i], colnames(df_genes))
}

marker_names <- c("SOX2", "PAX6", "DCX", "STMN2", "RBFOX3", "SYP", "GFAP", "C4B")

df_long <- df_final %>%
  select(sample, C4A, any_of(marker_names)) %>%
  pivot_longer(cols = any_of(marker_names), 
               names_to = "Marker", 
               values_to = "Marker_Value")

df_long$Marker <- factor(df_long$Marker, levels = intersect(marker_names, unique(df_long$Marker)))

p_corr <- ggplot(df_long, aes(x = Marker_Value, y = C4A)) +
  geom_point(alpha = 0.4, size = 1.8, color = "#2c3e50") +
  geom_smooth(method = "lm", se = TRUE, color = "#e74c3c", fill = "#fadbd8") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", 
           size = 4.5, fontface = "bold.italic") +
  facet_wrap(~ Marker, scales = "free", ncol = 4) +
  theme_bw(base_size = 15) +
  labs(
    title = "Stage 1 Correlation: C4A vs. Brain Markers",
    subtitle = "Raw Expression Analysis (log2 TPM+1) - No Batch Correction",
    x = "Marker Expression [log2(TPM+1)]",
    y = "C4A Expression [log2(TPM+1)]"
  ) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

print(p_corr)

ggsave("C4A_Stage1_Correlation_Raw.png", p_corr, width = 12, height = 7, dpi = 300)