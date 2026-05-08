library(tidyverse)
library(data.table)

file_path <- "fetal_onlychb_imputed_haps_refined.R5.txt"
df <- fread(file_path, header = TRUE)
all_haps <- c(df$H1, df$H2)
all_haps <- gsub("H_", "", all_haps)
all_haps <- gsub("_", "-", all_haps)
hap_counts <- as.data.frame(table(all_haps))
colnames(hap_counts) <- c("Haplotype", "Count")

hap_counts$Freq <- hap_counts$Count / sum(hap_counts$Count)
hap_counts <- hap_counts %>% arrange(desc(Freq))

top_n <- 8

if(nrow(hap_counts) > top_n) {
  top_haps <- hap_counts[1:top_n, ]

  others_count <- sum(hap_counts[(top_n + 1):nrow(hap_counts), "Count"])
  others_freq <- sum(hap_counts[(top_n + 1):nrow(hap_counts), "Freq"])
  others_row <- data.frame(Haplotype = "Others", Count = others_count, Freq = others_freq)
  plot_df <- bind_rows(top_haps, others_row)
} else {
  plot_df <- hap_counts
}

plot_df <- plot_df %>%
  mutate(
    Label = case_when(
      Haplotype == "2-1-1-1" ~ "2-1-1-1\n(AL-BS)",
      Haplotype == "2-1-1-2" ~ "2-1-1-2\n(AL-BL)",
      TRUE ~ as.character(Haplotype) 
    ),
    Label = factor(Label, levels = c(Label[Label != "Others"], "Others"))
  )

p <- ggplot(plot_df, aes(x = Label, y = Freq)) +

  geom_col(fill = "#9ecae1", color = "black", width = 0.55) + 

  geom_text(aes(label = sprintf("%.1f%%", Freq * 100)), 
            vjust = -0.8, size = 4, color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(color = "black", size = 11, lineheight = 0.9), 
    axis.text.y = element_text(color = "black", size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16, face = "plain"), 
    plot.margin = margin(t = 20, r = 20, b = 10, l = 10)
  ) +
  labs(
    title = "C4 Haplotype Frequencies (Prenatal Brain)", 
    x = "Haplotype Structure",
    y = "Frequency"
  )

print(p)

ggsave("C4_Haplotype_Frequencies_PrenatalBrain.png", p, width = 6, height = 5, dpi = 400)



library(tidyverse)
library(data.table)

file_path <- "adult_only_imputed_haps_refinedR5.txt"
df <- fread(file_path, header = TRUE)
all_haps <- c(df$H1, df$H2)
all_haps <- gsub("H_", "", all_haps)
all_haps <- gsub("_", "-", all_haps)
hap_counts <- as.data.frame(table(all_haps))
colnames(hap_counts) <- c("Haplotype", "Count")

hap_counts$Freq <- hap_counts$Count / sum(hap_counts$Count)
hap_counts <- hap_counts %>% arrange(desc(Freq))

top_n <- 10

if(nrow(hap_counts) > top_n) {

  top_haps <- hap_counts[1:top_n, ]
  others_count <- sum(hap_counts[(top_n + 1):nrow(hap_counts), "Count"])
  others_freq <- sum(hap_counts[(top_n + 1):nrow(hap_counts), "Freq"])
  
  others_row <- data.frame(Haplotype = "Others", Count = others_count, Freq = others_freq)
  plot_df <- bind_rows(top_haps, others_row)
} else {
  plot_df <- hap_counts
}

plot_df <- plot_df %>%
  mutate(
    Label = case_when(
      Haplotype == "2-1-1-1" ~ "2-1-1-1\n(AL-BS)",
      Haplotype == "2-1-1-2" ~ "2-1-1-2\n(AL-BL)",
      TRUE ~ as.character(Haplotype) 
    ),
    Label = factor(Label, levels = c(Label[Label != "Others"], "Others"))
  )

p <- ggplot(plot_df, aes(x = Label, y = Freq)) +
  geom_col(fill = "#9ecae1", color = "black", width = 0.55) + 
  geom_text(aes(label = sprintf("%.1f%%", Freq * 100)), 
            vjust = -0.8, size = 4, color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(color = "black", size = 11, lineheight = 0.9),
    axis.text.y = element_text(color = "black", size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16, face = "plain"),
    plot.margin = margin(t = 20, r = 20, b = 10, l = 10)
  ) +
  labs(
    title = "C4 Haplotype Frequencies (Postnatal Brain)", 
    x = "Haplotype Structure",
    y = "Frequency"
  )

print(p)

ggsave("C4_Haplotype_Frequencies_PostnatalBrain.png", p, width = 6.5, height = 5, dpi = 400)
