# 4. Plot
library(tidyverse)
library(ggplot2)

magma_file <- "SCZ_Heritability_Results.gsa.out"
magma_res <- read.table(magma_file, header = TRUE, stringsAsFactors = FALSE)


magma_res <- magma_res %>%
  mutate(
    Network = factor(VARIABLE, levels = c("Prenatal_Positive", "Prenatal_Negative", "Postnatal_Positive", "Postnatal_Negative")),
    Stage = ifelse(grepl("Prenatal", VARIABLE), "Prenatal", "Postnatal"),
    Stage = factor(Stage, levels = c("Prenatal", "Postnatal")),
    Sig_Label = case_when(
      P < 0.001 ~ "***",
      P < 0.01 ~ "**",
      P < 0.05 ~ "*",
      TRUE ~ ""
    ),

    Label_Y = BETA + SE + 0.004 
  )

p_magma <- ggplot(magma_res, aes(x = Network, y = BETA, fill = Stage)) +
  geom_col(color = "black", alpha = 0.85, width = 0.6) +

  geom_errorbar(aes(ymin = BETA - SE, ymax = BETA + SE), width = 0.15, color = "black", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 1) +
  geom_text(aes(y = Label_Y, label = Sig_Label), size = 8, fontface = "bold", color = "black") +
  scale_fill_manual(values = c("Prenatal" = "#FF9999", "Postnatal" = "#6BAED6")) +
  
  theme_classic(base_size = 16) +
  labs(
    title = "SCZ Heritability Enrichment in C4A Seed Networks",
    x = "Seed Network",
    y = "Effect Size (beta ± SE)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 17),
    axis.text.x = element_text(angle = 30, hjust = 1, color = "black", size = 16),
    axis.text.y = element_text(color = "black"),
    legend.position = "top",
    legend.title = element_blank()
  )

print(p_magma)

ggsave("MAGMA_SCZ_Heritability_R0.3.pdf", p_magma, width = 6.5, height = 5.5)
ggsave("MAGMA_SCZ_Heritability_R0.3.png", p_magma, width = 6.5, height = 5, dpi = 400)

