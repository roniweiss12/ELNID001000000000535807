library(EPICATAC)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggrepel)


full_count_table <- read.csv("~/ELNID001000000000535807/analysis/ATAC/countTable.csv")
row.names(full_count_table) <- full_count_table$X
full_count_table$X <- NULL

rownames(full_count_table) <- gsub("_", ":", rownames(full_count_table))
rownames(full_count_table) <- sub("^([^:]+:[^:]+):", "\\1-", rownames(full_count_table))

standard_chrs <- grepl("^chr([0-9]+|X|Y):[0-9]+-[0-9]+$", rownames(full_count_table))
full_count_table <- full_count_table[standard_chrs, ]
normalized_counts <- EPICATAC:::get_TPMlike_counts(full_count_table)

# # remove neutrophils from the reference as this is PBMC data
# custom_ref <- atacRef_PBMC
# custom_ref$refProfiles <- custom_ref$refProfiles[, colnames(custom_ref$refProfiles) != "Neutrophils"]
# custom_ref$refProfiles.var <- custom_ref$refProfiles.var[, colnames(custom_ref$refProfiles.var) != "Neutrophils"]
# custom_ref$sigGenes <- custom_ref$sigGenes[custom_ref$sigGenes %in% rownames(custom_ref$refProfiles)]
# 
# out_no_neutrophils <- EPIC_ATAC(bulk = normalized_counts, reference = custom_ref, ATAC = TRUE)

out <- EPIC_ATAC(bulk = normalized_counts, reference = atacRef_PBMC, ATAC = TRUE)

cell_fractions <- as.data.frame(out[["cellFractions"]])
cell_fractions$sample <- rownames(cell_fractions)
cell_fractions$group <- ifelse(grepl("High", cell_fractions$sample), "High fatigue", "Low fatigue")

df_long <- cell_fractions %>%
  pivot_longer(cols = -c(sample, group), names_to = "cell_type", values_to = "fraction")

ggplot(df_long, aes(x = group, y = fraction, fill = group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2) +
  geom_text_repel(aes(label = sample), size = 2.5, max.overlaps = Inf) +
  facet_wrap(~ cell_type, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c("High fatigue" = "#533AB7", "Low fatigue" = "#1D9E75")) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = NULL, y = "Cell fraction", fill = NULL) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 10))

ggsave("~/ELNID001000000000535807/analysis/ATAC/celltypeProportions_noNeu.png", plot = last_plot(), width = 8, height = 6, units = "in")
