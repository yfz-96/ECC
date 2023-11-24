## install and load necessary libraries for data analyses
#-------------------------------
p <- c("reshape2","ggplot2","pheatmap","combinat","randomForest", "gridExtra", "crossRanger", "RColorBrewer",
       "vegan", "biomformat", "doMC", "cowplot", "pROC", "ggsignif", "ggpubr", "dplyr", "viridis")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="https://cloud.r-project.org/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))
#-------------------------------
outpath = "../../Results/Figure_8/"
picked_feature <- read.table("../../data/1867/taxonomy/sig_differental_results/level_8/1867_H2H_clr_sig_feature.txt", 
                             sep = "\t", header = T, row.names = 1)

getLastElement <- function(s) {
  elements <- strsplit(s, ";")[[1]]
  return(tail(elements, 1))
}
Taxonomy <- read.table("../../data/1867/taxonomy/silva_taxonomy_collapse/level8/taxonomy.tsv", 
                    sep = "\t", header = T, row.names = 1)
Taxonomy$Taxon <- with(Taxonomy, paste(Taxon, ASV_ID, sep = "_"))
Taxonomy$Taxon <- sapply(Taxonomy$Taxon, getLastElement)
Taxonomy$Taxon <- trimws(Taxonomy$Taxon) # Remove leading and/or trailing whitespace from character strings.
 
ordered_sig_diff_asvs <- read.table("../Figure_3/ordered_spatial_related_features.txt", sep = "\t", header = T)
ordered_sig_diff_asvs <- unlist(ordered_sig_diff_asvs$x)

mean_logfc_heatmap <- function(wilcox_res_result, picked_feature) {
  wilcox_res_result_importance_score <- subset(wilcox_res_result, wilcox_res_result$feature %in% colnames(picked_feature))
  meta_637 <- read.table("../../data/637/taxonomy/637_taxonomy_metadata.tsv", sep = "\t", header = T, row.names = 1)
  mapping <- unique(meta_637[, c("Tooth_num2", "Position2")])
  mapping <- rbind(mapping, c("T0506", "T5161"))
  colnames(mapping) = c("Tooth_num2", "position")
  
  data <- merge(wilcox_res_result_importance_score, mapping, by.x = "dataset", by.y = "position", all.x = T, all.y = F)
  data$feature <-Taxonomy[data$feature, "Taxon"]
  data$feature <- factor(data$feature, levels = ordered_sig_diff_asvs, ordered = T)
  data <- data %>%
    arrange(feature, Tooth_num2)
  data$dataset <- factor(data$dataset, levels = unique(data$dataset), order = T)

  plot <- ggplot(data, aes(dataset, feature, fill = mean_logfc)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
    #scale_fill_gradient(low="#EFEFFF", high="blue") +  
    labs(fill = "ConfidentH\nVS\nCaries\nlog2(fold change)") +
    xlab("Tooth position") + 
    ylab("ASV") + 
    theme_bw() + 
    theme(strip.text = element_text(size = 14),
          legend.title = element_text(size = 14), 
          legend.text = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(size = 14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  return (plot)
}

wilcox_res_result <- read.table("../../Results/Figure_6/MiC/taxonomy_Position_crossRF_out/BetweenGroupTest_out_all.xls",
                            sep = "\t", header = T)
plot <- mean_logfc_heatmap(wilcox_res_result, picked_feature)
plot
ggsave(filename=paste0(outpath, "/Fig8A_log2fold_change_of_spatial_related_asvs_for_ECC_diagnosis.pdf"), plot=plot, width=8.5, height=3.5)

