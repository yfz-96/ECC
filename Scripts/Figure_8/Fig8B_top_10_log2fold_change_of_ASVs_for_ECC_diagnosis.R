#-------------------------------
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
outpath <- "../../Results/Figure_8/"

wilcox_res_result <- read.table("../../Results/Figure_6/MiC/taxonomy_Position_crossRF_out/BetweenGroupTest_out_all.xls",
                                sep = "\t", header = T)
wilcox_res_result <- wilcox_res_result[order(-wilcox_res_result$mean_logfc), ]
picked_asvs <- unique(wilcox_res_result$feature)[1:10]
wilcox_res_result <- subset(wilcox_res_result, feature %in% picked_asvs)
write.table(wilcox_res_result, paste0(outpath, "/Top_10_log2fold_change_asvs_for_ECC_diagnosis.txt"), 
            sep = "\t", quote = F, row.names = F, col.names = NA)
getLastElement <- function(s) {
  elements <- strsplit(s, ";")[[1]]
  return(tail(elements, 1))
}

mean_logfc_heatmap <- function(wilcox_res_result) {
  wilcox_res_result_importance_score <- wilcox_res_result
  
  meta_637 <- read.table("../../data/637/taxonomy/637_taxonomy_metadata.tsv", sep = "\t", header = T, row.names = 1)
  mapping <- unique(meta_637[, c("Tooth_num2", "Position2")])
  mapping <- rbind(mapping, c("T0506", "T5161"))
  colnames(mapping) = c("Tooth_num2", "position")
  
  data <- merge(wilcox_res_result_importance_score, mapping, by.x = "dataset", by.y = "position", all.x = T, all.y = F)
  taxon <- read.table("../../data/1867/taxonomy/silva_taxonomy_collapse/level8/taxonomy.tsv", 
                      sep = "\t", header = T, row.names = 1)
  taxon$Taxon <- with(taxon, paste(Taxon, ASV_ID, sep = "_"))
  taxon$Taxon <- sapply(taxon$Taxon, getLastElement)
  wilcox_res_result$Taxon <- taxon[wilcox_res_result$feature, "Taxon"]
  taxon <- unique(wilcox_res_result[, c("feature", "Taxon")])
  rownames(taxon) <- taxon$feature
  data <- data %>%
    arrange(Tooth_num2, mean_logfc)
  data$dataset <- factor(data$dataset, levels = unique(data$dataset), order = T)
  
  plot <- ggplot(data, aes(dataset, feature, fill = mean_logfc)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
    #scale_fill_gradient(low="#EFEFFF", high="blue") +  
    labs(fill = "ConfidentH\nVS\nCaries\nlog2(fold change)") +
    scale_y_discrete(labels = function(x) {taxon[x, "Taxon"]}) +
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

plot <- mean_logfc_heatmap(wilcox_res_result)
ggsave(filename=paste0(outpath, "/Fig8B_Top_10_log2fold_change_asvs_for_ECC_diagnosis.pdf"), plot=plot, width=9, height=4)

