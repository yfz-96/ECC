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

imp_score <- read.table("../../Results/Figure_6/MiC/taxonomy_HHRH_models_Position_crossRF_out/all_Taxon_OTU.RF_imps.tsv", 
                        sep = "\t", header = T)
imp_score <- imp_score[order(imp_score$imp_rank_mean), ]
top_32_features <- imp_score[1:32, ]

Taxonomy <- read.table("../../data/1867/taxonomy/silva_taxonomy_collapse/level8/taxonomy.tsv", 
                       sep = "\t", header = T, row.names = 1)

getLastElement <- function(s) {
  elements <- strsplit(s, ";")[[1]]
  return(tail(elements, 1))
}

Taxonomy$Taxon <- with(Taxonomy, paste(Taxon, ASV_ID, sep = "_"))
Taxonomy$Taxon <- sapply(Taxonomy$Taxon, getLastElement)
Taxonomy$Taxon <- trimws(Taxonomy$Taxon)

glogfc <- read.table("../../Results/Figure_8/generalized_log2_fold_change/HHRH_generalized_logfc_by_datasets.txt",
                     sep = "\t", header = T, row.names = 1)
glogfc <- subset(glogfc, rownames(glogfc) %in% top_32_features$feature)
glogfc$feature <- rownames(glogfc)
glogfc <- melt(glogfc)
glogfc$Taxon <- Taxonomy[glogfc$feature, "Taxon"]

glogfc_plot <- ggplot(glogfc, aes(variable, Taxon, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "gray") +
  #scale_fill_gradient(low="#EFEFFF", high="blue") +  
  labs(fill = "RelativeH\nVS\nConfidentH\ngeneralized\nlog2(fold change)") +
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
ggsave(filename="../../Results/Figure_8/generalized_log2_fold_change/HHRH_top32_glogfc.pdf", plot=glogfc_plot, width=10, height=10)

