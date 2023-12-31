# install and load necessary libraries for data analyses
#-------------------------------
p <- c("ggplot2", "RColorBrewer", "viridis", "palmerpenguins", "tidyverse", "pheatmap", 
       "colorspace", "plyr", "stats", "grDevices", "graphics", "biomformat", "reshape2",
       "corrplot", "GGally", "ggimage", "pdftools", "ggpubr", "DistVis")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="https://cloud.r-project.org/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))
source("../data_trimming_util.R")

outpath <- "../../Results/Figure_3/"
pc1_correlation_boxplot <- function(table, metadata, pre_ourdir) {
  identical(rownames(table), rownames(metadata))
  
  table$Var1 <- rownames(table)
  table$Var2 <- "PC1"
  table$value <- table$PC1
  table <- table[, -c(1, 2, 3)]
  
  table[, c("HostID", "Timepoint", "Position2", "HostGroup")] <- 
    metadata[table$Var1, c("HostID", "Timepoint", "Position2", "HostGroup")]
  
  table <- acast(table, HostID + Timepoint + Var2 ~ Position2)
  cor_res <- as.data.frame(cor(table, method = "spearman"))
  print(cor_res)
  
  mapping <- read.table("mapping.txt", sep = "\t", header = T, row.names = 1)
  cor_res$x <- rownames(cor_res)
  df <- melt(cor_res)
  colnames(df) <- c("Var1", "Var2", "value")
  df <- merge(df, mapping, by.x = 1, by.y = 0, all = F)
  df <- merge(df, mapping, by.x = 2, by.y = 0, all = F)
  
  
  df <- subset(df, Var1 != Var2)
  
  index_adj <- with(df, which(UL.x == UL.y & abs(UL_num.x - UL_num.y) == 1))
  index_lr <- with(df, which(LR.x != LR.y & LR_num.x == LR_num.y))
  index_ul <- with(df, which(UL.x != UL.y & UL_num.x == UL_num.y))
  index_all <- 1:nrow(df)
  index_other <- which(!(index_all %in% c(index_adj, index_lr, index_ul)))
  
  adj <- df[index_adj, ]
  adj$group <- "Adjacent"
  lr <- df[index_lr, ]
  lr$group <- "Left_VS_Right"
  ul <- df[index_ul, ]
  ul$group <- "Upper_VS_Lower"
  other <- df[index_other, ]
  other$group <- "Others"
  
  df <- rbind(adj, lr, ul, other)
  
  kruskal_test_result <-kruskal.test(data = df, value ~ group)
  print(kruskal_test_result)
  
  p <- ggboxplot(df, x = "group", y = "value",
                 color = "group", 
                 add = "jitter", shape = 1) +
    scale_colour_manual(values = viridis(4)) +
    # geom_signif(comparisons = list(c("Adjacent", "Left_VS_Right"), c("Left_VS_Right", "Upper_VS_Lower"),
    #                                c("Upper_VS_Lower", "Others")), test = "wilcox.test", textsize = 2,
    #             test.args = list(exact = FALSE, correct = FALSE, conf.int = TRUE, conf.level = 0.95)) +
    annotate("text", x = 2.5, y = 0.1, 
             label=paste("p-value = ", round(kruskal_test_result$p.value, 5)), 
             hjust = 0.5, size = 5) +
    # scale_colour_manual(values = brewer.pal(8, "Pastel2")) +
    xlab("Spatial group") + 
    ylab("PC1 correlation between\n two teeth positions") + 
    ylim(0, 1) + 
    theme_bw() + 
    theme(strip.text = element_text(size = 14),
          legend.position = "none",
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(size = 14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())

  ggsave(filename = paste0(pre_ourdir, "_correlation_boxplot.pdf"), plot = p, width = 3, height = 4)
}

metadata_file="../../data/637/taxonomy/637_taxonomy_metadata.tsv"
data_file="../../data/637/taxonomy/637_all_asv_phylo_rpca/phylo_rpca.biplot/PCs.txt"
table <- read.table(data_file, header = T, row.names = 1)
table <- table[order(rownames(table)), ]

metadata<-read.table(metadata_file,header=T,sep="\t",row.names=1, quote="", comment.char="")
metadata<-metadata[which(rownames(metadata) %in% rownames(table)), ]
metadata<-metadata[order(rownames(metadata)), ]
identical(rownames(table), rownames(metadata))
stat_metadata <- table(metadata[, c("HostID", "Timepoint")])
stat_metadata <- melt(stat_metadata)
stat_metadata <- subset(stat_metadata, value == 20)
stat_metadata$HostID_Timepoint <- paste(stat_metadata$HostID, stat_metadata$Timepoint, sep = "_")
print(nrow(stat_metadata))
metadata <- subset(metadata, HostID_Timepoint %in% stat_metadata$HostID_Timepoint)
table <- subset(table, rownames(table) %in% rownames(metadata))
identical(rownames(table), rownames(metadata))

metadata_h <- subset(metadata, HostGroup == "H2H")
table_h <- subset(table, rownames(table) %in% rownames(metadata_h))
pc1_correlation_boxplot(table_h, metadata_h, paste0(outpath, "/Fig3E1_637_taxonomy_H2H"))

metadata_h <- subset(metadata, HostGroup != "H2H")
table_h <- subset(table, rownames(table) %in% rownames(metadata_h))
pc1_correlation_boxplot(table_h, metadata_h, paste0(outpath, "/Fig3E2_637_taxonomy_C2C"))


metadata_file="../../data/637/function/637_function_metadata.tsv"
data_file="../../data/637/function/637_all_ko_rpca/rpca.biplot/PCs.txt"
table <- read.table(data_file, header = T, row.names = 1)
table <- table[order(rownames(table)), ]

metadata<-read.table(metadata_file,header=T,sep="\t",row.names=1, quote="", comment.char="")
metadata<-metadata[which(rownames(metadata) %in% rownames(table)), ]
metadata<-metadata[order(rownames(metadata)), ]
identical(rownames(table), rownames(metadata))
stat_metadata <- table(metadata[, c("HostID", "Timepoint")])
stat_metadata <- melt(stat_metadata)
stat_metadata <- subset(stat_metadata, value == 20)
stat_metadata$HostID_Timepoint <- paste(stat_metadata$HostID, stat_metadata$Timepoint, sep = "_")
print(nrow(stat_metadata))
metadata <- subset(metadata, HostID_Timepoint %in% stat_metadata$HostID_Timepoint)
table <- subset(table, rownames(table) %in% rownames(metadata))
identical(rownames(table), rownames(metadata))

metadata_h <- subset(metadata, HostGroup == "H2H")
table_h <- subset(table, rownames(table) %in% rownames(metadata_h))
pc1_correlation_boxplot(table_h, metadata_h, paste0(outpath, "/Fig3K1_637_function_H2H"))

metadata_h <- subset(metadata, HostGroup != "H2H")
table_h <- subset(table, rownames(table) %in% rownames(metadata_h))
pc1_correlation_boxplot(table_h, metadata_h, paste0(outpath, "/Fig3K2_637_function_C2C"))


metadata_file="../../data/1867/taxonomy/1867_taxonomy_metadata.tsv"
options(scipen = 4)
data_file="../../data/1867/taxonomy/1867_all_asv_phylo_rpca/phylo_rpca.biplot/PCs.txt"
table <- read.table(data_file, header = T, row.names = 1)
table <- table[order(rownames(table)), ]
metadata<-read.table(metadata_file,header=T,sep="\t",row.names=1, quote="", comment.char="")
metadata<-metadata[which(rownames(metadata) %in% rownames(table)), ]
identical(rownames(table), rownames(metadata))

stat_metadata <- table(metadata[, c("HostID", "Timepoint")])
stat_metadata <- melt(stat_metadata)
stat_metadata <- subset(stat_metadata, value == 9)
stat_metadata$HostID_Timepoint <- paste(stat_metadata$HostID, stat_metadata$Timepoint, sep = "_")
print(nrow(stat_metadata))
metadata$HostID_Timepoint <- paste(metadata$HostID, metadata$Timepoint, sep = "_")
metadata <- subset(metadata, HostID_Timepoint %in% stat_metadata$HostID_Timepoint)
table <- subset(table, rownames(table) %in% rownames(metadata))
identical(rownames(table), rownames(metadata))

metadata_h <- subset(metadata, HostGroup == "H2H")
table_h <- subset(table, rownames(table) %in% rownames(metadata_h))
pc1_correlation_boxplot(table_h, metadata_h, paste0(outpath, "/Fig3F1_1867_taxonomy_H2H"))

metadata_h <- subset(metadata, HostGroup != "H2H")
table_h <- subset(table, rownames(table) %in% rownames(metadata_h))
pc1_correlation_boxplot(table_h, metadata_h, paste0(outpath, "/Fig3F2_1867_taxonomy_H2C_C2C"))



metadata_file="../../data/1867/function/1867_function_metadata.tsv"
options(scipen = 4)
data_file="../../data/1867/function/1867_all_ko_rpca/rpca.biplot/PCs.txt"
table <- read.table(data_file, header = T, row.names = 1)
metadata<-read.table(metadata_file,header=T,sep="\t",row.names=1, quote="", comment.char="")
metadata<-metadata[which(rownames(metadata) %in% rownames(table)), ]
identical(rownames(table), rownames(metadata))

stat_metadata <- table(metadata[, c("HostID", "Timepoint")])
stat_metadata <- melt(stat_metadata)
stat_metadata <- subset(stat_metadata, value == 9)
stat_metadata$HostID_Timepoint <- paste(stat_metadata$HostID, stat_metadata$Timepoint, sep = "_")
print(nrow(stat_metadata))
metadata$HostID_Timepoint <- paste(metadata$HostID, metadata$Timepoint, sep = "_")
metadata <- subset(metadata, HostID_Timepoint %in% stat_metadata$HostID_Timepoint)
table <- subset(table, rownames(table) %in% rownames(metadata))
identical(rownames(table), rownames(metadata))

metadata_h <- subset(metadata, HostGroup == "H2H")
table_h <- subset(table, rownames(table) %in% rownames(metadata_h))
pc1_correlation_boxplot(table_h, metadata_h, paste0(outpath, "/Fig3L1_1867_function_H2H"))

metadata_h <- subset(metadata, HostGroup != "H2H")
table_h <- subset(table, rownames(table) %in% rownames(metadata_h))
pc1_correlation_boxplot(table_h, metadata_h, paste0(outpath, "/Fig3L21867_function_H2C_C2C"))

