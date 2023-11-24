# install and load necessary libraries for data analyses
#-------------------------------
p <- c("ggplot2", "RColorBrewer", "viridis", "palmerpenguins", "tidyverse", "colorspace", "plyr", "ggpubr")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="https://cloud.r-project.org/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))
source("../data_trimming_util.R")

outpath <- "../../Results/Figure_3/"
if (!dir.exists(outpath)) {
  dir.create(outpath)
} 

pc1_position_group_boxplot <- function(df, metadata, mapping = data.frame()) {
  df<-df[order(rownames(df)), ]
  metadata<-metadata[order(rownames(metadata)),]
  identical(rownames(df), rownames(metadata))
  
  df <- merge(metadata[, c("Tooth_num2", "Position2", "NicheICM", "HostGroup", "UL")], df, by="row.names", all=FALSE)
  if("T5161" %in% df$Position2) {
    x <- subset(df, Position2 == "T5161")
    x$Position2 <- "T51"
    df <- rbind(df, x)
    x$Position2 <- "T61"
    df <- rbind(df, x)
    df <- subset(df, Position2 != "T5161")
  }
  
  df <- rbind.fill(df, mapping)
  df <- df[order(df$Tooth_num2), ]
  table(df$Position2)
  df$Position2 <- factor(df$Position2, levels=unique(df$Position2), order = T)
  
  plot <- ggplot(df, aes(x = Position2, y = PC1)) +
    geom_boxplot(outlier.shape = NA) +
    scale_colour_manual(values = brewer.pal(8, "Pastel2")) +
    geom_jitter(aes(color=NicheICM), shape = 1, width = 0.2) +
    geom_smooth(mapping = aes(x = Position2, y = PC1, group = 1), data = subset(df, UL == "Upper_teeth"), method = "loess", se=TRUE) + 
    facet_wrap(~HostGroup, scales = "free_x") +
    xlab("Tooth Position") + 
    ylab("PC1") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.text = element_text(size = 14),
          legend.title = element_text(size = 14), 
          legend.text = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(size = 14),
          legend.position = "bottom")
  return(plot)
}

metadata_file="../../data/637/taxonomy/637_taxonomy_metadata.tsv"
options(scipen = 4)
data_file="../../data/637/taxonomy/637_all_asv_phylo_rpca/phylo_rpca.biplot/PCs.txt"
df<-read.table(data_file, header=T, row.names=1, quote="", comment.char = "")
metadata<-read.table(metadata_file,header=T,sep="\t",row.names=1, quote="", comment.char="")
metadata<-metadata[which(row.names(metadata) %in% row.names(df)), ]
plot <- pc1_position_group_boxplot(df, metadata)
ggsave(filename=paste(outpath, "/Fig3A_637_all_asv_phylo_rpca_PC1_position_grouphost.pdf", sep=""), plot=plot, width=10, height=3.5)

mapping <- unique(metadata[, c("Tooth_num2", "Position2", "NicheICM", "UL")])
mapping$HostGroup <- "H2H"
x <- mapping
x$HostGroup <- "H2C / C2C"
mapping <- rbind(mapping, x)
x <- mapping
x$HostGroup <- "H2C / C2C"
mapping <- rbind(mapping, x)

metadata_file="../../data/1867/taxonomy/1867_taxonomy_metadata.tsv"
options(scipen = 4)
data_file="../../data/1867/taxonomy/1867_all_asv_phylo_rpca/phylo_rpca.biplot/PCs.txt"
df<-read.table(data_file, header=T, row.names=1, quote="", comment.char = "")
metadata<-read.table(metadata_file,header=T,sep="\t",row.names=1, quote="", comment.char="")
metadata<-metadata[which(row.names(metadata) %in% row.names(df)), ]
df <- df[order(rownames(df)), ]
metadata <- metadata[order(rownames(metadata)), ]
metadata[which(metadata$HostGroup %in% c("H2C", "C2C")), "HostGroup"] = "H2C / C2C"
plot <- pc1_position_group_boxplot(df, metadata, mapping)
ggsave(filename=paste(outpath, "/Fig3B_1867_all_asv_phylo_rpca_PC1_position_grouphost.pdf", sep=""), plot=plot, width=10, height=3.5)


