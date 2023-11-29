## install and load necessary libraries for data analyses
#-------------------------------
p <- c("reshape2","ggplot2","pheatmap","combinat","randomForest", "gridExtra", "crossRanger", "RColorBrewer",
       "vegan", "biomformat", "doMC", "cowplot", "pROC", "ggsignif", "ggpubr")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="https://cloud.r-project.org/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))

rf.opts<-list(outdir=NA, ntree=5000, verbose=FALSE, nfolds=10)
#source("ranger_util_20190520.R")
source("../data_trimming_util.R")
#-------------------------------
datafile <- "../../data/1867/taxonomy/ML_features_table/1867_taxonomy_feature_table.biom"
sample_metadata <- "../../data/1867/taxonomy/1867_taxonomy_metadata.tsv"
feature_metadata <- "../../data/1867/taxonomy/ML_features_table/1867_taxonomy_feature_taxon.txt"
prefix_name<-"taxonomy"
s_category<-"Position2" 
c_category<-"Future_Status_Tooth"
AddTaxonomy=TRUE
p_cutoff=0.05
q_cutoff=0.05
p.adj.method="BH"
h=10
outpath <- paste("../../Results/Figure_8/generalized_log2_fold_change/", sep="")
dir.create(outpath, recursive = T)

g_logfc <- read.table(paste0(outpath, "/generalized_logfc_by_datasets.txt"),
            sep = "\t", header = T, row.names = 1)
g_logfc <- apply(g_logfc, 2, function(x) {rank(-x)})
status_related_features_g_logfc <- g_logfc[which(rownames(g_logfc) %in% status_related_features$feature), ]


getLastElement <- function(s) {
  elements <- strsplit(s, ";")[[1]]
  return(tail(elements, 1))
}
Taxonomy <- read.table("../../data/1867/taxonomy/silva_taxonomy_collapse/level8/taxonomy.tsv", 
                       sep = "\t", header = T, row.names = 1)
Taxonomy$Taxon <- with(Taxonomy, paste(Taxon, ASV_ID, sep = "_"))
Taxonomy$Taxon <- sapply(Taxonomy$Taxon, getLastElement)
Taxonomy$Taxon <- trimws(Taxonomy$Taxon) # Remove leading and/or trailing whitespace from character strings.

rownames(status_related_features_g_logfc) <- Taxonomy[rownames(status_related_features_g_logfc), "Taxon"]

status_related_features_g_logfc <- melt(status_related_features_g_logfc)
status_related_features_g_logfc$Var2 <- factor(status_related_features_g_logfc$Var2,
                                               levels = c("T55", "T54", "T5161", "T64", "T65", "T75", "T74", "T84", "T85"),
                                               ordered = T)
plot <- ggplot(status_related_features_g_logfc, aes(Var2, Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "gray") +
  #scale_fill_gradient(low="#EFEFFF", high="blue") +  
  labs(fill = "ConfidentH\nVS\nCaries\nrank of\ngeneralized\nlog2(fold change)") +
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
plot
ggsave(filename=paste0(outpath, "/rank_of_generalized_log2fold_change_of_status_related_asvs_for_ECC_diagnosis.pdf"), plot=plot, width=8.5, height=3.5)
