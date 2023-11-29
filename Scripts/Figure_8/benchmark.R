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

log.mat<-function(mat, base=2){
  if(any(mat == 0)) {
    if(sum(mat == 0) == length(unlist(mat))) {
      mat <- mat + 1e-5
    }
    else{
      v <- as.vector(mat)
      minval <- min(v[v > 0])/2
      mat <- mat + minval
    }
  }
  out<-log(mat, base)
  return(out)
}

mean_logfc <- function(x, y, base=2, positive_class=positive_class){
  if(nlevels(y)>2) levels(y)[levels(y)!=positive_class] <- "Others"
  logMeanAbd<-log.mat(t(apply(x,2,function(x) tapply(x, y, mean))), base=base)
  out<-logMeanAbd[, positive_class]-logMeanAbd[, colnames(logMeanAbd)!=positive_class]
  out
}

median_logfc <- function(x, y, base=2, positive_class=positive_class){
  if(nlevels(y)>2) levels(y)[levels(y)!=positive_class] <- "Others"
  logMedianAbd<-log.mat(t(apply(x,2,function(x) tapply(x, y, median))), base=base)
  out<-logMedianAbd[, positive_class]-logMedianAbd[, colnames(logMedianAbd)!=positive_class]
  out
}

new_quantile <- function(x) {
  quantile(x, probs = seq(0, 1, 0.1))
}
generalized_logfc <- function(x, y, base=2, positive_class=positive_class){
  if(nlevels(y)>2) levels(y)[levels(y)!=positive_class] <- "Others"
  QuantileAbd <- apply(x,2,function(x) tapply(x, y, new_quantile))
  results <- rep(0, ncol(x))
  names(results) <- colnames(x)
  for(i in 1:length(QuantileAbd)) {
    QuantileAbd[[i]] <- do.call(data.frame, QuantileAbd[[i]])
    logQuantileAbd <- log.mat(QuantileAbd[[i]], base = base)
    out<-logQuantileAbd[, positive_class]-logQuantileAbd[, colnames(logQuantileAbd)!=positive_class]
    out<-mean(out)
    results[i] <- out
  }
  results
}

AUROC <- function(x, y) {
  rocobj <- roc(y, x, smooth = F)       # 曲线是否光滑，当光滑时，无法计算置信区间
  auc<-auc(rocobj)[1]
  auc
}

AUROCs <- function(x, y, positive_class = positive_class){
  if(nlevels(y)>2) {
    levels(y)[levels(y)!=positive_class] <- "Others"
    y <- factor(y, levels = c("Others", positive_class), ordered = T)
  }
  else {
    y1 <- unique(y[which(y != positive_class)])
    y <- factor(y, levels = c(y1, positive_class), ordered = T)
  }
  auroc <- apply(x,2,function(x) {AUROC(x, y)})
  auroc
}

outpath <- "../../Results/Figure_8/benchmark_result/"
dir.create(outpath, recursive = T)

#-------------------------------
# Biom table input
#-------------------------------
clinical_features <- c("sum_s_dt", "sum_ns_dt", "sum_dmfs", "spatial_dist_weighted_mean_dmfs")
mic_dist_features <- c("spatial_dist_weighted_mean_md", "mean_md", "mean_md_H_T51", "mean_md_H_T52", 
                       "mean_md_H_T53", "mean_md_H_T54", "mean_md_H_T55", "mean_md_H_T61", "mean_md_H_T62", 
                       "mean_md_H_T63", "mean_md_H_T64", "mean_md_H_T65", "mean_md_H_T71", "mean_md_H_T72", 
                       "mean_md_H_T73", "mean_md_H_T74", "mean_md_H_T75", "mean_md_H_T81", "mean_md_H_T82", 
                       "mean_md_H_T83", "mean_md_H_T84", "mean_md_H_T85", "mean_md_H_T5161")
if(grepl("biom$", datafile)){
  biom <- read_biom(datafile)
  df <- data.frame(t(as.matrix(biom_data(biom))), check.names = FALSE)
  df <- df[, which((!colnames(df) %in% clinical_features) & (!colnames(df) %in% mic_dist_features))]
}else{
  df<-read.table(datafile, header=T, row.names=1, sep="\t", quote="", comment.char = "")
  df <- df[, which((!colnames(df) %in% clinical_features) & (!colnames(df) %in% mic_dist_features))]
}
df<-df[order(rownames(df)), ]
#df<-sweep(df, 1, rowSums(df), "/")
#-------------------------------
# Feature metadata input
#-------------------------------
if(!is.na(feature_metadata)){
  #fmetadata<-read.table(feature_metadata,header=T, sep="\t", fill = TRUE, comment.char = "")
  fmetadata<-read.table(feature_metadata,header=T, sep='\t', quote = "",
                        row.names = NULL,
                        stringsAsFactors = FALSE, comment.char = "")
  fmetadata <- subset(fmetadata, fmetadata[, 1] %in% colnames(df))
}

add_ann<-function(tab, fmetadata, tab_id_col=1, fmetadata_id_col=1){
  fmetadata_matched<-fmetadata[which(fmetadata[, fmetadata_id_col] %in% tab[, tab_id_col]),]
  out<-merge(tab, fmetadata_matched, by.x=tab_id_col, by.y=fmetadata_id_col)
  out
}
#-------------------------------
# Sample Metadata input
#-------------------------------
allmetadata<-read.table(sample_metadata,header=T,sep="\t",row.names=1, quote="", comment.char="")
if(length(allmetadata)==1){metadata<-data.frame(allmetadata[order(rownames(allmetadata)),])
all_group<-colnames(metadata)<-colnames(allmetadata)
}else{
  metadata<-allmetadata[order(rownames(allmetadata)),]
  all_group<-colnames(metadata)
  all_group_f<-colnames(metadata)[sapply(metadata,class)=="factor"]
  all_group_n<-colnames(metadata)[sapply(metadata,class)!="factor"]
}
# heatmap of contingency tables
cols <- colorRampPalette(brewer.pal(3, "Blues"))
tab_ori<-table(metadata[, c(s_category, c_category)])
pheatmap(tab_ori, cluster_rows = F, cluster_cols = F, display_numbers = TRUE,  number_format = "%.0f", color = cols(100),
         filename = paste(outpath, "contingency_table", ".",c_category, "_by_", s_category,".pdf",sep=""), width=3, height=4)

print(identical(rownames(df), rownames(metadata)))
print(identical(colnames(df), fmetadata[, 1]))

HHCC_metadata <- subset(metadata, Future_Status_Tooth %in% c("ConfidentH", "Caries"))
HHCC_df <- subset(df, rownames(df) %in% rownames(HHCC_metadata))
HHCC_df <- HHCC_df[, which(colSums(HHCC_df) != 0)]
auroc_results <- AUROCs(HHCC_df, HHCC_metadata$Future_Status_Tooth, positive_class = "Caries")
auroc_results <- data.frame(auroc_results)
write.table(auroc_results, paste0(outpath, "/HHCC_AUROC.txt"), sep = "\t", quote = F, row.names = T, col.names = NA)

Scatter_plot <- function(x, y, x_label, y_label) {
  data <- data.frame(x, y)
  correlation <- cor.test(x, y, method = "spearman")
  ggplot(data, aes(x = x, y = y)) + 
    geom_point(color = "blue", alpha = 0.2) +
    # geom_smooth(method = lm, linetype = 1, se = FALSE, span = 1) + # 趋势线
    # stat_cor(method = "spearman",label.x = min(x), label.y = min(y)) + # library(ggpubr)
    geom_text(x = min(x), y = max(y), label = paste("rho = ", round(correlation$estimate, 3)),
              hjust = 0, vjust = 1, color = "black", size = 4) +  # 添加文本标签
    expand_limits(x = c(min(x), max(x)), y = c(min(y), max(y))) + # 设置坐标轴范围
    xlab(x_label) + # 设置X轴名称
    ylab(y_label) + # 设置Y轴名称
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
}

g_logfc <- read.table("../../Results/Figure_8/generalized_log2_fold_change/HHCC_generalized_log_fold_change.txt",
                      sep = "\t", header = T, row.names = 1)

glogfc_auroc <- merge(auroc_results, g_logfc, by = 0, all = T)
glogfc_auroc_plot <- Scatter_plot(glogfc_auroc$glogfc, glogfc_auroc$auroc_results, "Generalized log2 fold change", "AUROC")
ggsave(paste0(outpath, "/glogfc_auroc_plot.pdf"), plot = glogfc_auroc_plot, width = 5, height = 5)


mean_logfc <- read.table("../../Results/Figure_8/mean_log2_fold_change/HHCC_mean_log_fold_change.txt",
                         sep = "\t", header = T, row.names = 1)
meanlogfc_auroc <- merge(auroc_results, mean_logfc, by = 0, all = T)
meanlogfc_auroc_plot <- Scatter_plot(meanlogfc_auroc$mean_logfc, meanlogfc_auroc$auroc_results, 
                                     "Mean log2 fold change", "AUROC")
ggsave(paste0(outpath, "/mean_logfc_auroc_plot.pdf"), plot = meanlogfc_auroc_plot, width = 5, height = 5)


median_logfc <- read.table("../../Results/Figure_8/median_log2_fold_change/HHCC_median_log_fold_change.txt",
                         sep = "\t", header = T, row.names = 1)
medianlogfc_auroc <- merge(auroc_results, median_logfc, by = 0, all = T)
medianlogfc_auroc_plot <- Scatter_plot(medianlogfc_auroc$median_logfc, medianlogfc_auroc$auroc_results, 
                                     "Median log2 fold change", "AUROC")
ggsave(paste0(outpath, "/median_logfc_auroc_plot.pdf"), plot = medianlogfc_auroc_plot, width = 5, height = 5)

