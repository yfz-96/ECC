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
outpath <- paste("../../Results/Figure_8/generalized_fold_change/", sep="")
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

log.mat<-function(mat, base=2){
  if(any(mat == 0)) {
    if(sum(mat == 0) == length(mat)) {
      mat <- mat + 1e-5
    }
    v <- as.vector(mat)
    minval <- min(v[v > 0])/2
    mat <- mat + minval
  }
  out<-log(mat, base)
  return(out)
}

new_quantile <- function(x) {
  quantile(x, probs = seq(0, 1, 0.1))
}
generalized_logfc <- function(x, y, base = 2, positive_class = positive_class){
  if(nlevels(y) > 2) levels(y)[levels(y) != positive_class] <- "Others"
  QuantileAbd <- apply(x, 2, function(x) tapply(x, y, new_quantile))
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

HHCC_metadata <- subset(metadata, Future_Status_Tooth %in% c("ConfidentH", "Caries"))
HHCC_df <- subset(df, rownames(df) %in% rownames(HHCC_metadata))
HHCC_df <- HHCC_df[, which(colSums(HHCC_df) != 0)]
glogfc <- generalized_logfc(HHCC_df, factor(HHCC_metadata$Future_Status_Tooth), positive_class = "Caries")
glogfc <- data.frame(glogfc)
write.table(glogfc, paste0(outpath, "/HHCC_generalized_log_fold_change.txt"), sep = "\t", quote = F, row.names = T, col.names = NA)

status_related_features <- read.table("../../Results/Figure_6/MiC/taxonomy_Position_crossRF_out/BetweenGroupTest_out_all.xls", 
                                      sep = "\t", header = T)
status_related_features <- subset(status_related_features, IfSig == "Sig")
status_related_features_glogfc <- subset(glogfc, rownames(glogfc) %in% status_related_features$feature)
