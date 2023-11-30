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
#-------------------------------
## input args
#-------------------------------
rf.opts<-list(outdir=NA, ntree=5000, verbose=FALSE, nfolds=10)
#source("ranger_util_20190520.R")
source("../../data_trimming_util.R")
#-------------------------------
datafile <- "../../../data/1867/function/ML_features_table/1867_function_feature_table.biom"
sample_metadata <- "../../../data/1867/function/1867_function_metadata.tsv"
feature_metadata <- "../../../data/1867/function/ML_features_table/1867_function_feature_taxon.txt"
prefix_name<-"function_HHRH_models"
s_category<-"Position2" 
c_category<-"Future_Status_Tooth"
Addfunction=TRUE
p_cutoff=0.05
q_cutoff=0.05
p.adj.method="BH"
h=10
outpath <- paste("../../../Results/Figure_6/MiC/",prefix_name,"_Position_crossRF_out/", sep="")
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
#'-------------------------------
#' Train data: filtering
#'-------------------------------
data_list<-filter_samples_by_sample_ids_in_metadata(df, metadata)
#'-------------------------------
#' Train data: filter out samples with null values in target_field
#'-------------------------------
data_list<-filter_samples_by_NA_in_target_field_of_metadata(data_list$data, data_list$metadata, target_field = s_category)
#'-------------------------------
#' Train data: filter out samples in particular groups
#'-------------------------------
data_HHRH_list<-filter_samples_by_groups_in_target_field_of_metadata(data_list$data, data_list$metadata,
                                                                    target_field = c_category, negate=FALSE,
                                                                    groups = c("ConfidentH", "RelativeH"))
#-------------------------------
# rf_clf.by_datasets
#-------------------------------
## "rf_clf.by_datasets" runs standard random forests with oob estimation for classification of 
## c_category in each the sub-datasets splited by the s_category. 
## The output includes a summary of rf models in the sub datasets
## and all important statistics for each of features.

res_file<-paste(outpath, prefix_name, "_rf_clf.by_datasets_res.RData", sep="")
if(file.exists(res_file)){
  load(res_file)
}else{
  rf_clf_res<-rf_clf.by_datasets(data_HHRH_list$data, data_HHRH_list$metadata, s_category, c_category, 
                                 positive_class="RelativeH", clr_transform = TRUE, rf_imp_pvalues = FALSE,
                                 p.adj.method=p.adj.method, q_cutoff=q_cutoff, nfolds=5, verbose=FALSE, ntree=500)
  save(rf_clf_res, file=res_file)
}
rf_clf_res$feature_imps_list<-lapply(rf_clf_res$feature_imps_list, function(x) {x[x$mean_all==0, "rf_imps"]<-NA; return(x)})

rf_clf_res.summ<-plot_clf_res_list(rf_clf_res, p_cutoff=0.05, p.adj.method = p.adj.method, q_cutoff=0.05, outdir=outpath)

wilcox_res<-rf_clf_res.summ$feature_res
rf_models<-rf_clf_res$rf_model_list
# Add feature annotations using feature metadata
rbind.na<-function(l){
  max_len<-max(unlist(lapply(l, length)))
  c_l<-lapply(l, function(x) {c(x, rep(NA, max_len - length(x)))})
  do.call(rbind, c_l)
}
expand_Taxon<-function(df, Taxon){
  taxa_df <- rbind.na(strsplit(as.character(df[, Taxon]), ';'))
  colnames(taxa_df) <- c("phylum","class","order","family","genus","species") #"kingdom", 
  data.frame(df, taxa_df)
}

if(Addfunction==TRUE){
  wilcox_res<-add_ann(wilcox_res, fmetadata)
  #wilcox_res<-expand_Taxon(wilcox_res, "Taxon")
}else{
  wilcox_res$Taxon=wilcox_res$feature
}
## modify the rf imp scores
#wilcox_res[wilcox_res$mean_all==0, "rf_imps"]<-NA
library("dplyr")
wilcox_res<-wilcox_res %>% dplyr::group_by(dataset) %>% mutate(rf_imps_rank=rank(-rf_imps, na.last = "keep"))

## RF classification performance VS number of features used
# Prediction performances at increasing number of features obtained by 
# retraining the random forest regressor on the top-ranking features identified 
# with a first random forest model training in a cross-validation setting
top_n_perf_list<-list()
for(n in 1:length(rf_clf_res$rf_model_list)){
  rf_imps<-rf_clf_res$feature_imps_list[[n]][, "rf_imps"]
  rf_imps_rank<-rank(-rf_imps, na.last = "keep")
  x<-rf_clf_res$x_list[[n]]
  max_n<-ncol(x)
  n_features<-c(2, 4, 8, 16, 32, 64, 128, 256, 512, 1024)
  n_features<-n_features[n_features < max_n]
  nrow<-length(n_features)+1
  top_n_perf<-matrix(NA, ncol=2, nrow=nrow)
  colnames(top_n_perf)<-c("n_features", "AUROC")
  rownames(top_n_perf)<-top_n_perf[,1]<-c(n_features, max_n)
  for(i in 1:length(n_features)){
    idx<-which(rf_imps_rank<=n_features[i])
    x_n<-x[, idx]
    y_n<-rf_clf_res$y_list[[n]]
    top_n_rf<-rf.out.of.bag(x_n, y_n, ntree=5000)
    rf_AUROC<-get.auroc(top_n_rf$probabilities[, "RelativeH"], y_n, positive_class="RelativeH")
    top_n_perf[i, 1]<-n_features[i]
    top_n_perf[i, 2]<-rf_AUROC
  }
  top_n_perf[nrow, ]<-c(max_n, rf_clf_res$rf_AUROC[n])
  top_n_perf_list[[n]]<-top_n_perf
}
names(top_n_perf_list)<-rf_clf_res$datasets
top_n_perf_list<-lapply(1:length(top_n_perf_list), 
                        function(x) data.frame(Dataset=rep(names(top_n_perf_list)[x], nrow(top_n_perf_list[[x]])), top_n_perf_list[[x]]))

top_n_perf_comb<-do.call(rbind, top_n_perf_list)
top_n_perf_comb$n_features<-as.numeric(as.character(top_n_perf_comb$n_features))
top_n_perf_comb_m<-melt(top_n_perf_comb, id.vars = c("n_features", "Dataset"))
breaks<-top_n_perf_comb_m$n_features
p<-ggplot(subset(top_n_perf_comb_m, variable=="AUROC"), aes(x=n_features, y=value)) + 
  xlab("# of features used")+
  ylab("AUROC")+
  scale_x_continuous(trans = "log",breaks=breaks)+
  ylim(c(0, 1)) +
  geom_point() + geom_line() +facet_wrap(~Dataset) +
  theme_bw()+
  theme(axis.line = element_line(color="black"),
        axis.title = element_text(size=18),
        strip.background = element_rect(colour = "white"), 
        panel.border = element_blank(),
        axis.text.x=element_text(angle=90, hjust=1))
ggsave(filename=paste(outpath,"AUROC__top_rankings.facets.scatterplot.pdf",sep=""), plot=p, width=7, height=6)
p<-ggplot(subset(top_n_perf_comb_m, variable=="AUROC"), aes(x=n_features, y=value)) + 
  xlab("# of features used")+
  ylab("AUROC")+
  scale_x_continuous(trans = "log",breaks=breaks)+
  ylim(c(0.5, 1))+
  geom_point(aes(color=Dataset)) + geom_line(aes(color=Dataset)) +#facet_wrap(~Dataset) +
  theme_bw()+
  theme(axis.line = element_line(color="black"),
        axis.title = element_text(size=18),
        strip.background = element_rect(colour = "white"), 
        panel.border = element_blank(),
        axis.text.x=element_text(angle=90, hjust=1))
ggsave(filename=paste(outpath,"AUROC__top_rankings.scatterplot.pdf",sep=""), plot=p, width=5, height=4)
sink(paste(outpath,"top_n_feature_perf.xls",sep=""));write.table(top_n_perf_comb,quote=FALSE,sep="\t", row.names = F);sink()



## Non-specific features across datasets (at least present in two of datasets)
res_spcf<-id_non_spcf_markers(wilcox_res, positive_class="RelativeH", other_class="ConfidentH", p.adj.method, outdir=outpath)
#summary(res_spcf)
wilcox_res_spcf<-res_spcf$feature_res_spcf


## keeps only OTUs which were significant in at least one datasets.
keep_sig_markers<-function(wilcox_res, q_cutoff){
  wilcox_res_sig<-wilcox_res[which(wilcox_res[, "feature"] %in% unique(as.character(subset(wilcox_res, non.param.test_p.adj < q_cutoff)[, "feature"]))), ]
  wilcox_res_sig
}

wilcox_res_spcf_sig<-keep_sig_markers(wilcox_res_spcf, q_cutoff=q_cutoff)
sink(paste(outpath,"BetweenGroupTest_out_all.xls",sep=""));write.table(wilcox_res_spcf,quote=FALSE,sep="\t", row.names = F);sink(NULL)
sink(paste(outpath,"BetweenGroupTest_out_sig_",p.adj.method,"_q_cutoff_",q_cutoff,".xls",sep=""));write.table(wilcox_res_spcf_sig,quote=FALSE,sep="\t", row.names = F);sink(NULL)
#-------------------------------
# Plotting 
#-------------------------------
## customised colors for Enr's 3 factors: XX-depleted, XX-enriched and Neutral
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
my3cols<-c("grey60", rev(gg_color_hue(2)))
#-------------------------------
# The performance of rf models
#-------------------------------
m_s<-split(data_HHRH_list$metadata, data_HHRH_list$metadata[, s_category])
rocobj_list<-list()
ciobj_list<-list()
for(i in 1:length(rf_models)){
  sub_group<-m_s[[i]][, c_category]
  #--------------------------------------------------
  #  ROC plot using "pROC" package
  #--------------------------------------------------
  pdf(paste(outpath, "rf_", names(m_s)[i],".pROC.ci.pdf",sep=""), width=4, height=4)
  rocobj_list[[i]]<-rocobj <- plot.roc(sub_group, rf_models[[i]]$probabilities[,2],main="", percent=TRUE,ci=TRUE) # print the AUROC (will contain the CI)
  ciobj_list[[i]]<-ciobj <- ci.se(rocobj,specificities=seq(0, 100, 5)) # over a select set of specificities
  plot(ciobj, type="shape", col="#1c61b6AA", print.auc=TRUE) # plot as a blue shape
  text(70,25, paste0(levels(sub_group), collapse = " VS "), pos=4)
  text(70,15, paste("AUROC = ",formatC(rocobj$auc,digits=2,format="f"),sep=""),pos=4)
  ci.lower<-formatC(rocobj$ci[1],digits=2,format="f")
  ci.upper<-formatC(rocobj$ci[3],digits=2,format="f")
  text(70,5, paste("95% CI: ",ci.lower,"-",ci.upper,sep=""),pos=4)
  dev.off()
}

#-------------------------------
# AUROC distribution
#-------------------------------
p<-ggplot(wilcox_res_spcf_sig, aes(AUROC, fill=dataset)) + geom_histogram(alpha=0.4) + 
  #scale_fill_manual(values = my3cols)+
  facet_wrap(~Enr, nrow=1) 
ggsave(filename=paste(outpath,"Markers_AUROC_distr_",p.adj.method, ".",c_category, "_by_", s_category,".pdf",sep=""),plot=p, width=8, height=3)
#-------------------------------
# Mean VS sd of markers in each of datasets
#-------------------------------
l_m<-wilcox_res_spcf_sig[, c("feature", "Taxon","dataset","mean_all", "sd_all", "IfSig", "IfSigEnr", "Enr")]
p<-ggplot(l_m, aes(x=log(mean_all), y=log(sd_all), color=Enr)) + scale_color_manual(values = my3cols) +
  geom_point(stat="identity", alpha=.5) + xlab("log10(mean)")+ ylab("log10(sd)")+
  geom_abline(slope = 1, intercept = 0)+
  coord_flip()+ # if want to filp coordinate 
  theme_bw()+
  facet_wrap(~ dataset, nrow=1)
ggsave(filename=paste(outpath,"Markers_mean_VS_sd_",p.adj.method,"_",c_category, "_by_", s_category,".pdf",sep=""),plot=p, width=8, height=4)

#-------------------------------
# Statistics summary of features in all datasets
#-------------------------------
rf_clf_res.summ_all<-plot_clf_res_list(rf_clf_res, p_cutoff=1, p.adj.method = p.adj.method, q_cutoff=1, outdir=outpath)
all_wilcox_res<-rf_clf_res.summ_all$feature_res
if(Addfunction==TRUE){
  all_wilcox_res<-add_ann(all_wilcox_res, fmetadata)
  #wilcox_res<-expand_Taxon(wilcox_res, "Taxon")
}else{
  all_wilcox_res$Taxon=all_wilcox_res$feature
}
#-------------------------------
# RF imps
#-------------------------------
if(any(grepl("package:plyr", search()))) detach("package:plyr") else message("plyr not loaded")
## plyr not loaded
library("dplyr")
all_wilcox_res$Taxon<-as.factor(all_wilcox_res$Taxon)
if(any(nchar(all_wilcox_res$feature)>20)){
  all_wilcox_res<-all_wilcox_res %>% dplyr::group_by(dataset) %>% mutate(asv_ids=paste("ASV_", seq(1, length(feature)),sep = ""))
  all_wilcox_res<-all_wilcox_res %>% dplyr::group_by(dataset) %>% mutate(rf_imps_rank=rank(-rf_imps, na.last = "keep"), 
                                                                         feature_Taxon=paste(asv_ids, Taxon, sep="__"))
}else{
  all_wilcox_res<-all_wilcox_res %>% dplyr::group_by(dataset) %>% mutate(rf_imps_rank=rank(-rf_imps, na.last = "keep"), 
                                                                         feature_Taxon=paste(feature, Taxon, sep="__"))
}
                                          
all_Taxon_df<- all_wilcox_res %>% 
  dplyr::group_by(Taxon, feature) %>% 
  dplyr::summarise(n=length(rf_imps), n_NA=sum(is.na(rf_imps)), n_neg=sum(rf_imps<0, na.rm = T), n_pos=sum(rf_imps>0, na.rm = T), 
                   imp_sum = sum(rf_imps, na.rm = T), imp_mean=mean(rf_imps, na.rm = T), 
                   imp_rank_sum= sum(rf_imps_rank, na.rm = T), imp_rank_mean=mean(rf_imps_rank, na.rm = T),
                   imp_rank_min= min(rf_imps_rank, na.rm = T) #, 
                   #imp_rank_min_Enr=ifelse(length(rf_imp_rank)>1 && !all(is.na(rf_imps_rank)), Enr[which.min(rf_imp_rank)], as.character(Enr))
                   #imp_rank_min_mean_logfc=mean_logfc[which.min(rf_imp_rank)]
                   )

write.table(all_Taxon_df, paste(outpath, "all_Taxon_OTU.RF_imps.tsv", sep=""), sep="\t", quote=FALSE, row.names = FALSE)

