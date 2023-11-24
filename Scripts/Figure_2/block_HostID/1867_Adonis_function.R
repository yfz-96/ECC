source("../../data_trimming_util.R")

#-------------------------------
# install and load necessary libraries for data analyses
#-------------------------------

# install.packages('devtools') # if devtools not installed
# devtools::install_github('shihuang047/crossRanger')

# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("biomformat")

p <- c("reshape2","ggplot2", "dplyr", "viridis")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="https://cloud.r-project.org/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))


sample_metadata_function_file <- "../../../data/1867/function/1867_function_metadata.tsv" 

dm_func_rpca_file <- "../../../data/1867/function/1867_all_ko_rpca/rpca.distance_matrix/distance-matrix.tsv"

outpath <- "../../../Results/Figure_2/1867_function_dist_block_hostid.PerMANOVA_out/"
dir.create(outpath)

#-------------------------------
# Sample Metadata input
#-------------------------------
sample_metadata_function <- read.table(sample_metadata_function_file,header=T,sep="\t",row.names=1, quote="", comment.char="")

#if(length(allmetadata)==1){
#  metadata<-data.frame(allmetadata[order(rownames(allmetadata)),])
#  all_group<-colnames(metadata)<-colnames(allmetadata)
#}else{
#  metadata<-allmetadata[order(rownames(allmetadata)), ]
#}

md_sample_metadata_function_summ<-check_metadata(metadata = sample_metadata_function)
write.table(md_sample_metadata_function_summ, paste0(outpath, "md_sample_metadata_function_summ.tsv"), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)


# remove two columns on sampling date 
idx = !grepl("date", md_sample_metadata_function_summ$metadata)
md_sample_metadata_function_summ <- md_sample_metadata_function_summ[idx, ]

file_list <- list(dm_func_rpca_file)
dm_list<-lapply(file_list, function(x) {
    table <- read.table(x, header=T, row.names=1, sep="\t", quote="", comment.char = "", check.names = FALSE)
    rownames(table) <- gsub("-", ".", rownames(table))
    colnames(table) <- gsub("-", ".", colnames(table))
    return(table)
    }
)

#dm_list<-lapply(dm_list, function(x) x[order(rownames(x)), order(colnames(x))])

lapply(dm_list, dim)

#for(i in 1:length(dm_list)){
#    rownames(dm_list[[i]])<-paste("X", rownames(dm_list[[i]]), sep="")
#}

#rownames(metadata)<-paste("X", rownames(metadata), sep="")

filter_samples_by_sample_ids_in_metadata <- function(data, dm=FALSE, metadata, ids_col=NA){
  data<-data[order(rownames(data)), order(colnames(data))]
  metadata<-metadata[order(rownames(metadata)),]
  if(is.na(ids_col)){
    shared_ids<-intersect(rownames(data), rownames(metadata))
    metadata_idx<-which(rownames(metadata) %in% shared_ids)
  }else{
    shared_ids<-intersect(rownames(data), metadata[, ids_col])
    metadata_idx<-which(metadata[, ids_col] %in% shared_ids)
  }
    if(dm==TRUE){
     data_matched<-data[shared_ids, shared_ids]
     data_matched<-data_matched[order(rownames(data_matched)),order(colnames(data_matched))]
     cat("The number of samples in distance matrix (after filtering out samples with no metadata): ", 
      nrow(data_matched) ,"\n")
    }else{
     data_matched<-data[shared_ids, ]
     data_matched<-data_matched[order(rownames(data_matched)),]
     cat("The number of samples in feature table (after filtering out samples with no metadata): ", 
      nrow(data_matched) ,"\n")
    }

  metadata_matched<-metadata[metadata_idx, ]
  
  cat("The number of samples metadata (after filtering out samples with no metadata): ", 
      nrow(metadata_matched) ,"\n")
  
  if(is.na(ids_col)){
    metadata_matched<-metadata_matched[order(rownames(metadata_matched)),]
    cat("The sample IDs are idenical in feature table and metadata: ", 
        identical(rownames(data_matched), rownames(metadata_matched)), "\n")
  }else{
    metadata_matched<-metadata_matched[order(as.character(metadata_matched[, ids_col])),]
    cat("The sample IDs are idenical in feature table and metadata: ", 
        identical(rownames(data_matched), as.character(metadata_matched[, ids_col])), "\n")
  }

  result<-list()
  result$data<-data_matched
  result$metadata<-metadata_matched
  return(result)
}


dm_list<-lapply(dm_list, function(x){filter_samples_by_sample_ids_in_metadata(x, dm=TRUE, sample_metadata_function)})

all_group <- c("Status_Host", "dmfs_Host", "Timepoint", "Status_Tooth",
               "dmfs_Tooth", "Position_Tooth", "HostID", "Age",
               "HostGroup", "Future_Status_Tooth")
all_group_f <- c("Status_Host", "Timepoint", "Status_Tooth",
                 "Position_Tooth", "HostID",
                 "HostGroup", "Future_Status_Tooth")
block_group <- "HostID"

stat_summ_list <- list()

for(i in 1:length(dm_list)){
    #--------------------------------
    # Statistical test: Adonis and Anosim
    #--------------------------------
    stat_summ <- data.frame(matrix(NA, nrow=length(all_group), ncol=9))
    rownames(stat_summ) <- all_group
    colnames(stat_summ) <- c("raw_sample_size", "filtered_sample_size", "num_class", "class_distribution", 
                           "Adonis.F", "Adonis.R2", "Adonis.P","Anosim.R","Anosim.P")
    #--------------------------------
    dm <- dm_list[[i]][[1]]
    metadata <- dm_list[[i]][[2]]
    
    cat(names(dm_list)[i],"\n")
    #suppressWarnings(
    for(group in all_group){
        stat_summ[group, 1] <- dim(dm)[1]
        # filter samples with NA in the metadata
        filtered_dm_list <- filter_dm_by_NA_in_target_field_of_metadata(dm, metadata, group)
        if(all(is.na(filtered_dm_list))){
            next
            stat_summ[group, 2] <- 0 
            stat_summ[group, 3:9] <- NA
        }else{
            dm_f <- filtered_dm_list$dm
            metadata_f <- filtered_dm_list$metadata
            y <- metadata_f[, group]
            stat_summ[group, 2] <- length(y)
            #--------------------------------
            if(nlevels(factor(y))==1){
                next 
                cat("All values are identical in ", group,"!\n")
            }else{
                if(is.element(group, all_group_f)){
                    y <- factor(y)
                    stat_summ[group, 3] <- nlevels(y)
                    stat_summ[group, 4] <- paste0(levels(y), collapse="|")
                    if(all(table(y)!=1)){
                        ano <- anosim(dm_f, y)
                        stat_summ[group, 8] <- ano.R <- ano$statistic
                        stat_summ[group, 9] <- ano.P <- ano$signif
                        cat("ANOSIM (", group, "): \n")
                        cat("--------------------------------")
                        print(ano)
                        ado <- adonis2(dm_f ~ y, strata = eval(parse(text = paste0("metadata_f$",block_group))))
			stat_summ[group, 5] <- ado$F[1]
                        stat_summ[group, 6] <- ado$R2[1]
                        stat_summ[group, 7] <- ado$`Pr(>F)`[1]
                        cat("ADONIS/PERMANOVA (",group,"): \n")
                        cat("--------------------------------\n")
                        print(ado)
                        cat("--------------------------------\n\n")
                    }else{
                        stat_summ[group, 5:9] <- NA
                    }
                }else{
                    y<-as.numeric(as.character(y))
                    stat_summ[group, 3]  <-  NA
                    stat_summ[group, 4]  <-  NA
                    ado <- adonis2(dm_f ~ y, strata = eval(parse(text = paste0("metadata_f$",block_group))))
                    stat_summ[group, 5] <- ado$F[1]
                    stat_summ[group, 6] <- ado$R2[1]
                    stat_summ[group, 7] <- ado$`Pr(>F)`[1]
                    cat("ADONIS/PERMANOVA (",group,"): \n")
                    cat("--------------------------------\n")
                    print(ado)
                    cat("--------------------------------\n\n")
                }
            }
        }
    stat_summ_list[[i]] <- stat_summ
    sink(paste(outpath, names(dm_list)[i], ".Beta_diversity_summ.xls",sep=""));
    cat("\t");
    write.table(stat_summ, quote=FALSE,sep='\t',row.names=TRUE);
    sink(NULL)
}
}

dist_name=c("Function")
stat_summ_list2<-lapply(1:length(stat_summ_list), function(x) data.frame(feature=rownames(stat_summ_list[[x]]), 
                                                                         dist=dist_name[x], 
                                                                         stat_summ_list[[x]]))

stat_summ_all<-data.frame(do.call("rbind", stat_summ_list2))

write.table(stat_summ_all, paste(outpath, "All_metadata", length(all_group),".Beta_diversity_summ.xls",sep=""), 
            quote=FALSE,sep='\t',row.names=FALSE)

