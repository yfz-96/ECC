rm(list = ls())

library(crossRanger)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(biomformat)
library(stringr)

source("../../data_trimming_util.R")
outdir <- "../../../Results/Figure_6/"

mic_abd_637 <- read_hdf5_biom("../../../data/637/function/KO/637_function_kos_abd.biom")
mic_abd_637 <- biom(mic_abd_637)
mic_abd_637 <- as.matrix(biom_data(mic_abd_637))
mic_abd_637 <- mic_abd_637[, order(colnames(mic_abd_637))]
clinical_features_637 <- read.table("../../../Results/Figure_S5/637_clinical_features_function.txt",
                                    sep = "\t", header = T, row.names = 1)
clinical_features_637 <- as.matrix(clinical_features_637)
clinical_features_637 <- t(clinical_features_637)
dist_features_637 <- read.table("../../../Results/Figure_S5/637_all_function_dist_features.txt", 
                                sep = "\t", header = T, row.names = 1)
dist_features_637 <- as.matrix(dist_features_637)
dist_features_637 <- t(dist_features_637)
identical(colnames(mic_abd_637), colnames(clinical_features_637))
identical(colnames(mic_abd_637), colnames(dist_features_637))
all_features_table_637 <- rbind(mic_abd_637, clinical_features_637, dist_features_637)

taxo_mic_637 <- read.table("../../../data/637/function/KO/637_function_taxonomy.txt", 
                           sep = "\t", header = T, row.names = 1, quote = "", comment.char = "")
taxo_mic_637 <- taxo_mic_637[, c("Taxon", "Confidence")]
taxo_new_features_637 <- read.table("taxo_new_function_feature.txt", sep = "\t", header = T, row.names = 1)
taxo_new_features_637 <- subset(taxo_new_features_637, rownames(taxo_new_features_637) %in% rownames(all_features_table_637))
taxo_all_features_637 <- rbind(taxo_mic_637, taxo_new_features_637)

all_features_table_637 <- all_features_table_637[order(rownames(all_features_table_637)), ]
taxo_all_features_637 <- taxo_all_features_637[order(rownames(taxo_all_features_637)), ]
identical(rownames(all_features_table_637), rownames(taxo_all_features_637))

meta_637 <- read.table("../../../data/637/function/637_function_metadata.tsv", sep = "\t", header = T, row.names = 1)
identical(colnames(all_features_table_637), rownames(meta_637))

all_features_table_biom_637 <- make_biom(all_features_table_637)
write_biom(all_features_table_biom_637, paste0(outdir, "/637_function_feature_table.biom"))
write.table(taxo_all_features_637, paste0(outdir, "/637_function_feature_taxon.txt"), 
            sep = "\t", quote = F, row.names = T, col.names = NA)



mic_abd_1867 <- read_hdf5_biom("../../../data/1867/function/KO/1867_function_kos_abd.biom")
mic_abd_1867 <- biom(mic_abd_1867)
mic_abd_1867 <- as.matrix(biom_data(mic_abd_1867))
col_sums <- colSums(mic_abd_1867)
mic_abd_1867 <- sweep(mic_abd_1867, 2, col_sums, FUN="/")
mic_abd_1867 <- mic_abd_1867[, order(colnames(mic_abd_1867))]
clinical_features_1867 <- read.table("../../../Results/Figure_5/1867_clinical_features_function.txt",
                                     sep = "\t", header = T, row.names = 1)
clinical_features_1867 <- as.matrix(clinical_features_1867)
clinical_features_1867 <- t(clinical_features_1867)
dist_features_1867 <- read.table("../../../Results/Figure_5/1867_all_function_dist_features.txt", 
                                 sep = "\t", header = T, row.names = 1)
dist_features_1867 <- as.matrix(dist_features_1867)
dist_features_1867 <- t(dist_features_1867)
identical(colnames(mic_abd_1867), colnames(clinical_features_1867))
identical(colnames(mic_abd_1867), colnames(dist_features_1867))
all_features_table_1867 <- rbind(mic_abd_1867, clinical_features_1867, dist_features_1867)

taxo_mic_1867 <- read.table("../../../data/1867/function/KO/1867_function_taxonomy.txt", 
                            sep = "\t", header = T, row.names = 1, quote = "", comment.char = "")
taxo_mic_1867 <- taxo_mic_1867[, c("Taxon", "Confidence")]
taxo_new_features_1867 <- read.table("taxo_new_function_feature.txt", sep = "\t", header = T, row.names = 1)
taxo_new_features_1867 <- subset(taxo_new_features_1867, rownames(taxo_new_features_1867) %in% rownames(all_features_table_1867))
taxo_all_features_1867 <- rbind(taxo_mic_1867, taxo_new_features_1867)

all_features_table_1867 <- all_features_table_1867[order(rownames(all_features_table_1867)), ]
taxo_all_features_1867 <- taxo_all_features_1867[order(rownames(taxo_all_features_1867)), ]
identical(rownames(all_features_table_1867), rownames(taxo_all_features_1867))

meta_1867 <- read.table("../../../data/1867/function/1867_function_metadata.tsv", sep = "\t", header = T, row.names = 1)
identical(colnames(all_features_table_1867), rownames(meta_1867))

all_features_table_biom_1867 <- make_biom(all_features_table_1867)
write_biom(all_features_table_biom_1867, paste0(outdir, "/1867_function_feature_table.biom"))
write.table(taxo_all_features_1867, paste0(outdir, "/1867_function_feature_taxon.txt"), 
            sep = "\t", quote = F, row.names = T, col.names = NA)

