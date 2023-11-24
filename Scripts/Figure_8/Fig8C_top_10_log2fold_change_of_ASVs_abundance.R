library(crossRanger)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(biomformat)
library(ggpubr)
library(viridis)

differencial <- function(x, y, meta, group, outdir) {
  if(!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  result_clr <- BetweenGroup.test(x, y, clr_transform=TRUE)
  write.table(result_clr, paste0(outdir, "/1867_", group, "_clr_differential_result.txt"), sep = "\t", quote = F, row.names = T, col.names = NA)
  
  sig_feature <- result_clr[which(result_clr$IfSig == "Sig"), ]
  write.table(sig_feature, paste0(outdir, "/1867_", group, "_clr_sig_feature.txt"), sep = "\t", quote = F, row.names = T, col.names = NA)
  
  t_sig_feature <- t(sig_feature)
  write.table(t_sig_feature, paste0(outdir, "/1867_", group, "_clr_sig_feature_reversed.txt"), sep = "\t", quote = F, row.names = T, col.names = NA)
}

meta <- read.table("../../data/1867/taxonomy/1867_taxonomy_metadata.tsv", 
                       sep = "\t", header = T, row.names = 1)

abd_table <- read_hdf5_biom("../../data/1867/taxonomy/silva_taxonomy_collapse/level8/feature-table.biom")
abd_table <- biom(abd_table)
abd_table <- data.frame(t(as.matrix(biom_data(abd_table))), check.names = FALSE)
abd_table <- abd_table / rowSums(abd_table)

abd_table <- abd_table[order(rownames(abd_table)), ]
meta <- meta[order(rownames(meta)), ]
identical(rownames(meta), rownames(abd_table))

features <- read.table(paste0(outpath, "/Top_10_log2fold_change_asvs_for_ECC_diagnosis.txt"),
                      sep = "\t", header = T)
features <- unique(features$feature)
taxon <- read.table("../../data/1867/taxonomy/silva_taxonomy_collapse/level8/taxonomy.tsv", 
                    sep = "\t", header = T, row.names = 1)
taxon$Taxon <- with(taxon, paste(Taxon, ASV_ID, sep = "\n"))
getLastElement <- function(s) {
  elements <- strsplit(s, ";")[[1]]
  return(tail(elements, 1))
}
taxon$Taxon <- sapply(taxon$Taxon, getLastElement)

abd_table <- abd_table[, which(colnames(abd_table) %in% features)]
colnames(abd_table) <- taxon[colnames(abd_table), "Taxon"]
abd_table$samples <- rownames(abd_table)
melt_abd_table <- melt(abd_table)
melt_abd_table$Future_Status_Tooth <- meta[melt_abd_table$samples, "Future_Status_Tooth"]
melt_abd_table$Future_Status_Tooth <- factor(melt_abd_table$Future_Status_Tooth, 
                                                          levels = c("ConfidentH", "C_H", "RelativeH", "Caries"))
melt_abd_table[which(melt_abd_table$value == 0), "value"] <- 1e-7
p <- ggplot(melt_abd_table, aes(x = Future_Status_Tooth, y = log10(value))) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color=Future_Status_Tooth), position=position_jitterdodge(jitter.width= 0.2,dodge.width = 0.8),size=1,alpha=0.4) +
  ylim(min(log10(melt_abd_table$value)), 2) + 
  scale_colour_manual(values = viridis(4)) +
  ylab("Abundance")+ xlab("The real status of a tooth")+
  geom_signif(comparisons = list(c("ConfidentH", "RelativeH"),
                                 c("ConfidentH", "Caries")),
              map_signif_level = function(p) {if(p < 0.01) {p = "**"} else if(p < 0.05) {p = "*"} else {p = "NS"}},
              test = "wilcox.test", textsize = 5, step_increase = 0.2,
              test.args = list(exact = FALSE, correct = FALSE, conf.int = TRUE, conf.level = 0.95)) +
  facet_wrap(~variable, scales="free_x", ncol = 2)+
  theme_bw()+
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        legend.position = "bottom",
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
  
ggsave(filename=paste("Fig8C_top_10_log2fold_change_of_ASVs_abundance.pdf"),plot=p, width=7, height=20)

