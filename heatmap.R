library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(matrixStats)
library(tidyverse)
library(circlize)
library(cluster)
library(RColorBrewer)
library(ggpubr)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Specify which experiment to plot and output directory", call.=FALSE)
} else if (length(args)==1) {
  stop("2 arguments required. Specify which experiment to plot and output directory", call.=FALSE)
}

experiment_types <- c(
  "grouping_all" = file.path("data", "grouping1_firstAttempt"),
  "v1" = file.path("data", "v1"),
  "walk" = file.path("data", "walk"),
  "ap" = file.path("data", "ap")
)
experiment <- args[1]
path<-args[2]
mstats_quant <- readRDS(file.path(path, 'step3_protein_abundance_woHist.rds'))
geneOfInterest <- unlist(lapply(c("ACTL6A", "ACTR6", "ANAPC1", "ANAPC2", "anapc4", 
                                  "anapc5", "bard1", "baz1b", "cdc16","cdc23", "cdc27", "dek", 
                                  "dido1", "dmap1", "dnttip1", "hmgn2", "kdm2a", "kdm2b", 
                                  "kif20b","mbip", "mdc1", "parp1", "rad50", "ran", "rcc1", "rnf2", 
                                  "ruvbl1", "ruvbl2", "smarca5", "tada2a", 
                                  "trip12", "vps72", "vrk1", "wdr76", "yeats2", "yeats4", "zzz3"), toupper))

prot.mstats <- mstats_quant$ProteinLevelData %>% 
  dplyr::select(Abundance, Mixture, Channel, Protein) %>%
  dplyr::group_by(Protein) %>%
  dplyr::mutate(present_mixture = length(unique(Mixture))) %>%
  pivot_wider(names_from=c("Channel", "Mixture"), values_from="Abundance", names_sep="_", names_sort = T)
design <- read.csv(file.path(path, 'experiment.csv'))
design <- design %>% dplyr::filter(!str_detect(BioReplicate, "Norm"))


# log2 intensities -----------------------------------------------------------
if (experiment == "grouping_all"){
  prot.mstats <- prot.mstats %>%
    filter(present_mixture >= 2)
  prot.mstats.fc <- matrix(,nrow=dim(prot.mstats)[1],ncol=0)
  wt_list <- list()
  inten.order <- c()
  for (mix in 1:4){
    wt <- paste('channel.', 14:(14+2), '_', mix, sep='')
    wt_med <- rowMedians(base::as.matrix(prot.mstats[, wt]), na.rm = T)
    wt_list[[mix]] <- wt_med
    #for (sample in seq(2,13,3)){
    for (sample in sequence(4,2,3)){
      inten.order <- append(inten.order, c(paste('channel.', sample:(sample+2), '_', mix, sep='')))
      grouping <- c(paste('channel.', sample:(sample+2), '_', mix, sep=''))
      for (i in 1:3){
        prot.mstats.fc<-cbind(prot.mstats.fc, wt_med)
      }
    }
  }
} else{
  prot.mstats.fc <- matrix(,nrow=dim(prot.mstats)[1],ncol=0)
  wt_list <- list()
  inten.order <- c()
  wt <- paste('channel.', 14:(14+2), '_', 1, sep='')
  wt_med <- rowMedians(base::as.matrix(prot.mstats[, wt]), na.rm = T)
  wt_list[[1]] <- wt_med
  #for (sample in seq(2,13,3)){
  mix <- 1
  for (sample in sequence(4,2,3)){
    inten.order <- append(inten.order, c(paste('channel.', sample:(sample+2), '_', mix, sep='')))
    grouping <- c(paste('channel.', sample:(sample+2), '_', mix, sep=''))
    for (i in 1:3){
      prot.mstats.fc<-cbind(prot.mstats.fc, wt_med)
    }
  }
  
}
#prot.mstats.fc <- matrix(,nrow=dim(prot.mstats)[1],ncol=0)
# 61A 2_2, 90A 5_3, 92A 2_4, 64A 5_2, 56A 8_1, 113A 8_4	
# inten.order <- c("channel.2_2","channel.3_2","channel.4_2",
#                  "channel.5_3","channel.6_3","channel.7_3",
#                  "channel.2_4","channel.3_4","channel.4_4",
#                  "channel.5_2","channel.6_2","channel.7_2",
#                  "channel.8_1","channel.9_1","channel.10_1",
#                  "channel.8_4","channel.9_4","channel.10_4")
# prot.mstats.fc <- cbind(prot.mstats.fc, wt_list[[2]],wt_list[[2]],wt_list[[2]],
#                         wt_list[[3]],wt_list[[3]],wt_list[[3]],
#                         wt_list[[4]],wt_list[[4]],wt_list[[4]],
#                         wt_list[[2]],wt_list[[2]],wt_list[[2]],
#                         wt_list[[1]],wt_list[[1]],wt_list[[1]],
#                         wt_list[[4]],wt_list[[4]],wt_list[[4]])
# 92K 5_4, 56Q 11_1, 113Q		11_4
# 92A 2_4, 56A 8_1,  113A   8_4
# inten.order <- c("channel.5_4","channel.6_4","channel.7_4",
#                  "channel.2_4","channel.3_4","channel.4_4",
#                 "channel.11_1","channel.12_1","channel.13_1",
#                 "channel.8_1","channel.9_1","channel.10_1",
#                   "channel.11_4","channel.12_4","channel.13_4",
#                 "channel.8_4","channel.9_4","channel.10_4")
# prot.mstats.fc <- cbind(prot.mstats.fc, wt_list[[4]],wt_list[[4]],wt_list[[4]],
#                         wt_list[[4]],wt_list[[4]],wt_list[[4]],
#                         wt_list[[1]],wt_list[[1]],wt_list[[1]],
#                         wt_list[[1]],wt_list[[1]],wt_list[[1]],
#                         wt_list[[4]],wt_list[[4]],wt_list[[4]],
#                         wt_list[[4]],wt_list[[4]],wt_list[[4]])
# inten.order <- c("channel.5_4","channel.6_4","channel.7_4",
#                  "channel.11_1","channel.12_1","channel.13_1",
#                  "channel.11_4","channel.12_4","channel.13_4")
# prot.mstats.fc <- cbind(prot.mstats.fc, wt_list[[4]],wt_list[[4]],wt_list[[4]],
#                         wt_list[[1]],wt_list[[1]],wt_list[[1]],
#                         wt_list[[4]],wt_list[[4]],wt_list[[4]])
  # 61A 2_2, 90A 5_3, 92A 2_4, 92K		5_4
  # 56A 8_1, 56Q 11_1, 113A 8_4, 113Q		11_4
  # inten.order <- c("channel.2_2","channel.3_2","channel.4_2",
  #                  "channel.5_3","channel.6_3","channel.7_3",
  #                  "channel.2_4","channel.3_4","channel.4_4",
  #                  "channel.5_4","channel.6_4","channel.7_4")
  # prot.mstats.fc <- cbind(prot.mstats.fc, wt_list[[2]],wt_list[[2]],wt_list[[2]],
  #                         wt_list[[3]],wt_list[[3]],wt_list[[3]],
  #                         wt_list[[4]],wt_list[[4]],wt_list[[4]],
  #                         wt_list[[4]],wt_list[[4]],wt_list[[4]])
# inten.order <- c("channel.8_1","channel.9_1","channel.10_1",
#                  "channel.11_1","channel.12_1","channel.13_1",
#                  "channel.8_4","channel.9_4","channel.10_4",
#                  "channel.11_4","channel.12_4","channel.13_4")
# prot.mstats.fc <- cbind(prot.mstats.fc, wt_list[[1]],wt_list[[1]],wt_list[[1]],
#                         wt_list[[1]],wt_list[[1]],wt_list[[1]],
#                         wt_list[[4]],wt_list[[4]],wt_list[[4]],
#                         wt_list[[4]],wt_list[[4]],wt_list[[4]])



# z score stuff -----------------------------------------------------------


# prot.mstats.avg <- matrix(,nrow=295,ncol=0)
# prot.mstats.sds <- matrix(,nrow=295,ncol=0)
# inten.order <- c()
# for (mix in 1:4){
#   wt <- paste('channel', 14:(14+2), '_', mix, sep='')
#   for (sample in seq(2,13,3)){
#     inten.order <- append(inten.order, c(paste('channel', sample:(sample+2), '_', mix, sep='')))
#     grouping <- append(c(paste('channel', sample:(sample+2), '_', mix, sep='')), wt )
#     for (i in 1:3){
#       prot.mstats.avg<-cbind(prot.mstats.avg, rowMeans(prot.mstats[,grouping]))
#       prot.mstats.sds<-cbind(prot.mstats.sds, rowSds(as.matrix(prot.mstats[,grouping])))
#   }
#     }
# }
# 
# 

# heat map stuff -----------------------------------------------------------
data_matrx <- data.matrix(prot.mstats)
inten <- data_matrx[,inten.order]

#prot.z.scores <- (inten-prot.mstats.avg)/prot.mstats.sds
#prot.z.scores[is.na(prot.z.scores)] <- 0
prot.log2.fc <- inten-prot.mstats.fc
rownames(prot.log2.fc) <- prot.mstats$Protein
# prot.log2.fc <- tibble::rownames_to_column(as.data.frame(prot.log2.fc), "ProteinName") 
# prot.log2.fc <- arrange(prot.log2.fc, ProteinName)
# write.csv(prot.log2.fc, file.path(path,"xcel_step4_volcano.csv"), row.names = F, quote = F, na= "")


#clustering stuff
gower_cluster <- daisy(prot.log2.fc, metric = "gower", warnType = T)
gower_cluster[is.na(gower_cluster)] <- 0

#prot.log2.fc.t <- t(prot.log2.fc)

#plot(hclust(daisy(matrix(as.vector(prot.log2.fc.t), nrow = 16), metric = "gower", warnType = T)))
# plot(hclust(daisy(prot.log2.fc.t, metric = "gower", warnType = T)))
# 
# 
# temp <- ckmeans(prot.log2.fc.t, 2, mustLink = matrix(c(43,44),nrow = 1), cantLink = matrix(c(45,47),nrow = 1))
# plot(hclust(temp))

#vertical clust
# group.medians <- matrix(,nrow=dim(prot.mstats)[1],ncol=0)
# for (mix in sequence(6,1,3)){ #full: sequence(16,1,3)
#   group.medians <- cbind(group.medians, rowMeans(prot.log2.fc[,c(mix,mix+1,mix+2)]) )
# }
# v.clust <- hclust(daisy(t(group.medians), metric = "gower", warnType = T))
# plot(v.clust)
# 
# ord <- v.clust$order
# ord.func <- function(x) c((x-1)*3+1, (x-1)*3+2, (x-1)*3+3)
# 
# prot.log2.fc <- prot.log2.fc[, unlist(lapply(ord, ord.func))]


gene.names<- prot.mstats$Protein
prot.names3 <- c()
prot.names4 <- c()
for (p in gene.names){
  if(p%in% geneOfInterest){
    prot.names3 <- c(prot.names3,p)
    prot.names4 <- c(prot.names4,p)
  }
  else{
    prot.names3<- c(prot.names3,"")
  }
}
prot.names3.pos <- which(is.na(str_detect("",prot.names3))==FALSE)

label <- design[design$Class2 %in% colnames(prot.log2.fc),] %>% 
  arrange(factor(Class2, levels= colnames(prot.log2.fc))) %>%
  select(Condition)

label


label <- as.list(label)
col_fun <- colorRamp2(seq(-5, 5, length=256), rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)))
if (experiment == "grouping_all"){
ha <- HeatmapAnnotation(mutants= c(label[[1]]), show_annotation_name = F,
                        col=list(mutants = c("AP"="#ff0000","Q24A"='#808080',
                                             "E56A"='#0000ff',"E56Q"='#00ffff',
                                             "E61A"='#ff5dff',	"E64A"="#ffff00",	"N68A"="#b7b7ff",	"D72A"="#ffbf00",
                                             "N89A"="#bcffbc",	"D90A"="#00bfff",	"E91A"="#ab7942",
                                             "E92A"="#00ac00",	"E92K"="#00ff00",
                                             "Q47A"="#bcbc00",	"E113A"="#aa00ff",	"E113Q"="#ffb7e7") ))
} else{
  ha <- HeatmapAnnotation(mutants= c(label[[1]]), show_annotation_name = F)
}
har <- rowAnnotation(gene.name = anno_mark(at = prot.names3.pos, labels = prot.names4, labels_gp = gpar(fontsize = 5)))
hm <- Heatmap(prot.log2.fc, height = unit(30, "cm"),
              row_names_gp = gpar(fontsize = 3),
              #right_annotation = har,
              name='log2(mutant/wt)', cluster_columns = F, cluster_rows = hclust(gower_cluster),
              col=col_fun,
              top_annotation, show_column_names = F, na_col = "darkgrey",
              show_row_names = T,
              heatmap_legend_param = list(direction = "horizontal"),
              row_dend_reorder = T, row_split= 4, column_title = sprintf("%s proteins", nrow(prot.log2.fc)))
pdf(file = file.path(path, 'full_heatmap.pdf'), height = 15)
draw(hm,  heatmap_legend_side = "bottom")
dev.off()
#write.csv(prot.mstats[,c('Protein',inten.order)], 'heatmap.csv')
#"WT"	"AP"	"Q24A"	"E56A"	"E56Q"	"E61A"	"E64A"	"N68A"	"D72A"	
#"N89A"	"D90A"	"E91A"	"E92A"	"E92K"	"Q47A"	"E113A"	"E113Q"
#000000	#ff0000	#808080	#0000ff	#00ffff	#ff5dff	#ffff00	#b7b7ff	
#ffbf00	#bcffbc	#00bfff	#ab7942	#00ac00	#00ff00	#bcbc00	#aa00ff	#ffb7e7


compare_excel_final <- readRDS(file.path(path,"step4_compare_woHist.rds"))$filtered
compare_excel_final <- compare_excel_final %>% select("Protein", "Label", "pvalue", "log2FC","adj.pvalue", "log2FC") %>% 
  pivot_wider(names_from=c("Label"), values_from=c("pvalue", "log2FC", "adj.pvalue"), names_sep="_", names_sort = T, values_fn = unique)

#heatmap
p_val_anova <- compare_excel_final[grep("^pvalue_|Protein", names(compare_excel_final))]
AnovaPVal <- manova(cbind(`pvalue_AP-WT`,`pvalue_Q24A-WT`,`pvalue_E56A-WT`,`pvalue_E56Q-WT`,
                          `pvalue_E61A-WT`,`pvalue_E64A-WT`,`pvalue_N68A-WT`,`pvalue_D72A-WT`,
                          `pvalue_N89A-WT`,`pvalue_D90A-WT`,`pvalue_E91A-WT`,`pvalue_Q47A-WT`,
                          `pvalue_E92A-WT`, `pvalue_E92K-WT`, `pvalue_E113A-WT`,`pvalue_E113Q-WT`) ~ Protein, 
                    data = p_val_anova)
cols <- data.frame(GeneName = rownames(prot.log2.fc))


anova_p_genes <- AnovaPVal$model %>% dplyr::rename(GeneName = Protein)
stats_pval <- merge(cols, anova_p_genes, by =c("GeneName"), all = T)
stats_pval <- stats_pval %>%
  dplyr::mutate(MultipleGenes = str_detect(GeneName, ";")) %>%
  arrange(MultipleGenes, GeneName) %>%
  select(-MultipleGenes)
write.csv(stats_pval, file.path(path,"xcel_step5_pval.csv"), row.names = F, quote = F, na= "")

#write.csv(pval_excel_final, file.path(path,"step5_pval.csv"), row.names = F, quote = F, na="")

q_val_anova <- compare_excel_final[grep("^adj.pvalue_|Protein", names(compare_excel_final))]
AnovaPVal2 <- manova(cbind(`adj.pvalue_AP-WT`,`adj.pvalue_Q24A-WT`,`adj.pvalue_E56A-WT`,`adj.pvalue_E56Q-WT`,
                           `adj.pvalue_E61A-WT`,`adj.pvalue_E64A-WT`,`adj.pvalue_N68A-WT`,`adj.pvalue_D72A-WT`,
                           `adj.pvalue_N89A-WT`,`adj.pvalue_D90A-WT`,`adj.pvalue_E91A-WT`,`adj.pvalue_Q47A-WT`,
                           `adj.pvalue_E92A-WT`, `adj.pvalue_E92K-WT`, `adj.pvalue_E113A-WT`,`adj.pvalue_E113Q-WT`) ~ Protein, 
                     data = q_val_anova)
anova_p_genes <- AnovaPVal2$model %>% dplyr::rename(GeneName = Protein)
stats_qval <- merge(cols, anova_p_genes, by =c("GeneName"), all = T)
stats_qval <- stats_qval %>%
  dplyr::mutate(MultipleGenes = str_detect(GeneName, ";")) %>%
  arrange(MultipleGenes, GeneName) %>%
  select(-MultipleGenes)
write.csv(stats_qval, file.path(path,"xcel_step5_qval.csv"), row.names = F, quote = F, na= "")


fc_anova <- compare_excel_final[grep("^log2FC_|Protein", names(compare_excel_final))]
AnovaFC2 <- manova(cbind(`log2FC_AP-WT`,`log2FC_Q24A-WT`,`log2FC_E56A-WT`,`log2FC_E56Q-WT`,
                           `log2FC_E61A-WT`,`log2FC_E64A-WT`,`log2FC_N68A-WT`,`log2FC_D72A-WT`,
                           `log2FC_N89A-WT`,`log2FC_D90A-WT`,`log2FC_E91A-WT`,`log2FC_Q47A-WT`,
                           `log2FC_E92A-WT`, `log2FC_E92K-WT`, `log2FC_E113A-WT`,`log2FC_E113Q-WT`) ~ Protein, 
                     data = fc_anova)
anova_p_genes <- AnovaFC2$model %>% dplyr::rename(GeneName = Protein)
stats_fc <- merge(cols, anova_p_genes, by =c("GeneName"), all = T)
stats_fc <- stats_fc %>%
  dplyr::mutate(MultipleGenes = str_detect(GeneName, ";")) %>%
  arrange(MultipleGenes, GeneName) %>%
  select(-MultipleGenes)
write.csv(stats_fc, file.path(path,"xcel_step5_fc.csv"), row.names = F, quote = F, na= "")


#anova for q val


anova2 <- rownames_to_column( data.frame(t(prot.log2.fc)), "condition")
design2 <- design %>% 
  select(Class2, Condition) %>% 
  dplyr::rename(condition = Class2) %>% 
  unique()
anova2 <- merge(anova2, design2, by = "condition", all = T)
anova2 <- anova2 %>% 
  select( -condition) %>% 
  dplyr::rename(condition = Condition)
anova2$condition <- factor(anova2$condition)
formulae <- lapply(colnames(anova2)[1:ncol(anova2)-1], function(x) as.formula(paste0(x, " ~ condition")))
res <- lapply(formulae, function(x) summary(aov(x, data = anova2)))
names(res) <- format(formulae)
res
p <- unlist(lapply(res, function(x) x[[1]]$"Pr(>F)"[1]))
p
pval <- data.frame(
  Gene = sub(' ~ condition', '', names(p)),
  pvalue = p)
pval <-arrange(pval, Gene)
write.csv(pval, file.path(path,"xcel_step6_hmpval.csv"), row.names = F, quote = F, na= "")

pvaladj <- data.frame( 
  Gene = pval$Gene,
  padj = p.adjust(pval$pvalue))
write.csv(pvaladj, file.path(path,"xcel_step6_hmpvaladj.csv"), row.names = F, quote = F, na= "")
