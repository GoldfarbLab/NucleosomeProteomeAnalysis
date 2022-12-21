library(dplyr)
library(tidyverse)
library(MSstatsTMT)
library(ggfortify)
library(broom)
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(data.table)
library(factoextra)
library(cluster)


args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Specify output directory", call.=FALSE)
} 



# condition normalization -------------------------------------------------
experiment <- "grouping_all"
path<-args[1]
step1_path <- file.path(path, 'step1_msstats_psm.csv')
input <- read.csv(step1_path)

#protein inference
peptides_across_mix <- 1


fn <- input %>%
  dplyr::filter(!str_detect(GeneName, ";")) %>% #unique to gene, concern for P62805
  
  dplyr::group_by(ProteinName) %>%
  dplyr::mutate(psm_1 = uniqueN(PSM[Mixture == 1]),psm_2 = uniqueN(PSM[Mixture == 2]),
                psm_3 = uniqueN(PSM[Mixture == 3]),psm_4 = uniqueN(PSM[Mixture == 4])) %>%
  dplyr::ungroup() %>%
  dplyr::rowwise() %>%
  dplyr::mutate(psm_greater1 = sum(c_across(paste0("psm_",1:4)) > 1)) %>%
  dplyr::filter(psm_greater1 >= 2) %>% 

  dplyr::group_by(ProteinName) %>%
  dplyr::mutate(pep_1 = uniqueN(PeptideSequence[Mixture == 1]),pep_2 = uniqueN(PeptideSequence[Mixture == 2]),
                pep_3 = uniqueN(PeptideSequence[Mixture == 3]),pep_4 = uniqueN(PeptideSequence[Mixture == 4])) %>%


  # dplyr::ungroup() %>%
  # dplyr::rowwise() %>%
  # dplyr::mutate(pep_greater1 = sum(c_across(paste0("psm_",1:4)) > 1)) %>%
  # dplyr::filter(pep_greater1 ==4) %>%
  # dplyr::select(-pep_greater1) %>%
  # group_by(GeneName, Mixture) %>%
  # mutate(num_pep = n_distinct(PeptideSequence)) %>%
  # filter(num_pep >= 2) %>%
  
  dplyr::group_by(PSM) %>%   #total intensity if same psm, across all channels
  dplyr::mutate(total = sum(Intensity)) %>%
  dplyr::group_by(PeptideSequence, Charge, Run) %>%
  dplyr::filter(total == max(total)) %>%      
  dplyr::mutate(Mixture = as.factor(Mixture),
                TechRepMixture = as.factor(TechRepMixture),
                Channel = as.factor(Channel),
                BioReplicate = as.factor(BioReplicate),
                Condition = as.factor(Condition)) %>%
  dplyr::select(-total) %>%
  
  dplyr::group_by(GeneName, Mixture) %>%
  dplyr::mutate(psm =n_distinct(PSM)) %>%
  dplyr::filter(psm_greater1 >= 2 || psm > 1) %>%
  
  dplyr::group_by(Mixture, BioReplicate) %>%
  dplyr::mutate(median_run = median(Intensity)) %>%     #normalize
  dplyr::group_by(Mixture, Condition) %>%
  dplyr::mutate(median_condition = median(Intensity))%>%
  dplyr::mutate(Intensity = Intensity / median_run * median_condition)
write.csv(fn, file.path(path,"step2_normalized_msstats.csv"), row.names = F, quote = F, na= "")
#without histone
fn_without_hist <- fn %>%
  filter(!str_detect(ProteinName, ';H2A.D|;H2B.C|;H4|;xH|;H3|^H2A.D|^H2B.C|^H4|^xH|^H3'))

#excel export
exp <- read.csv(file.path(path, 'experiment.csv'))
norm_xcel <- fn %>% select("Intensity", "BioReplicate","GeneName", "ProteinName", "Mixture") %>%
  dplyr::group_by(GeneName, ProteinName, BioReplicate) %>%
  dplyr::mutate(total = sum(Intensity)) %>%
  dplyr::select(-Intensity) %>%
  pivot_wider(names_from=c("BioReplicate", "Mixture"), values_from="total", names_sep="_", 
              names_sort = T, values_fn = unique) %>%
  dplyr::select(c("ProteinName","GeneName",exp$Class)) %>%
  dplyr:: select(-contains("Norm")) %>%
  group_by(GeneName, ProteinName) %>%
  summarise_all(list(~ sum(., na.rm = T)))
temp <- fn %>% 
  ungroup() %>% 
  select(ProteinName, GeneName, matches("psm_"),matches("pep_")) %>% 
  unique()
norm_xcel <- merge(norm_xcel, temp, by = c("GeneName", "ProteinName"), all.x = T)
cols <- read.csv(file.path(path, 'xcel_step1_raw_noise_corrected.csv'))
cols <- cols %>%  select(GeneName, ProteinName)
norm_xcel2 <- merge(cols, norm_xcel, by =c("GeneName", "ProteinName"), all.y = T)
norm_xcel2 <- norm_xcel2 %>%
  dplyr::mutate(MultipleGenes = str_detect(GeneName, ";")) %>%
  arrange(MultipleGenes, GeneName) %>%
  select(-MultipleGenes)
write.csv(norm_xcel2, file.path(path,"xcel_step2_normalized_msstats.csv"), row.names = F, quote = F, na= "")

#use gene
fn_genes <- fn_without_hist %>% dplyr::rename(Protein=ProteinName)  %>% dplyr::rename(ProteinName=GeneName) 

prot.mstats_gene <- proteinSummarization(fn_genes,
                                         method="msstats",
                                         global_norm=F,
                                         reference_norm=T, remove_norm_channel=T,
                                         MBimpute = T, verbose = TRUE)

get_accession <- function(string) {
  y <- NULL
  for (i in str_split(string[[1]][1], ";")[[1]]){
    if (grepl("|",i, fixed = T)){
      y<-c(y, strsplit(i, "|", fixed = T)[[1]][2])}
    #else{
    #y=c(y,i)}
  }
  return(paste(y, collapse=";"))
}
prot.mstats$temp2$GeneName <- lapply(unlist(prot.mstats$temp2[,'Protein']), get_accession)
saveRDS(prot.mstats_gene, file = file.path(path,'step3_protein_abundance_woHist.rds'))



# profile plot ------------------------------------------------------------
prot.mstats_gene2 <- prot.mstats_gene
prot.mstats_gene2$FeatureLevelData <- prot.mstats_gene2$FeatureLevelData %>%
  mutate(Run = paste("Mixture" ,Mixture, sep = " " ))
prot.mstats_gene2$ProteinLevelData <- prot.mstats_gene2$ProteinLevelData %>%
  mutate(Run = paste("Mixture" ,Mixture, sep = " " ))

.plotProfileTMT2 <- MSstatsTMT:::.plotProfileTMT
dataProcessPlotsTMT2 <- MSstatsTMT::dataProcessPlotsTMT
prot.mstats_gene3 <- prot.mstats_gene2


prot.mstats_gene3$FeatureLevelData <- prot.mstats_gene3$FeatureLevelData %>%
  filter(ProteinName == "KDM2A")
prot.mstats_gene3$ProteinLevelData <- prot.mstats_gene3$ProteinLevelData %>%
  filter(Protein == "KDM2A")


dataProcessPlotsTMT2(prot.mstats_gene3,type= "ProfilePlot",
                     address=path, legend.size=10, width = 21, height =7)



# pca ---------------------------------------------------------------------


prot_for_pca <- prot.mstats_gene$ProteinLevelData  %>% 
  select(Abundance, Mixture, Channel, Protein) %>%
  pivot_wider(names_from=c("Channel", "Mixture"), values_from="Abundance", names_sep="_", names_sort = T)



meta <- unique( select(prot.mstats_gene$ProteinLevelData, c('Condition','Mixture','Channel')) )
meta['Class'] <- str_c(meta$Channel, '_', meta$Mixture)

prot_for_pca <- prot_for_pca %>%
  drop_na() %>%
  column_to_rownames("Protein")

pca<-prcomp(t(prot_for_pca), center = T, scale = F)

pca_fit <- augment(pca) %>%
  rename_at(vars(starts_with(".fitted")),
            list(~str_replace(.,".fitted",""))) 

pca_fit <- merge(pca_fit, meta, by.x = ".rownames", by.y = "Class")

pca.var <- facto_summarize(pca, "var", c("coord", "contrib", "cos2"), axes = c(1,2))
colnames(pca.var)[2:3] <-  c("x", "y")
pca.ind <- get_pca_ind(pca)
ind <- data.frame(pca.ind$coord[, c(1,2), drop=FALSE], stringsAsFactors = TRUE)
colnames(ind) <- c("x", "y")

r <- min(
  (max(ind[,"x"])-min(ind[,"x"])/(max(pca.var[,"x"])-min(pca.var[,"x"]))),
  (max(ind[,"y"])-min(ind[,"y"])/(max(pca.var[,"y"])-min(pca.var[,"y"])))
)
pca.var <- pca.var %>% 
  mutate(x = x*r*-0.75, y = y*r*-0.75) %>% 
  mutate(r = sqrt(x^2+y^2))
pca.var <- cbind.data.frame(pca.var, xstart = 0, ystart = 0, stringsAsFactors = TRUE)
options(ggrepel.max.overlaps = 10)
pca.var2 <- pca.var %>% 
  filter(r > 17)
pca_fit2 <- pca_fit[match(unique(pca_fit$Condition), pca_fit$Condition),]
pdf(file = file.path(path, 'PCA.pdf'))
ggplot(pca_fit, aes(x=PC1, y=PC2)) + #sometimes we need shape = Mixture
  geom_point(size=3, aes(x=PC1, y=PC2, color=Condition, shape = Mixture)) + 
  geom_vline(xintercept = 0, linetype="dotted") +
  geom_hline(yintercept = 0, linetype="dotted") +
  lemon::scale_y_symmetric(mid=0) +
  geom_label_repel(data= pca_fit2, aes(x=PC1, y=PC2, label = Condition, color = Condition),
                   min.segment.length = Inf, fill=alpha(c("white"), 0.25)) +
  #scale_shape_discrete(name = "Condition") + 
  #scale_color_discrete(guide='none') + 
  xlab(paste("PC1 (", summary(pca)$importance[2,1]*100, "% explained variance)", sep=""))+ 
  ylab(paste("PC2 (", summary(pca)$importance[2,2]*100, "% explained variance)", sep=""))+
  
  theme_clean()+
  theme(panel.grid.major.y = element_blank()) +
  scale_color_manual(breaks = c("WT","AP","Q24A","E56A","E56Q",
                                "E61A","E64A","N68A","D72A",	"N89A",
                                "D90A","E91A","E92A","E92K","Q47A",
                                "E113A","E113Q"),
                     values=c('#000000','#ff0000','#808080','#0000ff','#00ffff',
                              "#ff5dff","#ffff00",	"#b7b7ff","#ffbf00",	"#bcffbc",
                              "#00bfff",	"#ab7942",	"#00ac00",	"#00ff00","#bcbc00",
                              "#aa00ff",	"#ffb7e7","#ffffff"))+
  ggpubr::geom_exec(geom_segment, data = pca.var2,
                    x = "xstart", y = "ystart", xend = "x", yend = "y",
                    arrow = grid::arrow(length = grid::unit(0.2, 'cm')),
                    color = "steelblue", alpha = 1, size = 0.5) +
  geom_text_repel(data = pca.var2,mapping = aes( x=x,y=y, label = name), 
                  color = "steelblue")

dev.off()

# other pca ---------------------------------------------------------------------

summarization <- prot.mstats_gene$ProteinLevelData %>% 
  select("Abundance", "Channel","Protein", "Mixture") %>%
  dplyr::rename(ProteinName = Protein) %>%
  pivot_wider(names_from=c("Channel", "Mixture"), values_from="Abundance", names_sep="_", names_sort = T, values_fn = unique)
names(summarization)[1] <- "GeneName"
summarization <- merge(cols, summarization, by =c("GeneName"), all.y = T)
summarization <- summarization %>%
  dplyr::mutate(MultipleGenes = str_detect(GeneName, ";")) %>%
  arrange(MultipleGenes, GeneName) %>%
  select(-MultipleGenes)
write.csv(summarization, file.path(path,"xcel_step3_abundance_woHist.csv"), row.names = F, quote = F, na="")



write.csv(pca$x, file.path(path,"pca_coord.csv"), row.names = T, quote = F, na="")
write.csv(pca$rotation, file.path(path,"eigen.csv"), row.names = T, quote = F, na="")
write.csv(pca.var, file.path(path,"eigen_scaled.csv"), row.names = T, quote = F, na="")

