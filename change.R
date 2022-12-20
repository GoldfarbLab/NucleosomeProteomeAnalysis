library(MSstatsTMT)
library(MSstats)
library(EnhancedVolcano)
library(dplyr)
library(tidyverse)
library(UpSetR)
library(ComplexHeatmap)
library(conclust)
path<-"C:/Users/anh.h.nguyen/Documents/nucleosome/raw/rts/grouping1_firstAttempt/out/"
prot.mstats <- readRDS(file.path(path, 'step3_protein_abundance_woHist.rds'))
#prot.mstats$ProteinLevelData <- prot.mstats$temp2 #no histone



# heat map---------------------------------------------------------------------
geneOfInterest <- unlist(lapply(c("ACTL6A", "ACTR6", "ANAPC1", "ANAPC2", "anapc4", 
                                  "anapc5", "bard1", "baz1b", "cdc16","cdc23", "cdc27", "dek", 
                                  "dido1", "dmap1", "dnttip1", "hmgn2", "kdm2a", "kdm2b", 
                                  "kif20b","mbip", "mdc1", "parp1", "rad50", "ran", "rcc1", "rnf2", 
                                  "ruvbl1", "ruvbl2", "smarca5", "tada2a", 
                                  "trip12", "vps72", "vrk1", "wdr76", "yeats2", "yeats4", "zzz3"), toupper))
no_norms = levels(prot.mstats$ProteinLevelData$Condition)[levels(prot.mstats$ProteinLevelData$Condition) != "Norm"]
compare_matrix <- cbind(1 * diag(length(no_norms)-1), rep(-1, length(no_norms)-1))
rownames(compare_matrix) <- str_c(no_norms[1:16], "-WT")
colnames(compare_matrix) <- no_norms[no_norms != "Norm"]

test.maxquant.pairwise <- groupComparisonTMT(
  data = prot.mstats,
  contrast.matrix = compare_matrix,
  moderated = T,
  # do moderated t test
  adj.method = "BH",
  remove_norm_channel = T,
  remove_empty_channel = T
)
test.maxquant.pairwise$filtered <- test.maxquant.pairwise$ComparisonResult %>% 
  filter(is.na(issue))  #FIXME: Only include the one that is full plate

saveRDS(test.maxquant.pairwise, file = file.path(path,'step4_compare_woHist.rds'))

#test.maxquant.pairwise$ComparisonResult <- test.maxquant.pairwise$ComparisonResult[complete.cases(test.maxquant.pairwise$ComparisonResult[ , 3:7]),]
#Volcano Plot
test.maxquant.pairwise <- readRDS(file.path(path,'step4_compare_woHist.rds'))
condition <- unique(test.maxquant.pairwise$filtered[,'Label'][[1]])
for (c in condition ){
  subDf <- test.maxquant.pairwise$filtered %>% filter(str_detect(c,Label))
  pdf(file = file.path(path, paste(c,'.pdf',sep='')))
  plot(EnhancedVolcano(subDf,
                  lab = levels(subDf$Protein)[subDf$Protein],
                  x = 'log2FC',
                  y = 'adj.pvalue',
                  title = c,
                  subtitle ='',
                  selectLab = geneOfInterest,
                  xlab = bquote(~Log[2]~ 'fold change'),
                  ylab = bquote(~-Log[10]~ 'P(adj.pvalue)'),
                  pCutoff =  0.05,
                  FCcutoff = 1.0,
                  pointSize = 4.0,
                  labSize = 3.0,
                  labCol = 'black',
                  labFace = 'bold',
                  boxedLabels = TRUE,
                  colAlpha = 4/5,
                  legendPosition = 'bottom',
                  legendLabSize = 14,
                  legendIconSize = 4.0,
                  gridlines.major = FALSE,
                  gridlines.minor = FALSE,
                  drawConnectors = TRUE,
                  widthConnectors = 1.0,
                  colConnectors = 'black'))
  dev.off()
}


change <- test.maxquant.pairwise$filtered %>% select(Protein, Label, log2FC) %>%
  pivot_wider(names_from="Label", values_from=c("log2FC"))
change$Protein <- as.character(change$Protein)
write.csv(change, file = file.path(path,'step4_fc_woHist.csv'), quote = F, row.names = F, na = "")

#export xcel
condition <- unique(read.csv(file.path(path, 'experiment.csv'))$Condition)
condition <- rep(append(condition[2:5], condition[7:18]), each = 3)
tag <- rep(c("pvalue_", "adj.pvalue_", "log2FC_"), times = 16)
condition <- paste(tag, condition,sep="")

stats_vc <- test.maxquant.pairwise$filtered %>% 
  select(Protein, pvalue,`adj.pvalue`, log2FC, Label) %>% 
  mutate(Label = str_remove(Label, "-WT")) %>%
  pivot_wider(names_from = Label, values_from = c('pvalue','adj.pvalue', "log2FC")) %>%
  dplyr::select(c("Protein", condition)) %>%
  dplyr::rename(GeneName = Protein)
cols <- read.csv(file.path(path, 'xcel_step1_raw_noise_corrected.csv'))
cols <- cols %>%  select(GeneName, ProteinName)
stats_vc2 <- merge(cols, stats_vc, by =c("GeneName"), all = T)
stats_vc2 <- stats_vc2 %>%
  dplyr::mutate(MultipleGenes = str_detect(GeneName, ";")) %>%
  arrange(MultipleGenes, GeneName) %>%
  select(-MultipleGenes)
write.csv(stats_vc2, file.path(path,"xcel_step4_volcano.csv"), row.names = F, quote = F, na= "")


mix_set <- test.maxquant.pairwise$filtered %>% 
  select(Label, log2FC, pvalue, Protein, adj.pvalue ) %>%
  mutate(Label = str_remove(Label, "-WT"))
mix_set$binding <- ifelse(mix_set$log2FC > 1 & mix_set$adj.pvalue < 0.05,"increased", 
                          ifelse(mix_set$log2FC < -1 & mix_set$adj.pvalue < 0.05,"decreased",
                                 "unchanged")) #FIXME: BUG HERE
pdf(file = file.path(path, 'barplot.pdf'))
ggplot(mix_set, aes(fill=binding, x=Label)) + 
  geom_bar(position='stack', stat='count',colour="black")+
  geom_text(stat = 'count',aes(label=..count..))+
  scale_fill_manual(values = c(unchanged="white", 
                               increased="darkgrey",
                      decreased="black") )+
  theme_clean() +
  labs(y = "number of proteins", x = "mutant")
dev.off()

mix_set_2 <- mix_set %>% 
  group_by(binding, Label) %>% 
  summarize(count = n())

write.csv(mix_set_2, file.path(path,"stats_barplot.csv"), row.names = F, quote = F, na= "")


full_plate <- mix_set %>% 
  select(Protein, Label) %>%
  group_by(Protein) %>%
  dplyr::summarise(count = n_distinct(Label)) %>%
  filter(count == 16)
full_set <- mix_set %>%
  filter(Protein %in% full_plate$Protein)


upset_set <- full_set %>%
  filter(full_set$log2FC < -1 & full_set$adj.pvalue < 0.05)
upset_set <- split(upset_set,upset_set$Label,
                 drop = T)
upset_set <- map(upset_set, ~ (.x %>% select(Protein)))
upset_set <- map(upset_set, unlist)
upset_set <- map(upset_set, unique)
upset_matrix <- make_comb_mat(upset_set)
upset_matrix <-upset_matrix[comb_size(upset_matrix) >1 ]
pdf(file = file.path(path, 'upset.pdf'))
UpSet(upset_matrix, comb_order = order(-comb_size(upset_matrix)))#, order.by = "freq", set_size.show=T)
dev.off()
