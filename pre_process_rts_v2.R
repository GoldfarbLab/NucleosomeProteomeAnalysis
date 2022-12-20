library(dplyr)
library(tidyverse)
library(TMTPurityCorrection)
library(MSstatsHelper)
library(MSstatsTMT)
library(plyr)
library(stringr)

### all the paths
out_dir_path<-"C:/Users/anh.h.nguyen/Documents/nucleosome/raw/rts/ap/out"
#search_path <- "C:/Users/anh.h.nguyen/Documents/nucleosome/raw/rts/ap/merged.tsv"
search_path <- "C:/Users/anh.h.nguyen/Documents/nucleosome/raw/rts/ap/MC033_HEK293T_Nucleosomes_TMTpro_Grouping2_AP_realtimesearch.tsv"
impurity_path <- "C:/Users/anh.h.nguyen/Documents/nucleosome/WB314414_2.csv"
inten_path  <- "C:/Users/anh.h.nguyen/Documents/nucleosome/raw/rts/ap/intensity.txt"
noise_path <- "C:/Users/anh.h.nguyen/Documents/nucleosome/raw/rts/ap/noise.txt"
gene_path <- "C:/Users/anh.h.nguyen/Documents/nucleosome/raw/gene_list_proteome.tsv"
#percolator_path <- "C:/Users/anh.h.nguyen/Documents/nucleosome/raw/rts/ap/percolator.target.psms.txt"
percolator_path <- "C:/Users/anh.h.nguyen/Documents/nucleosome/raw/percolator/ap/percolator.target.psms.txt"
annot_path <- file.path("C:/Users/anh.h.nguyen/Documents/nucleosome/raw/rts/ap", "annotations.csv")

### all df
impurities_df<- read.csv(impurity_path)
intensities_df <- read_tsv(inten_path) 
noise_df <- read_tsv(noise_path) 
geneName_df <- read_tsv(gene_path)
geneName_df <- geneName_df %>% dplyr::rename(ProteinName = Entry, GeneName = `Gene names`)
search_df <- read_tsv(search_path)

#names(search_df)[1] <- "Raw file"
search_df <- search_df %>% dplyr::rename(`Scan number`=`Scan Number`)
msms <- search_df %>% right_join(intensities_df, by=c("Scan number"))#, "Raw file")) #CHANGE THIS IF MERGED THEN RAW FILE SHOULD BE ON THERE

corrected_df <- correctImpurities.msms(msms, impurities_df, intensities_df, noise_df, remove.missing.rows = F)
percolator_df <- read_tsv(percolator_path)
annot_df <- read.csv(annot_path, header = TRUE)

### filter round 1: no CON, no REV, signal to noise >= 50, has bridge
filtered_df <- corrected_df %>% 
  filter(across(everything(), ~ !grepl("CON_", .))) %>% 
  filter(across(everything(), ~ !grepl("REV_", .))) %>% 
  mutate(SN_sum = rowSums(.[grep("SN \\d+", names(.))], na.rm = TRUE )) %>% 
  filter(SN_sum >= 50) #%>%
#filter(`Is missing 0` == 0) #WATCH OUT READ THE FILTERING CRITERIA!!! 


### filter in percolator search q <= 0.01
percolator_filtered <- percolator_df %>% 
  filter(`percolator q-value` <= 0.01) %>% 
  select(file_idx, scan, sequence, charge) %>%
  dplyr::rename(`Raw.file`=file_idx, `Scan.number`=scan, Peptide=sequence, Charge.State=charge) %>%
  mutate(Raw.file = as.numeric(Raw.file))
percolator_filtered$Peptide <- gsub("\\[[^)]*]", "", percolator_filtered$Peptide)


### filter round 2: tryptic sequence reformat
filtered_df_rd2 <- filtered_df %>% transform(Peptide=str_sub(Peptide,3,-3))
#runs <- unique(filtered_df_rd2[["Raw.file"]])     #run1 = 1, run2=3, run3= 0, run4=2
#idex <- as.double(c(1,3,0,2))
#filtered_df_rd2$Raw.file <- mapvalues(filtered_df_rd2$Raw.file, from=runs, to=idex)
#filtered_df_rd2 <- filtered_df_rd2 %>% mutate(Raw.file = as.numeric(Raw.file))
filtered_df_rd2$Peptide <- gsub("\\[[^)]*]", "", filtered_df_rd2$Peptide)

filtered_df_final <- filtered_df_rd2 %>% merge(percolator_filtered, by=c("Scan.number", "Peptide", "Charge.State"))#,"Raw.file")) #sometimes need to add Raw.file here
filtered_df_final$Peptide <- gsub("\\[[^)]*]", "", filtered_df_final$Peptide)

#protein group inference
# epiphany_df <- read.csv("C:/Users/anh.h.nguyen/Documents/nucleosome/raw/epiphany/out3/TOPPAS_out/003-TextExporter-out/proteingroup.txt",
#                         header=F)
# pg_df <- epiphany_df %>% select(sequence, accessions, charge,score) %>%
#   dplyr::rename(`Peptide`=sequence, ProteinName = accessions, Charge.State = charge) %>%
#   filter(score <= 0.01)
# pg_df$Peptide <- gsub("\\[[^)]*]", "", pg_df$Peptide)
# filtered_df_final3 <- filtered_df_final %>% merge(pg_df, by=c("Peptide", "Charge.State"))


#protein group inference, gene stuff starts here


get_accession <- function(string, delim = " : ", conditional_delim = "_") {
  y <- NULL
  for (i in str_split(string[[1]][1], delim)[[1]]){
    if (grepl(conditional_delim,i, fixed = T)){
      canonical <- strsplit(strsplit(i, conditional_delim, fixed = T)[[1]][1], "-",fixed = T)[[1]][1]
      y<-c(y, canonical)}
    else{
      y=c(y,i)}
  }
  return(paste(unique(y), collapse=";"))
}

get_accession2 <- function(string, delim = " : ", conditional_delim = "|") {
  y <- NULL
  for (i in str_split(string[[1]][1], delim)[[1]]){
    if (grepl(conditional_delim,i, fixed = T)){
      canonical <- strsplit(strsplit(i, conditional_delim, fixed = T)[[1]][2], "-",fixed = T)[[1]][1]
      y<-c(y, canonical)}
    else{
      y=c(y,i)}
  }
  return(paste(unique(y), collapse=";"))
}


filtered_df_final$ProteinName <- lapply(filtered_df_final[,'Protein.ID'], get_accession2)

filtered_df_final_expanded <- filtered_df_final %>% tidyr::separate_rows(ProteinName, sep=';')
filtered_df_final_expanded <- merge(filtered_df_final_expanded, geneName_df, by = "ProteinName", all.x=T)
filtered_df_final_expanded <- filtered_df_final_expanded %>% dplyr::group_by(Scan.number, Peptide, Charge.State) %>% #, Raw.file) %>% 
  dplyr::summarise(ProteinName = paste(ProteinName, collapse = ";"), GeneName = paste(GeneName, collapse = ";"))
filtered_df_final<- merge(filtered_df_final, filtered_df_final_expanded, by=c('Scan.number', 'Peptide', 'Charge.State'),all.x=T)#, 'Raw.file'), all.x=T)
filtered_df_final <- filtered_df_final %>% select(-ProteinName.y) %>% dplyr::rename(ProteinName=ProteinName.x)

#psms msstats
# annot.maxquant <- annot_df %>% 
#   filter(grepl('Bridge', BioReplicate) ) %>% 
#   mutate(BioReplicate= str_c(Mixture, '_', 'Norm' ) ) %>% 
#   full_join(annot_df) %>% 
#   filter(!grepl('Bridge', BioReplicate)) %>%
#   mutate(Condition = replace(Condition, Condition == "Bridge", "Norm"))
annot.maxquant <- annot_df
#filtered_df_final$Raw.file <- mapvalues(filtered_df_final$Raw.file, from=idex, to=runs)

msstats_tmt_df <- filtered_df_final %>%
  select( ProteinName, Peptide, Charge.State, matches("Reporter.intensity.corrected"), GeneName, Scan.number) %>% #, Raw.file) %>% #sometimes we need to add raw.file
  mutate(PSM = paste(Peptide, Charge.State, Scan.number,sep  = "_" )) %>%  #,Raw.file, sep="_")) %>% #Sometimes we need raw file here
  pivot_longer(cols=matches("Reporter.intensity.corrected"), names_to='Channel', values_to='Intensity') %>%
  mutate(Channel = str_replace(Channel,"Reporter.intensity.corrected.", "channel.")) %>%
  dplyr::rename(PeptideSequence = Peptide, Charge=Charge.State) %>% #, Run = Raw.file) %>% #SOmetimes we have raw file here
  merge(annot.maxquant, by = c('Channel')) %>%  #, "Run")) %>% #sometimes we need run here
  select(-Scan.number)




#export
msstats_tmt_df$ProteinName = as.character(msstats_tmt_df$ProteinName)
msstats_tmt_df$GeneName = as.character(msstats_tmt_df$GeneName)
msstats_tmt_df$GeneName <- gsub("NA;","",as.character(msstats_tmt_df$GeneName))
msstats_tmt_df$GeneName <- gsub(";NA","",as.character(msstats_tmt_df$GeneName))
msstats_tmt_df$GeneName <- gsub("^NA$","",as.character(msstats_tmt_df$GeneName))
msstats_tmt_df <- msstats_tmt_df %>%
  dplyr::mutate(GeneName = sapply(strsplit(GeneName, ";"), function(x) paste(unique(x), collapse = ";")))
write.csv(msstats_tmt_df, file.path(out_dir_path,"step1_msstats_psm.csv"), row.names = F, quote = F, na= "")

#export experiment design
meta <- unique( select(msstats_tmt_df, c('Condition','Mixture','Channel', 'BioReplicate')) )
meta <- meta %>%
  arrange(Mixture,factor(as.numeric((gsub("[a-zA-Z]+\\.", "", Channel)))))
meta['Class'] <- str_c(meta$BioReplicate, '_', meta$Mixture)
meta['Class2'] <- str_c(meta$Channel, '_', meta$Mixture)
write.csv(meta, file.path(out_dir_path,"experiment.csv"), row.names = F, quote = F, na= "")

#excel export
raw_noise_corrected_excel <- msstats_tmt_df %>% select("Intensity", "BioReplicate","GeneName", "ProteinName", "Mixture") %>%
  dplyr::group_by(GeneName, ProteinName, BioReplicate, Mixture) %>%
  dplyr::mutate(total = sum(Intensity)) %>%
  dplyr::select(-Intensity) %>%
  pivot_wider(names_from=c("BioReplicate", "Mixture"), values_from="total", names_sep="_", names_sort = T, values_fn = unique) %>%
  dplyr::select(c(meta$Class,"ProteinName","GeneName")) %>%
  dplyr:: select(-contains("Norm")) %>%
  dplyr::mutate(MultipleGenes = str_detect(GeneName, ";")) %>%
  arrange(MultipleGenes, GeneName) %>%
  select(-MultipleGenes)
write.csv(raw_noise_corrected_excel, file.path(out_dir_path,"xcel_step1_raw_noise_corrected.csv"), row.names = F, quote = F, na= "")



