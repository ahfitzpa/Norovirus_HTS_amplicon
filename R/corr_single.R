library(dplyr)
library(tidyr)
library(stringr)

single_spike_samples <- c("AG1.LP2", "AG2.LP2", "AG3.LP2", "AG4.LP2","AG5.LP2", "AG12.LP2",
                          "KH1.LP2", "KH2.LP2", "KH3.LP2","KH4.LP2","KH5.LP2","KH12.LP2",
                          "KR12.LP2","KR1.LP2","KR2.LP2", "KR2.LP2","KR3.LP2","KR4.LP2", "KR5.LP2")

################################################################################################################
# read in biom file 
format_biom <- function(feature_table,pipeline){
  biomformat::read_biom(feature_table) %>%
    biomformat::biom_data(.)%>%
    as.matrix(.)%>%
    as_tibble(.,rownames = NA)%>%
    tibble::rownames_to_column("OTU") %>%
    mutate(OTU= paste0(OTU, ".", pipeline))%>%
    rename_at(vars(-OTU),function(x) paste0(x,".", pipeline)) %>% 
    mutate_if(is.numeric , replace_na, replace = 0)
}

files <- c("raw_data/LP1/all.otutab.biom","raw_data/LP2/all.otutab.biom", full.names = TRUE)
biom_list <- list()  
name <- c()
output <- for (i in seq_along(files)) {
  name[i] <- str_replace(files[[i]], "raw_data/", "") %>%
    str_replace(.,"/all.otutab", "") %>%
    str_replace(., ".biom", "")
  biom_list[[i]] <-  format_biom(files[[i]], name[i]) 
}
biom <- data.table::rbindlist(biom_list,fill=T) %>%
  mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%
  filter(purrr::reduce(across(where(is.numeric), ~.x!=0), `|`)) 


#############################################################################################################################################
format_metadata <- function(table, experiment){
  read.table(table,sep ="\t", header=T)%>%
    mutate(SampleID= paste0(SampleID, ".", experiment))
}

files <- c("raw_data/LP1/metadata.tsv","raw_data/LP2/metadata.tsv")
metadata_list <- list()  
name <- c()
output <- for (i in seq_along(files)) {
  name[i] <- str_replace(files[[i]], "raw_data/", "") %>%
    str_replace(.,"/metadata", "") %>%
    str_replace(., ".tsv", "")
  metadata_list[[i]] <-  format_metadata(files[[i]], name[i]) 
}

metadata <- data.table::rbindlist(metadata_list,fill=T,idcol=T)%>%
  mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%
  #filter(purrr::reduce(across(where(is.numeric), ~.x!=0), `|`)) %>%
  select(.id,SampleID, DNA_polymerase, RTase, EXP_ID, quant_reading, qPCR)

#############################################################################################################
#############################################################################################################################################
## CLASSIFICATION DATA ##
format_classified <- function(table,pipeline){
  read.table(table, sep ="\t", header=T) %>%
    select(OTU.ID,blast, capsid_type, genotype_capsid)%>%
    mutate(blast=str_remove(blast, "Caliciviridae Norovirus ")) %>%
    dplyr::rename(OTU=OTU.ID,genogroup=blast, genotype=capsid_type)%>%
    mutate(source =pipeline)
}

files <- c("raw_data/LP1/noronet.tsv","raw_data/LP2/noronet.tsv")
classified_list <- list()  
name <- c()
output <- for (i in seq_along(files)) {
  name[i] <- str_replace(files[[i]], "raw_data/", "") %>%
    str_replace(.,"/noronet", "") %>%
    str_replace(., ".tsv", "")%>%
    gsub(".*/","", .)
  classified_list[[i]] <-  format_classified(files[[i]], name[i]) 
}

classified <- data.table::rbindlist(classified_list,fill=T,idcol=T) %>%
  mutate(genotype_capsid =str_replace(genotype_capsid, " ", "_")) %>%
  unite(taxonomy,c("genogroup", "genotype", "genotype_capsid"), sep=";", remove=T) %>%
  select(.id, OTU, taxonomy) %>%
  mutate(source="observed") %>%
  select(.id, OTU,taxonomy, everything()) %>%
  mutate(taxonomy=ifelse(is.na(taxonomy)==T, OTU, taxonomy)) 
#############################################################################################################################################
format_expected <- function(table, experiment){
  read.table(table,sep ="\t", header=T)%>%
    rename_at(vars(-taxonomy),function(x) paste0(x,".", experiment)) %>%
    rename_at(vars(-taxonomy),function(x) str_remove(x, "^X"))
}

files <- c("raw_data/LP1/expected_composition.tsv","raw_data/LP2/expected_composition.tsv")
expected_list <- list()  
name <- c()
output <- for (i in seq_along(files)) {
  name[i] <- str_replace(files[[i]], "raw_data/", "") %>%
    str_replace(.,"/expected_composition", "") %>%
    str_replace(., ".tsv", "")
  expected_list[[i]] <-  format_expected(files[[i]], name[i]) 
}
expected <- data.table::rbindlist(expected_list,fill=T, idcol=T)%>%
  mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%
  mutate(across(where(is.numeric), ~ . *100)) %>%
  filter(purrr::reduce(across(where(is.numeric), ~.x!=0), `|`)) %>%
  mutate(.id=.id/100)%>%
  mutate(source = "expected")%>%
  dplyr::mutate(OTU=taxonomy) %>%
  mutate(taxonomy=strex::str_after_nth(taxonomy, ";", 4)) %>%
  mutate(expected_taxonomy=taxonomy)%>%
  select(!contains(single_spike_samples)) %>%
  select(-source, -OTU) %>%
  pivot_longer(!c(".id", "taxonomy", "expected_taxonomy"),names_to="SampleID", values_to="RA") %>%
  filter(RA>0) %>%
  select(-RA)

###############################################################################################################################################
biom_observed_expected <- biom %>%
  separate(OTU, c("OTU", "library"), "\\.") %>% select(-library) %>%
  pivot_longer(!OTU,names_to="SampleID", values_to="abundance") %>%
  filter(!SampleID %in% single_spike_samples) %>%
  filter(abundance >1) %>%
  full_join(., classified, by=c("OTU")) %>%
  right_join(., expected) %>%
  drop_na() %>%
  select(.id, OTU,taxonomy,source,everything()) %>%
  mutate(taxonomy=ifelse(is.na(taxonomy)==T, OTU, taxonomy)) %>%
  left_join(., metadata) %>%
  group_by(SampleID)%>%
  filter(taxonomy %in% expected_taxonomy) %>%
  ungroup() %>% distinct() %>%
  select(.id, SampleID, taxonomy, expected_taxonomy, abundance, qPCR, quant_reading, DNA_polymerase, RTase) %>%
  unite(source, c("RTase", "DNA_polymerase"), sep="_", remove=F) %>%
  group_by(SampleID) %>%
  mutate(reads=sum(abundance)) %>%
  select(-abundance, -taxonomy, -DNA_polymerase, -RTase) %>% distinct() %>%
  mutate(qPCR=log10(qPCR))%>%
  mutate(reads=log10(reads)) %>%
  ungroup()

############################################################################################################################
## CORRELATION ANALYSIS
cor_qPCR <- biom_observed_expected %>% select(.id,SampleID,qPCR, reads, source)%>%
  group_by(.id,source)%>%
  summarize(cor=stats::cor.test(qPCR, reads, method="kendall")$estimate, 
            p.value_qPCR=cor.test(qPCR, reads, method="kendall")$p.value,
            z_qqPCR=cor.test(qPCR, reads, method="kendall")$statistic) %>%
  dplyr::rename("Kendall correlation between gc/g of norovirus detected by RT-qPCR and reads after QC"=cor)

cor_qubit <- biom_observed_expected %>% select(.id,SampleID, quant_reading, reads, source)%>%
  group_by(.id,source)  %>%
  summarize(cor=stats::cor.test(quant_reading, reads, method="kendall")$estimate, 
            p.value_qubit=cor.test(quant_reading, reads, method="kendall")$p.value,
            z_qubit=cor.test(quant_reading, reads, method="kendall")$statistic) %>%
  dplyr::rename("Kendall correlation between ng/ul of DNA detected following semi-nested PCR and library preparation and reads after QC"=cor)

corr <- full_join(cor_qPCR, cor_qubit)

write.table(corr, "results/kendall_single.tsv", sep ="\t", quote=F, row.names = F)

