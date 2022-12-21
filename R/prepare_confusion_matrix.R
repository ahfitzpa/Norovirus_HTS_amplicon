library(dplyr)
library(tidyr)
library(qiime2R)
library(phyloseq)
library(stringr)
library(Biostrings)

###############################################################################################################################
# read in biom file 
format_biom <- function(feature_table,pipeline){
  biomformat::read_biom(feature_table) %>%
    biomformat::biom_data(.)%>%
    as.matrix(.)%>%
    as_tibble(.,rownames = NA)%>%
    tibble::rownames_to_column("OTU") %>%
    mutate(OTU= paste0(OTU, ".", pipeline))%>%
    rename_at(vars(-OTU),function(x) paste0(x,".", pipeline)) %>% 
    janitor::adorn_percentages("col", na.rm=F) %>%
    mutate_if(is.numeric , replace_na, replace = 0)%>%
    mutate_if(., is.numeric, ~ . * 100)
}

files <- list.files(path ="raw_data",pattern="*.biom", full.names = TRUE,recursive = TRUE)
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
format_expected <- function(table, experiment){
  read.table(table,sep ="\t", header=T)%>%
    rename_at(vars(-taxonomy),function(x) paste0(x,".", experiment)) %>%
    rename_at(vars(-taxonomy),function(x) str_remove(x, "^X"))
}

files <- list.files(path = "raw_data",pattern = "expected_composition.tsv",full.names = TRUE,recursive = TRUE)
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
  mutate(.id= .id/100) %>%
  filter(purrr::reduce(across(where(is.numeric), ~.x!=0), `|`)) %>%
  mutate(source = "expected")%>%
  dplyr::rename(OTU=taxonomy) %>%
  mutate(source="expected") 

#############################################################################################################################################
format_metadata <- function(table, experiment){
  read.table(table,sep ="\t", header=T)%>%
    mutate(SampleID= paste0(SampleID, ".", experiment))
}

files <- list.files(path = "raw_data",pattern = "metadata.tsv",full.names = TRUE,recursive = TRUE)
metadata_list <- list()  
name <- c()
output <- for (i in seq_along(files)) {
  name[i] <- str_replace(files[[i]], "raw_data/", "") %>%
    str_replace(.,"/metadata", "") %>%
    str_replace(., ".tsv", "")
  metadata_list[[i]] <-  format_metadata(files[[i]], name[i]) 
}

metadata <- data.table::rbindlist(metadata_list,fill=T,idcol=T)%>%
  filter(.id!=4)%>%
  mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%
  filter(purrr::reduce(across(where(is.numeric), ~.x!=0), `|`)) %>%
  filter(!grepl("pool", SampleName))%>%
 # filter(grepl("pool", SampleName))%>%
  select(.id,Sample_or_Control,SampleID, DNA_polymerase, RTase)

#############################################################################################################################################
## CLASSIFICATION DATA ##
format_classified <- function(table,pipeline){
  read.table(table, sep ="\t", header=T) %>%
    select(OTU.ID,blast, capsid_type, genotype_capsid)%>%
    mutate(blast=str_remove(blast, "Caliciviridae Norovirus ")) %>%
    dplyr::rename(OTU=OTU.ID,genogroup=blast, genotype=capsid_type)%>%
    mutate(source =pipeline)
}

files <- list.files(path = "raw_data",pattern = "noronet.tsv",full.names = TRUE,recursive = TRUE)
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
full_data <- biom %>%
  separate(OTU, c("OTU", "library"), "\\.") %>% select(-library) %>%
  full_join(., classified, by=c("OTU")) %>%
  select(-OTU, -source) %>%
  separate(taxonomy, c("genogroup", "genotype", "genotype_capsid"), sep=";")%>%
  select(-genogroup, -genotype_capsid) %>%
  pivot_longer(!c("genotype"),names_to="SampleID", values_to="observed") %>%
  select(SampleID, genotype, observed) %>%
  right_join(., metadata, by="SampleID") %>%
  filter(is.na(genotype)!=T)
  
observed <- full_data %>% select(.id,Sample_or_Control,genotype, SampleID,observed) %>%
  distinct()%>%
  filter(observed >0) %>%
  drop_na(genotype) %>%
  mutate(observed=ifelse(observed>0, "present", "absent")) %>%
  select(-Sample_or_Control)%>%
  distinct()
  
expected_df <- expected %>%
  separate(OTU, c("a", "b", "c", "d", "genogroup", "genotype", "genotype_capsid"), sep=";")%>%
  select(-a,-b,-c,-d, -genogroup, -genotype_capsid) %>%
  pivot_longer(!c(".id","source", "genotype"), names_to= "SampleID", values_to="expected")%>%
  right_join(., metadata) %>%
  drop_na(genotype) %>%
  filter(expected >0)%>%
  mutate(expected=ifelse(expected >0, "present", "absent")) %>%
  select(-source,-Sample_or_Control) 
  
observed_expected <- full_join(observed, expected_df) %>%
  right_join(., metadata)  %>% select(-Sample_or_Control) %>%
  drop_na(genotype) %>%
  mutate(source=paste0(DNA_polymerase, "_", RTase))%>%
  mutate(observed=ifelse(is.na(observed)==T, "absent", observed)) %>%
  distinct()

write.table(observed_expected , file = "results/pooled_input_confusion_matrix.tsv", sep="\t", quote = FALSE, row.names = F)
write.table(observed_expected , file = "results/input_confusion_matrix.tsv", sep="\t", quote = FALSE, row.names = F)
