library(dplyr)
library(tidyr)
library(qiime2R)
library(phyloseq)
library(stringr)
library(Biostrings)
# combine fasta and table output for all pipelines
################################################################################################################################
fasta_table <- function(fasta){
  fastaFile <- readDNAStringSet(fasta, format = "fasta",use.names = TRUE)
  OTU <- names(fastaFile)
  sequence <-  paste(fastaFile)
  df <- data.frame(OTU, sequence)
  return(df)
}

#functions
writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

################################################################################################################################
all_sequences <- data.table::rbindlist(mapply(
    c,
    (
      list.files(
        path = "raw_data",
        pattern = "*.fasta",
        full.names = TRUE,
        recursive = TRUE
      ) %>%
        lapply(
          fasta_table
        )
    ),
    (
      list.files(
        path = "raw_data",
        pattern = "*.fasta",
        full.names = TRUE,
        recursive = TRUE
      ) %>%
        #basename() %>%
        as.list()
    ),
    SIMPLIFY = FALSE
  ),
  fill = T)

otus <- all_sequences %>% 
  separate(V1, into=c("folder", "run", "file_type"), sep= "/") %>%
  #unite(OTU_id, c("OTU", "run"),sep= ".") %>%
  mutate(source = ifelse(grepl("expected", file_type), "expected", "observed")) %>%
  #mutate(OTU_id=ifelse(grepl("expected", source), paste0(OTU_id, ".", "expected"), OTU_id)) %>%
  dplyr::rename(name=OTU, seq=sequence)

writeFasta(otus, "data/expected_observed.fasta")
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
biom <- data.table::rbindlist(biom_list,fill=T,idcol=T) %>%
  mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%
  filter(purrr::reduce(across(where(is.numeric), ~.x!=0), `|`)) 
write.table(biom, "data/biom.tsv", sep ="\t", quote=F, row.names = F)

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
  filter(purrr::reduce(across(where(is.numeric), ~.x!=0), `|`)) %>%
  mutate(.id=.id/100)%>%
  mutate(source = "expected")%>%
  dplyr::mutate(OTU=taxonomy) %>%
  mutate(source="expected")%>%
  mutate(taxonomy=strex::str_after_nth(taxonomy, ";", 4))
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
  mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%
  filter(purrr::reduce(across(where(is.numeric), ~.x!=0), `|`)) %>%
  select(.id,SampleID, DNA_polymerase, RTase, EXP_ID, quant_reading, qPCR, Sample_or_Control, season, genogroup)
write.table(metadata, "data/metadata.tsv", sep ="\t", quote=F, row.names = F)

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
write.table(classified, "data/classified.tsv", sep ="\t", quote=F, row.names = F)
#############################################################################################################################################
biom_observed_expected <- biom %>%
  separate(OTU, c("OTU", "library"), "\\.") %>% select(-library, -.id) %>%
  pivot_longer(!OTU,names_to="SampleID", values_to="abundance") %>%
  filter(abundance >1) %>%
  mutate(Session = row_number()) %>%
  pivot_wider(names_from="SampleID", values_from="abundance") %>%
  select(-Session) %>%
  full_join(., classified, by=c("OTU")) %>%
  full_join(., expected) %>%
  select(.id, OTU,taxonomy,source,everything()) %>%
  mutate(taxonomy=ifelse(is.na(taxonomy)==T, OTU, taxonomy)) %>%
  select(-OTU) %>%
  mutate_if(is.numeric , replace_na, replace = 0) %>%
  filter(purrr::reduce(across(where(is.numeric), ~.x!=0), `|`)) 

write.table(biom_observed_expected, "data/biom_observed_expected.tsv", sep ="\t", quote=F, row.names = F)

#############################################################################################################################################
# data for unifrac
dir.create("data/unifrac_input")

unifrac_otus <- otus %>% select(name,run, seq) %>%
  dplyr::rename(OTU=name) %>%
  mutate(run=str_replace(run, 'LP',''))%>%
  mutate(run=str_replace(run, '4','3'))%>%
  mutate(run=str_replace(run, '6','4'))%>%
  mutate(run=as.numeric(run))%>%
  dplyr::rename(.id=run)

unifrac <- biom_observed_expected %>%
  filter(.id != '4') %>%
  left_join(., unifrac_otus)%>%
  pivot_longer(!c(.id,OTU,seq,source,taxonomy),names_to="SampleID", values_to="abundance") %>%
  filter(abundance >1) %>%
  distinct()%>%
  left_join(., metadata) %>%
  filter(Sample_or_Control=="True Sample") %>%
  select(-EXP_ID, -quant_reading, -qPCR, -Sample_or_Control, -season, -genogroup,-Technical_replicate)%>%
  mutate(source=ifelse(source =="observed", paste0(DNA_polymerase, ".", RTase), source)) %>%
  mutate(DNA_polymerase=ifelse(source =="expected", paste0("expected"), DNA_polymerase)) %>%
  mutate(RTase=ifelse(source =="expected", paste0("expected"), RTase))%>%
  unite(OTU_ID, c("OTU",".id"),sep=".", remove =F) %>%
  unite(Sample_ID, c("SampleID", "source"),sep=".", remove =F) %>%
  select(-SampleID)%>% distinct() 

unifrac_biom <- unifrac %>%
  select(-taxonomy,-OTU, -seq,-source, -DNA_polymerase, -RTase)%>%
  pivot_wider(names_from=Sample_ID,values_from=abundance) %>%
  mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%
  filter(purrr::reduce(across(where(is.numeric), ~.x!=0), `|`))

write.table(unifrac_biom, "data/unifrac_input/unifrac_biom.tsv",col.names=T,row.names=F, sep="\t", quote=F) 

unifrac_biom <- unifrac %>%
  select(-taxonomy,-OTU, -seq,-source, -DNA_polymerase, -RTase)%>%
  group_split(.id, .keep=F) %>%
  setNames(sort(unique(unifrac$.id)))%>%
  purrr::map(~ .x %>% pivot_wider(names_from=Sample_ID,values_from=abundance)) %>%
  purrr::map(~ .x %>% mutate(across(where(is.numeric), ~replace(., is.na(.), 0)))) %>%
  purrr::map(~ .x %>% filter(purrr::reduce(across(where(is.numeric), ~.x!=0), `|`)))

for(i in 1:length(unifrac_biom)) {                              # Head of for-loop
  write.table(unifrac_biom[[i]],                              # Write CSV files to folder
              paste0("data/unifrac_input/",names(unifrac_biom[i]), "_unifrac_biom.tsv"),
              col.names=T,row.names=F, sep="\t", quote=F)
}

unifrac_metadata <-  unifrac %>%
  select(.id, Sample_ID,source, DNA_polymerase, RTase)%>%
  group_split(.id, .keep=F) %>%
  setNames(sort(unique(unifrac$.id)))%>%
  purrr::map( ~ .x %>% distinct(.)) 

write.table(unifrac_metadata, "data/unifrac_input/unifrac_long_metadata.tsv",col.names=T,row.names=F, sep="\t", quote=F)  
for(i in 1:length(unifrac_metadata)) {                              # Head of for-loop
  write.table(unifrac_metadata[[i]],                              # Write CSV files to folder
              paste0("data/unifrac_input/",names(unifrac_metadata[i]), "_long_metadata.tsv"),
              col.names=T,row.names=F, sep="\t", quote=F)
}

unifrac_seq <- unifrac %>% select(.id, OTU_ID, seq)%>%
  dplyr::rename(name=OTU_ID)%>%
  group_split(.id, .keep=F) %>%
  setNames(sort(unique(unifrac$.id)))%>%
  lapply(., distinct)

writeFasta(unifrac_seq,"data/unifrac_input/unifrac_expected_observed.fasta")
for(i in 1:length(unifrac_seq)) {                              # Head of for-loop
  writeFasta(unifrac_seq[[i]],                              # Write CSV files to folder
             paste0("data/unifrac_input/",names(unifrac_seq[i]), "_", "expected_observed.fasta"))
}
