# split fasta file for making batchdb for expected sequences
# need fasta file for each sim/pipelien
library(readxl)
library(dplyr)
library(magrittr)
library(tidyr)
library(stringr)

fasta_table <- function(fasta){
  fastaFile <- Biostrings::readDNAStringSet(fasta, format = "fasta",use.names = TRUE)
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

metadata <- read.table("data/unifrac_input/unifrac_long_metadata.tsv", header=T, sep="\t") 
sequences <- fasta_table("data/unifrac_input/unifrac_expected_observed.fasta")
data <- read.table("data/unifrac_input/unifrac_biom.tsv", header=T, sep="\t") 
pipelines<- c("dada2", "deblur", "vsearch", "unoise3", "FROGS")
data_OTU <- data %>% 
  dplyr::rename(OTU=OTU_ID)%>%
  pivot_longer(!c(".id", "OTU"), names_to="Sample_ID", values_to="abundance") %>%
  mutate(Sample_ID=str_remove(Sample_ID, "^X")) %>%
  left_join(., metadata )%>%
  left_join(., sequences)  %>% 
  mutate(OTU_ID= strex::str_before_last(OTU, "\\.")) %>%
  mutate(lib= strex::str_after_last(OTU, "\\."))%>% 
  select(-OTU) %>%
  filter(abundance >0) %>%
  unite(source_id, c("source",".id"), remove=F) %>%
  unite(OTU, c("source",".id","OTU_ID")) %>%
  dplyr::rename(name=OTU, seq=sequence) %>%
  split(f = as.factor(.$source_id))

#############################################################################################################################################
# data for blastdb
dir.create("data/blastdb")
list_filenames <- paste0("data/blastdb" , paste0(names(data_OTU)), ".fasta")
sapply(names(data_OTU), 
       function (x) writeFasta(data_OTU[[x]], filename=paste0("data/blastdb/",x, ".fasta")))
