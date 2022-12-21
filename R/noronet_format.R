noronet_files <- list.files(path= "raw_data/LP2/noronet_files/", patter="^result", full.names= T)
#noronet_files <- list.files(path= "raw_data/LP6/", patter="^result", full.names= T)

for (file in noronet_files){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    dataset <- read.csv(file, header=TRUE)
  }
  # if the merged dataset does exist, append to it
     if (exists("dataset")){
         temp_dataset <-read.csv(file, header=TRUE)
         dataset<-rbind(dataset, temp_dataset)
         rm(temp_dataset)
       }
  }
  
#dataset <-  read.csv("raw_data/LP1/results_noronet.csv", header=TRUE)
noronet_data <- dataset %>%
  janitor::clean_names() %>%
  select(name,blast, capsid_type, capsid_subtype,  blast_score) %>%
  mutate(capsid_subtype = str_replace(capsid_subtype, "Could not assign", "")) %>%
  mutate(capsid_subtype = str_replace(capsid_subtype, " ", "Could not assign")) %>%
  mutate(capsid_type=str_trim(capsid_type, "right")) %>%
  unite(genotype_capsid, c("capsid_type", "capsid_subtype"), sep = " ", na.rm = TRUE,remove=F) %>%
  mutate(genotype_capsid=str_trim(genotype_capsid, "right")) %>%
  mutate(genotype_capsid= str_replace(genotype_capsid, "GII.4_", "GII.4")) %>% 
  mutate(capsid_type = case_when(capsid_type == "" ~ "Could not assign", 
                            TRUE ~as.character(capsid_type ))) %>%
  dplyr::rename(OTU.ID=name) %>%
  distinct()

write.table(noronet_data, "raw_data/LP2/noronet.tsv", sep="\t", quote = FALSE, row.names = F )
