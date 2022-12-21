#load qiime2 output into R 
library(dplyr)
library(tidyverse)
library(Biostrings)

expected_data <- read.delim(file = "expected_composition2.tsv", sep = "\t", header =T,check.names=F)

output1 <- expected_data  %>%
  tidyr::separate(taxonomy, c("taxonomy","library"), " ") 

library <- output1$library[[1]]
output <- output1 %>% select(-library) %>%
        dplyr::rename_at(vars(-taxonomy), ~paste0(., ".", library))
        
write.table(output,"expected_composition.tsv", sep= "\t", row.names=F, col.names=T, quote=F)
