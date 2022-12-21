#plot results of each sequencng library

library(dplyr)
library(tidyr)
library(phyloseq)
library(stringr)
library(Biostrings)

biom <- read.table("data/biom.tsv",  sep ="\t", header=T) 
classified <- read.table("data/classified.tsv", sep="\t", header=T)
metadata <- read.table("data/metadata.tsv", sep="\t", header=T)

check <- biom %>%   
  separate(OTU, c("OTU", "library"), "\\.")  %>%
  pivot_longer(!c(OTU, .id, library),names_to="SampleID", values_to="abundance") %>%
  filter(abundance >1) %>%
  left_join(., classified, by=c(".id", "OTU")) %>%
  mutate(SampleID=str_replace(SampleID, "^X", "")) %>%
  filter(!grepl('NEG', SampleID)) %>%
  separate(taxonomy, c("genogroup", "genotype", "genotype_capsid"), ";") %>%
  left_join(., metadata, by="SampleID") %>%
  unite(enzyme, c("RTase","DNA_polymerase"))%>%
  group_by(enzyme,library, genotype_capsid) %>%
  summarise() %>%
  mutate(bars_by_foo = paste0(genotype_capsid, collapse = ",")) %>%
  ungroup() %>%
  select(-genotype_capsid) %>%
  distinct()

write.table(check, "data/genotypes_enzyme.tsv", sep ="\t", quote=F, row.names = F)
library(RColorBrewer)
library(ggplot2)
# Define the number of colors you want
nb.cols <- 14
mycolors <-colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)

ggplot(check,mapping= aes(x=SampleID, y=abundance,fill=genotype_capsid)) +
  geom_bar(stat="identity", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()+
  facet_wrap(~ .id, scales="free") +
  scale_fill_manual(values= mycolors)
