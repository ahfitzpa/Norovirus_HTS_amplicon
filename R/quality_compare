library(dplyr)
library(tidyr)
library(stringr)
library(Biostrings)
# Prepare quality sequencing data
###############################################################################################################################
# read in quality control file 
format_quality <- function(table,experiment){
  read.table(table,sep ="\t", header=T)%>%
    mutate(ID= paste0(experiment))
}

files <- list.files(path ="raw_data",pattern="*.stats", full.names = TRUE,recursive = TRUE)
quality_list <- list()  
name <- c()
output <- for (i in seq_along(files)) {
  name[i] <- str_replace(files[[i]], "raw_data/", "") %>%
    str_replace(.,"3.quality_stats/", "") %>%
    str_replace(., ".stats", "")
  quality_list[[i]] <-  format_quality(files[[i]], name[i]) 
}

quality <- data.table::rbindlist(quality_list,fill=T)%>%
  separate(ID, c("library", "SampleID"), sep="/")

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
  filter(Sample_or_Control=="True Sample") %>%
  select(.id,SampleID, DNA_polymerase, RTase, EXP_ID, quant_reading, qPCR) %>%
  unite(source, c("RTase", "DNA_polymerase"), sep="_", remove=F) 

##############################################################################################################################
library(tidyverse)
library(ggpubr)
library(rstatix)
quality_meta <- right_join(quality, metadata) %>%
  filter(library=="LP4") %>%
  mutate(Mean_Q =log10(Mean_Q+0.001)) %>%
  mutate(Mean_EE =log10(Mean_EE+0.001)) %>%
  mutate(Mean_Pe =log10(Mean_Pe+0.001))

quality_meta %>%
  group_by(DNA_polymerase) %>%
  rstatix::kruskal_test(RTase ~ Mean_Q)

quality_meta %>%
  group_by(RTase) %>%
  rstatix::kruskal_test(DNA_polymerase ~ Mean_Q)

kruskal.test(quality_meta$source ~ quality_meta$Mean_Q)
kruskal.test(quality_meta$DNA_polymerase~ quality_meta$Mean_EE)
kruskal.test(quality_meta$RTase ~ quality_meta$Mean_Q)

quality_meta %>% 
  summarize(cor=stats::cor.test(source,Mean_Q method="kendall")$estimate, 
            p.value_qPCR=cor.test(qPCR, reads, method="kendall")$p.value,
            z=cor.test(qPCR, reads, method="kendall")$statistic)

MeanQ_posthoc <- rstatix::dunn_test(quality_meta, Mean_Q~ source, p.adjust.method = "bonferroni", detailed =T) %>%
  filter(p.adj <0.05)
MedEE_posthoc <- rstatix::dunn_test(quality_meta, Med_EE~ source, p.adjust.method = "bonferroni", detailed =T) %>%
  filter(p.adj <0.05)
Mean_Pe_posthoc <- rstatix::dunn_test(quality_meta, Med_Pe~ source, p.adjust.method = "bonferroni", detailed =T) %>%
  filter(p.adj <0.05)

library(ggdist)
library(ggthemes)
library(tidyquant)
library(ggpubr)
enzyme_pal <- c("Luna_AmplitaqGold"="#FFC142", "SSII_AmplitaqGold"="#5793FF", "SSIV_AmplitaqGold"= "#1f9e61",
                "SSIV_KapaRobust"="#FFA15C", "SSII_KapaRobust"="#2d9685","Luna_KapaRobust"="#C100FF",
                "SSII_KapaHiFi"="#F24162", "Luna_KapaHiFi"="#400040", "SSIV_KapaHiFi"="#E83938")

ggplot(quality_meta, aes(x=reorder(source,-Mean_Q, na.rm=T), y =Mean_Q, fill = reorder(source,-Mean_Q, na.rm=T))) +
  scale_fill_manual(values=enzyme_pal)+
  scale_colour_manual(values=enzyme_pal)+
  geom_violin()+
  geom_boxplot(colour="white",width=.1)+
  labs(x="RTase and DNA polymerase",
       y = "Mean Q score"
  )  +
  stat_summary(aes(color="Black"),fun.y=median, geom="point", alpha=1,
                   size=3.2, shape = 21, fill = "black", stroke = 1) +
  theme_classic(base_size=30)+
  theme(axis.text.x = element_text(angle =35, hjust=1, size=20))+
  guides(color = "none", fill="none")
ggsave("plots/QC_violin_meanee.png", dpi=600)
