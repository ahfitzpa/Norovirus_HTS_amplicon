#OVERALL UNIFRAC 
library(dplyr)
library(tidyverse)
library(phyloseq)
library(qiime2R)
library(vegan)

#import imported unifrac data from qiime2, performed using long format data
unifrac_unweighted <- list.files(path ="data/unifrac_input", pattern = "unweighted_unifrac_long$", include.dirs=T, full.names=T) %>%
  lapply(., function(x) read.table(paste0(x,"/distance-matrix.tsv"),  sep="\t", header=T) %>% tibble::column_to_rownames("X")) %>%
  setNames(c("1", "2", "3"))%>%
  lapply(.,function(x) x %>% rename_with(~str_remove(., 'X')))

meta <- list.files(path ="data/unifrac_input", pattern="metadata", full.names=T) %>%
  lapply(., function(x) read.table(x,  sep="\t", header=T))%>%
  setNames(c("1", "2", "3"))

pc.uuf <- lapply(unifrac_unweighted, cmdscale, k=2, eig=TRUE)

funcPoints1 <- lapply(pc.uuf ,function(x) data.frame(x$points))
eigs <- lapply(pc.uuf , eigenvals)
eigs <- lapply(eigs, function(x) x/sum(x))
eigs[1:2]

keys <- unique(c(names(funcPoints1), names(meta)))
funcPoints2 <- lapply(setNames(keys, keys), function(x) {cbind(meta[[x]], funcPoints1[[x]])}) %>%
  lapply(., setNames, c('sampleid', 'source', 'DNA_polymerase', 'RTase', 'PC1', 'PC2' )) 
funcPoints <- funcPoints2 %>% data.table::rbindlist(.,  idcol=T) 
library(ggplot2)
library(ggforce)
library(ggthemes)
library(ggpubr)
pal_designer <- c("AmpliTaq_Gold"="#DB073D","Kapa_HiFi"="#DBA507","Kapa_Robust"="#43A047",
                  "expected"="#07485B")
RTase_list <-  c("expected", "SSII")
funcPoints_DNA <- funcPoints %>% filter(RTase %in% RTase_list)
ggplot(funcPoints_DNA, aes(x = PC1, y = PC2,fill = DNA_polymerase)) +
  geom_mark_ellipse(aes(fill = DNA_polymerase,label=DNA_polymerase)) +
  geom_point(size = 3, pch = 21)+
  scale_fill_manual(values =pal_designer)+
  theme_classic()+
  xlim(-0.7, 0.7)+ylim(-0.6, 0.7)+
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle =45), text = element_text(size = 20)) 
ggplot2::ggsave("plots/unifrac_DNApol.png") 

pal_RTase<- c("Luna"="#E67E25", "SSII"="#E74C3B", "SSIV"="#23B99A")
DNA_list<-  c("expected", "AmpliTaq_Gold")
funcPoints_RTase <- funcPoints %>% filter(DNA_polymerase %in% DNA_list)
ggplot(funcPoints_RTase, aes(x = PC1, y = PC2,fill = RTase)) +
  geom_mark_ellipse(aes(fill = RTase,label=RTase)) +
  geom_point(size = 3, pch = 21)+
  scale_fill_manual(values =pal_RTase)+
  theme_classic()+
  xlim(-0.7, 0.8)+ylim(-0.5, 0.8)+
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle =45), text = element_text(size = 20)) 
ggplot2::ggsave("plots/unifrac_RTase.png")


############################################################################################################################################################
meta_DNA <- list.files(path ="data/unifrac_input", pattern="metadata", full.names=T) %>%
  lapply(., function(x) read.table(x,  sep="\t", header=T))%>%
  setNames(c("1", "2", "3")) %>%
  purrr::map(~ .x %>% filter(RTase %in% RTase_list))
unifrac_DNA <- mapply(function(x,y) x %>% select(y$Sample_ID) %>% filter(row.names(x) %in% y$Sample_ID),unifrac_unweighted, meta_DNA, SIMPLIFY =F)
AnT_DNA <- mapply(function(x,y) anosim(x, y$DNA_polymerase,permutations = 999, distance = "jaccard"), unifrac_DNA,meta_DNA, SIMPLIFY =F)
AnT_df_DNA <- lapply(AnT_DNA, function(x) data.frame(R_statistic=x$statistic, p_significance=x$signif)) %>%
  data.table::rbindlist(.,  idcol=T)

meta_RTase <- list.files(path ="data/unifrac_input", pattern="metadata", full.names=T) %>%
  lapply(., function(x) read.table(x,  sep="\t", header=T))%>%
  setNames(c("1", "2", "3")) %>%
  purrr::map(~ .x %>% filter(DNA_polymerase %in% DNA_list))
unifrac_RTase <- mapply(function(x,y) x %>% select(y$Sample_ID) %>% filter(row.names(x) %in% y$Sample_ID),unifrac_unweighted, meta_RTase, SIMPLIFY =F)
AnT_RTase <- mapply(function(x,y) anosim(x, y$RTase,permutations = 999, distance = "jaccard"), unifrac_RTase, meta_RTase, SIMPLIFY =F)
AnT_df_RTase <- lapply(AnT_RTase, function(x) data.frame(R_statistic=x$statistic, p_significance=x$signif)) %>%
  data.table::rbindlist(.,  idcol=T)

permanova_DNA <- mapply(function(x,y) adonis2(x ~y$DNA_polymerase,permutations = 999, distance = "jaccard"), unifrac_DNA, meta_DNA, SIMPLIFY=F)
permanova_df_DNA <- lapply(permanova_DNA, function(x) data.frame(df=x[1,1], sumofsquares=x[1,2], R_sq=x[1,3],F_stat = x[1,4], pr_f=x[1,5]))%>%
  data.table::rbindlist(.,  idcol=T)

permanova_RTase <- mapply(function(x,y) adonis2(x ~y$RTase,permutations = 999, distance = "jaccard"), unifrac_RTase, meta_RTase, SIMPLIFY=F)
permanova_df_RTase <- lapply(permanova_RTase, function(x) data.frame(df=x[1,1], sumofsquares=x[1,2], R_sq=x[1,3],F_stat = x[1,4], pr_f=x[1,5]))%>%
  data.table::rbindlist(.,  idcol=T)

# Write Output ---
write.table(permanova_df_DNA, "results/stats_abundance_adonis2_unifrac_DNA.tsv", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(permanova_df_RTase, "results/stats_abundance_adonis2_unifrac_RTase.tsv", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

######## Post-hoc Test for Permanova ##########
library(pairwiseAdonis)
library(compositions)
permtst_DNA<- mapply(function(x,y) pairwise.adonis(x, y$DNA_polymerase), unifrac_DNA, meta_DNA, SIMPLIFY =F)
permtst_RTase <- mapply(function(x,y) pairwise.adonis(x, y$RTase), unifrac_RTase, meta_RTase, SIMPLIFY =F)

df_padjust <- permtst_RTase %>%
  data.table::rbindlist(.,  idcol=T) %>%
  filter_all(any_vars(str_detect(., pattern = "expected"))) %>%
  tidyr::separate(pairs, c("pipeline1","vs", "pipeline2"), sep=" ") %>%
  select(-vs, -pipeline2, -Df, -F.Model,-sig, -SumsOfSqs, -p.value) 

R2 <- df_padjust %>% select(-p.adjusted) %>% tidyr::pivot_wider(names_from=pipeline1, values_from=R2)
df_padjust %>% select(-R2) %>% tidyr::pivot_wider(names_from=pipeline1, values_from=p.adjusted)

write.table(R2, file = "results/posthoc_unifrac_rtase.tsv", sep="\t", quote = FALSE, row.names = F)

