library(dplyr)
library(tidyr)
library(stringr)
library(magrittr)

fasta_table <- function(fasta){
  fastaFile <- Biostrings::readDNAStringSet(fasta, format = "fasta",use.names = TRUE)
  OTU <- names(fastaFile)
  sequence <-  paste(fastaFile)
  df <- data.frame(OTU, sequence)
  return(df)
}

otu <- read.table("data/biom_observed_expected.tsv", header=T, sep="\t") 
metadata <- read.table("data/metadata.tsv", header=T, sep="\t") 
#fasta <- fasta_table("data/expected_observed.fasta") 

###############################################################################################################################
data_cluster <- otu %>%
  filter(.id !='4') %>%
  pivot_longer(!c(.id, source,taxonomy),names_to="SampleID", values_to="abundance") %>%
  mutate(SampleID=str_remove(SampleID, "^X"))%>%
  left_join(., metadata) %>%
  #filter(RTase=="SSII") %>%
  filter(DNA_polymerase=="AmpliTaq_Gold") %>%
  #mutate(source=ifelse(source =="observed", paste0(RTase), source)) %>%
  mutate(source=ifelse(source =="observed", paste0(DNA_polymerase, "_", RTase), source)) %>%
  mutate(DNA_polymerase=ifelse(source =="expected", paste0("expected"), DNA_polymerase)) %>%
  mutate(RTase=ifelse(source =="expected", paste0("expected"), RTase)) %>%
  mutate(abundance=ifelse(abundance >0, 1,0)) %>%
  select(.id, taxonomy,SampleID, abundance,source, DNA_polymerase, RTase) %>%
  mutate(row = row_number()) %>%
  pivot_wider(names_from="SampleID", values_from="abundance") %>%
  select(-row)%>%
  distinct() %>%
  group_by(.id,taxonomy) %>%
  mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%
  filter(purrr::reduce(across(where(is.numeric), ~.x!=0), `|`))%>%
  ungroup() 
  
##########################################################################################################
## Perform GROUPED vegdist bray 
library(vegan)
df <- data_cluster %>% 
  distinct() %>%
  #tibble::column_to_rownames("OTU_sample") %>%
  select(-2,-3,-4,-5) %>%
  na.omit() %>%
  group_split(.id, .keep =F) %>%
  lapply(.,  select, where(~ any(. != 0))) %>%
  setNames(sort(unique(data_cluster$.id)))

meta <- data_cluster %>%
  na.omit() %>%
  select(1,2,3,4,5) %>%
  #tibble::column_to_rownames("OTU_sample") %>%
  group_split(.id,.keep =F) %>%
  lapply(.,  select, where(~ any(. != 0)))%>%
  setNames(sort(unique(data_cluster$.id)))

complete_vegdist<- lapply(df, vegdist, "jaccard", binary=T)

funcPcoa <-lapply(complete_vegdist,cmdscale, k=2, eig=TRUE)
funcPoints1 <- lapply(funcPcoa,function(x) data.frame(x$points))
eigs <- lapply(funcPcoa, eigenvals)
eigs <- lapply(eigs, function(x) x/sum(x))
eigs[1:2]

keys <- unique(c(names(funcPoints1), names(meta)))
funcPoints2 <- lapply(setNames(keys, keys), function(x) {cbind(meta[[x]], funcPoints1[[x]])})  %>%
 #lapply(., setNames, c('taxonomy', 'source','enzyme','PC1', 'PC2' )) 
  lapply(., setNames, c('taxonomy', 'source','DNA_polymerase','RTase','PC1', 'PC2' )) 
funcPoints <- funcPoints2 %>% data.table::rbindlist(.,  idcol=T, fill=T)

library(ggplot2)
library(ggforce)
library(ggthemes)
library(ggpubr)
# box plot based 1 dimensional scaling
pal_designer <- c("AmpliTaq_Gold"="#DB073D","Kapa_HiFi"="#DBA507","Kapa_Robust"="#43A047",
"expected"="#07485B")
pal_RTase =c ("Luna"="#E67E25", "SSII"="#E74C3B", "SSIV"="#23B99A")

ggplot(funcPoints, aes(x = PC1, y = PC2,fill = RTase)) +
  geom_mark_ellipse(aes(fill = RTase,label=RTase)) +
  geom_point(size = 3, pch = 21)+
  scale_fill_manual(values =pal_RTase)+
  theme_classic()+
  xlim(-0.7, 0.7)+ylim(-0.7, 0.7)+
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle =90), text = element_text(size = 30)) 
ggplot2::ggsave("plots/jaccard_RTase.png") 

##############################################################################
AnT <- mapply(function(x,y) anosim(x, y$RTase,permutations = 999, distance = "jaccard"), complete_vegdist, meta, SIMPLIFY =F)
AnT_df <- lapply(AnT, function(x) data.frame(R_statistic=x$statistic, p_significance=x$signif)) %>%
  data.table::rbindlist(.,  idcol=T)

y_permanova <- mapply(function(x,y) adonis2(x ~y$RTase,permutations = 999, distance = "jaccard"), complete_vegdist, meta, SIMPLIFY=F)
y_permanova_df <- lapply(y_permanova, function(x) data.frame(df=x[1,1], sumofsquares=x[1,2], R_sq=x[1,3],F_stat = x[1,4], pr_f=x[1,5]))%>%
  data.table::rbindlist(.,  idcol=T)
# Write Output ---
write.table(AnT_df, "results/jaccard_DNApol_anosim.tsv", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(y_permanova_df, "results/jaccard_rtase.tsv", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

######## Post-hoc Test for Permanova ##########
permtst <- mapply(function(x,y) RVAideMemoire::pairwise.perm.manova(resp=x ,fact=y$RTase,  test = "Pillai", nperm =999, progress = TRUE, p.method = "none"), complete_vegdist, meta, SIMPLIFY =F)

df <- lapply(permtst,function(x) reshape2::melt(x$p.value)) %>%
  lapply(., setNames,  c("source", "enzyme", "pvalue"))

df_padjust <- mapply(function(x) x %>% dplyr::mutate(p.adjust= p.adjust(x$pvalue, method = "bonferroni", n = length(x$pvalue))), df, SIMPLIFY =F)%>%
  data.table::rbindlist(.,  idcol=T) %>%
  filter_all(any_vars(str_detect(., pattern = "expected"))) %>%
  unite(enzyme, c("source", "enzyme")) %>%
  mutate(enzyme= gsub('expected_', '', enzyme))%>%
  mutate(enzyme= gsub('_expected', '', enzyme))%>%
  select(-pvalue) %>%
  filter(is.na(p.adjust)==F) %>%
  pivot_wider(names_from=enzyme, values_from=p.adjust)

write.table(df_padjust, file = "results/posthoc_jaccard_rtase.tsv", sep="\t", quote = FALSE, row.names = F)
