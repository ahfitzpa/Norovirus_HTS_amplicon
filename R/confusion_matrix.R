library(dplyr)
library(tidyr)
library(qiime2R)
library(phyloseq)
library(stringr)
library(Biostrings)
########################################################################################################################################
# import data
observed_expected <- read.table( "results/input_confusion_matrix.tsv", sep="\t", header=T)%>%
  select(.id,SampleID,genotype,observed, expected, source) %>%
  mutate(source=factor(source, levels=c("AmpliTaq_Gold_SSII","Kapa_HiFi_SSII","Kapa_Robust_SSII","AmpliTaq_Gold_Luna",
 "AmpliTaq_Gold_SSIV","Kapa_HiFi_Luna","Kapa_HiFi_SSIV","Kapa_Robust_Luna","Kapa_Robust_SSIV" ))) %>%
  mutate(expected=factor(expected, levels=c("present", "absent"))) %>%
  mutate(observed=factor(observed, levels=c("present", "absent")))%>%
  unite(lib_id, c(".id","source"))%>% 
  mutate(lib_id=as.factor(lib_id))%>%
  group_by(lib_id) %>% 
  group_split() %>% 
  setNames(sort(unique(observed_expected$lib_id)))
  
#########################################################################################################################################
`%nin%` = Negate(`%in%`)
confusion_table <- lapply(observed_expected , function(df) {
  df  %>%
    select(-SampleID) %>%
    add_row(expected="absent", observed="absent") %>%
    mutate(expected=factor(expected,levels=c("present", "absent"))) %>%
    mutate(observed=factor(observed,levels=c("present", "absent")))
}
)

library(yardstick)
cm <-lapply(confusion_table, yardstick::conf_mat, truth =expected,estimate =observed)

metrics_cm <- lapply(cm, function(df) {df %>%summary() %>%
    select(-.estimator) %>%
    filter(.metric %in%
             c("spec", "sens", "accuracy", "precision", "recall", "f_meas"))})


library(kableExtra)
output_cm <- metrics_cm %>% 
  dplyr::bind_rows(., .id = "lib_id")%>%
  mutate_if(is.numeric, round, digits = 2) %>% 
  tidyr::pivot_wider(names_from= .metric, values_from= .estimate) 

results <- output_cm  %>% arrange(desc(f_meas))%>%
  kbl(caption = "Accuracy of genotypic characterisation using various RTases and DNA polymerases based on RDP taxonomy") %>%
  kable_classic(full_width = F, html_font = "Cambria")
results
write.table(output_cm, file = "results/pooled_confusion_matrix.tsv", sep="\t", quote = FALSE, row.names = F)

