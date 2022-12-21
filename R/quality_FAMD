# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # #  K mode Cluster Analysis # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
library(cluster)    # clustering algorithms
library(Rtsne)
library(tibble)
library(FactoMineR)
library(factoextra) # clustering algorithms & visualization
# https://uc-r.github.io/kmeans_clustering
# https://dabblingwithdata.wordpress.com/2016/10/10/clustering-categorical-data-with-r/

clean <-  quality_meta  %>%
  mutate(row_id=paste0(Pos, "_", SampleID, "_", source)) %>%
  dplyr::select(row_id,DNA_polymerase, RTase, Mean_Q) %>%
  #select( -source, -Pos, -Recs, -PctRecs, -DNA_polymerase, -RTase, -EXP_ID,-.id, -library, -SampleID, -qPCR, -quant_reading)%>%
  column_to_rownames(., var = "row_id") %>%
  distinct() %>%
  mutate(source=as.factor(source)) %>%
  mutate(DNA_polymerase=as.factor(DNA_polymerase)) %>%
  mutate(RTase=as.factor(RTase)) 
 # mutate(SampleID=as.factor(SampleID))

# k-mode clustering, each cluster is represented by its center (i.e, centroid) 
# which corresponds to the mean of points assigned to the cluster
# Compute FAMD
# Factor analysis for mixed data, PCA for quantitative data and MCA for qualtiative data
# Tales incto account both quantitative and qualtiative variables
# https://rpkgs.datanovia.com/factoextra/reference/fviz_famd.html

# Compute Gower distance
df <- clean 
gower_dist <- daisy(df, metric = "gower")
gower_mat <- as.matrix(gower_dist)
tsne_obj <- Rtsne(gower_dist, is_distance = TRUE)

# Most similar and dissimilar clients according to Gower distance:
# Print most similar clients
df[which(gower_mat == min(gower_mat[gower_mat clean <-  quality_meta  %>%
                                      mutate(row_id=paste0(Pos, "_", SampleID, "_", source)) %>%
                                      #dplyr::select(row_id, Mean_Q, Mean_EE, Mean_EE) %>%
                                      select( -source, -Pos, -Recs, -PctRecs, -DNA_polymerase, -RTase, -EXP_ID,-.id, -library, -SampleID, -qPCR, -quant_reading)%>%
                                      column_to_rownames(., var = "row_id") %>%
                                      distinct() != min(gower_mat)]), arr.ind = TRUE)[1, ], ]
# Print most dissimilar clients
df[which(gower_mat == max(gower_mat[gower_mat != max(gower_mat)]), arr.ind = TRUE)[1, ], ]
set.seed(123)

# function to compute total within-cluster sum of square 
wss <- function(k) {
  kmeans(df, k, nstart = 10 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:15

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)
plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

gap_stat <- clusGap(df, FUN = kmeans, nstart = 25,
                    K.max = 10, B = 50)
fviz_gap_stat(gap_stat)


final <- kmeans(df, 4, nstart = 25)
fviz_cluster(final, data = df) 


ggplot(k9_pam_fit, aes(x = X, y = Y, color = cluster, fill = cluster))+
  geom_point(size = 4, shape = 21)+
  theme_bw()+
  ggsci::scale_color_nejm()+
  ggsci::scale_fill_nejm()

check <- df %>%
  as_tibble() %>%
  mutate(cluster = k2$cluster,
         row_id = row.names(clean)) %>%
  left_join(., clean)
  ggplot(aes(UrbanPop, Murder, color = factor(cluster), label = state)) +
  geom_text()

#save cluster of 7 as recommnded from checking optimum clustering
kmode8 <- visualise_tsne_obj(k9_pam_fit)
kmodeoutput<-  readxl::read_excel('results/kmode-edited.xlsx')
tablekmode <-   ggtexttable(kmodeoutput, rows = NULL, theme =  ttheme('classic'))
ggarrange(tablekmode, kmode7,ncol=1, nrow = 2, heights = c(0.2, 0.8))
ggsave('figures/kmode.png')

name_vector  <- c('c1','c2', 'c3', 'c4', 'c5', 'c6', 'c7')
kmoderesults  <-results_pam(k7_pam_fit)
names(kmoderesults )<-name_vector  
check<- kmoderesults %>% lapply(.,as_tibble, .name_repair = ~ .)%>% lapply(.,melt) 

map(names(kmoderesults),
    function(x) write.xlsx(check[[x]],'results/kmoderesults.xlsx', 
                           sheet=x, append = T))
fviz_cluster(final, data = df)

res.famd <- FAMD(clean, graph = FALSE)
sum_famd <- facto_summarize(res.famd,"var")
fviz_contrib(res.famd, choice = "var", axes = 1, top = 10)
fviz_cos2(res.famd, choice = "var", axes = 1:2)

write.table(sum_famd, "results/famd_meanEE.tsv", sep ="\t", quote=F, row.names = F)
var <- get_famd_var(res.famd)
# Coordinates of variables
head(var$coord)
# Cos2: quality of representation on the factore map
head(var$cos2)
# Contributions to the  dimensions
head(var$contrib)

# Eigenvalues/variances of dimensions
fviz_screeplot(res.famd)
# Graph of variables
fviz_famd_var(res.famd)
# Quantitative variables
fviz_famd_var(res.famd, "quanti.var", repel = TRUE,col.var ='contrib',gradient.cols='nejm')
# Qualitative variables
fviz_famd_var(res.famd, "quali.var", col.var ='contrib',gradient.cols='nejm')
# Graph of individuals colored by cos2
library(ggsci)

# box plot based 1 dimensional scaling
pal_designer <- c("AmplitaqGold"="#DB073D","KapaHiFi"="#DBA507","KapaRobust"="#43A047",
                  "expected"="#07485B","Luna"="#E67E25", "SSII"="#E74C3B", "SSIV"="#23B99A")

fviz_ellipses(res.famd,c('RTase', 'DNA_polymerase'), geom = "point", repel = TRUE, palette =pal_designer,
  addEllipses = TRUE,ellipse.type = "convex",labelsize = 20,pointsize = 2) +
  theme(text = element_text(size =20),
  axis.title = element_text(size = 30),
  axis.text = element_text(size =30)) + 
  ggtitle(label='') 

ggsave('plots/FAMD_compareMeanPe.png', dpi=600)
