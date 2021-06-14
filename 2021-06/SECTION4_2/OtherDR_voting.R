####################
#Real Data Analysis#
#setwd('D:\\ATMCS2020\\voting_data\\GCC_politics_voting_19902016')
library(coRanking)
yrs=2006
retainedD=2 #How many dimensional reduced components are retained in the result?
if(yrs==1990){year_leading_vec=c(8,1,10)}
if(yrs==1998){year_leading_vec=c(5,2,4)}
if(yrs==2006){year_leading_vec=c(2,4,3)}
for(i in year_leading_vec[1:retainedD]){
  destfile<-paste0('us.congress.',yrs,'.votes_CircularCoordinates_',lambda,'_',i,'.txt')
  if(!file.exists(destfile)){
    break
  }
  N_cocycle<-i
}
#retainedD=N_cocycle
filein=paste0('us.congress.',yrs,'.votes.txt')
sm1<-read.table(filein,header=F)
sm4=as.matrix(sm1)
storage.mode(sm4) <- "double"
#enforce the sm4 is a double matrix.
#Column denotes different voting reults (1=Yes,-1=No,0=other) while rows denote different senators, in accordance with
#In python/R, we usually use COLUMN denote features and ROW for samples, this is a significant difference.
#################### 
#  True Clustering #
clus<-read.csv(paste0('us.congress.',yrs,'.keys.csv'),header = F)
#clus<-clus[,3]
#plot(sm4[,1],sm4[,2],col=as.factor(clus))

clcl<-c()
for(i in 1:NROW(data_clus)){
  if(clus[i,3]=='D'){clcl[i]<-1;next}
  if(clus[i,3]=='R'){clcl[i]<-2;next}
  if(clus[i,3]=='I'){clcl[i]<-3;next}
}
data_clus<-cbind(sm4[,1],sm4[,2],clcl)
data_clus<-as.data.frame(data_clus)
colnames(data_clus)<-c('v1','v2','grp')
#                  #
####################
######Original Data
library(ggplot2)
######Plot functionality
plot_cluster=function(data, var_cluster, palette, maintx){
  ggplot(data, aes_string(x="V1", y="V2", color=var_cluster)) +
    geom_point(size=1) +
    guides(colour=guide_legend(override.aes=list(size=6))) +
    xlab("axis1") + ylab("axis2") +
    ggtitle(maintx) +
    theme_light(base_size=10) +
    theme(plot.title = element_text(hjust = 0.5),
          aspect.ratio=1,
          legend.direction = "vertical", 
          legend.position = "none",
          legend.box = "vertical") + 
    labs(color = "Cluster")+
    scale_colour_brewer(palette = palette) 
}
######Number of clusters
K=2
set.seed(123)
######t-SNE
#tsne part is adapted from https://github.com/pablo14/post_cluster_tsne
library(caret)
library(Rtsne)
tsne_model = Rtsne(sm4, check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=0.5, dims=retainedD)
#In Rtsne, the data matrix is stored in such a way that each row is a sample, each column is a feature.
#Like TDA, the major computation bottle-neck is the sample size, not dimensionality.
tsne_data = as.data.frame(tsne_model$Y)
#Plotting the results without clustering
tsne_original=tsne_data
tsne_original$cl_TRUTH = factor(data_clus$grp)
#Create k-means clustering model, and assigning the result to the data used to create the tsne
tsne_cluster_kmeans=kmeans(scale(tsne_data), K)
tsne_original$cl_kmeans = factor(tsne_cluster_kmeans$cluster)
#Create hierarchical cluster model, and assigning the result to the data used to create the tsne
tsne_cluster_hierarchical=hclust(dist(scale(tsne_data)))
tsne_original$cl_hierarchical = factor(cutree(tsne_cluster_hierarchical, k=K))
#Plotting
tsne_plot_TRUTH=plot_cluster(tsne_original, "cl_TRUTH","Set1","TRUTH (t-SNE)")
tsne_plot_k=plot_cluster(tsne_original, "cl_kmeans", "Set1","K-means (t-SNE)")
tsne_plot_h=plot_cluster(tsne_original, "cl_hierarchical", "Set1","Hierarchical (t-SNE)")

######UMAP
library(umap)
umap_model = umap(sm4)
umap_data= umap_model$layout
umap_model_conf<-umap_model$config
umap_model_conf$n_components<-retainedD
umap_model = umap(d=sm4,config=umap_model_conf)
umap_original=as.data.frame(umap_data)
umap_original$cl_TRUTH = factor(data_clus$grp)
#Create k-means clustering model, and assigning the result to the data used to create the tsne
umap_cluster_kmeans=kmeans(scale(umap_data), K)
umap_original$cl_kmeans = factor(umap_cluster_kmeans$cluster)
#Create hierarchical cluster model, and assigning the result to the data used to create the tsne
umap_cluster_hierarchical=hclust(dist(scale(umap_data)))
umap_original$cl_hierarchical = factor(cutree(umap_cluster_hierarchical, k=K))
#Plotting
umap_plot_TRUTH=plot_cluster(umap_original, "cl_TRUTH","Set1","TRUTH (UMAP)")
umap_plot_k=plot_cluster(umap_original, "cl_kmeans", "Set1","K-means (UMAP)")
umap_plot_h=plot_cluster(umap_original, "cl_hierarchical", "Set1","Hierarchical (UMAP)")

######PCA
pca_model = prcomp(as.matrix(t(sm4)), retx=T, scale=F)
pca_data=t(pca_model$x)[,1:retainedD]
colnames(pca_data)<-c('V1','V2')
pca_original=as.data.frame(pca_data)
pca_original$cl_TRUTH = factor(data_clus$grp)
#Create k-means clustering model, and assigning the result to the data used to create the tsne
pca_cluster_kmeans=kmeans(scale(pca_data), K)
pca_original$cl_kmeans = factor(pca_cluster_kmeans$cluster)
#Create hierarchical cluster model, and assigning the result to the data used to create the tsne
pca_cluster_hierarchical=hclust(dist(scale(pca_data)))
pca_original$cl_hierarchical = factor(cutree(pca_cluster_hierarchical, k=K))
#Plotting
pca_plot_TRUTH=plot_cluster(pca_original, "cl_TRUTH","Set1","TRUTH (PCA)")
pca_plot_k=plot_cluster(pca_original, "cl_kmeans", "Set1","K-means (PCA)")
pca_plot_h=plot_cluster(pca_original, "cl_hierarchical", "Set1","Hierarchical (PCA)")

######Laplacian Map
#devtools::install_github("gdkrmr/dimRed")
library(dimRed)
#dat <- loadDataSet("3D S Curve")
laplace_data <- as.dimRedData(.~.,sm4)
leim <- LaplacianEigenmaps()
leim@stdpars$ndim<-retainedD
emb <- leim@fun(laplace_data, leim@stdpars)
laplace_data <- emb@data@data
laplace_original=as.data.frame(laplace_data)
colnames(laplace_original)<-c('V1','V2')
laplace_original$cl_TRUTH = factor(data_clus$grp)
#Create k-means clustering model, and assigning the result to the data used to create the tsne
laplace_cluster_kmeans=kmeans(scale(laplace_data), K)
laplace_original$cl_kmeans = factor(laplace_cluster_kmeans$cluster)
#Create hierarchical cluster model, and assigning the result to the data used to create the tsne
laplace_cluster_hierarchical=hclust(dist(scale(laplace_data)))
laplace_original$cl_hierarchical = factor(cutree(laplace_cluster_hierarchical, k=K))
#Plotting
laplace_plot_TRUTH=plot_cluster(laplace_original, "cl_TRUTH","Set1","TRUTH (Laplacian)")
laplace_plot_k=plot_cluster(laplace_original, "cl_kmeans", "Set1","K-means (Laplacian)")
laplace_plot_h=plot_cluster(laplace_original, "cl_hierarchical", "Set1","Hierarchical (Laplacian)")


######
pdf(paste0('OtherDR',yrs,'.pdf'),height=20,width=10)

tsne_plot_TRUTH
tsne_plot_k
tsne_plot_h
umap_plot_TRUTH
umap_plot_k
umap_plot_h
pca_plot_TRUTH
pca_plot_k
pca_plot_h
laplace_plot_TRUTH
laplace_plot_k
laplace_plot_h
gridExtra::grid.arrange(pca_plot_TRUTH, pca_plot_k, pca_plot_h, 
                        tsne_plot_TRUTH, tsne_plot_k, tsne_plot_h, 
                        umap_plot_TRUTH, umap_plot_k, umap_plot_h, 
                        laplace_plot_TRUTH, laplace_plot_k,laplace_plot_h,ncol=3)
dev.off()
pdf(paste0('OtherDR',yrs,'_Combined.pdf'),height=10,width=20)
gridExtra::grid.arrange(pca_plot_TRUTH, pca_plot_k, pca_plot_h, 
                        tsne_plot_TRUTH, tsne_plot_k, tsne_plot_h, 
                        umap_plot_TRUTH, umap_plot_k, umap_plot_h, 
                        laplace_plot_TRUTH, laplace_plot_k,laplace_plot_h,ncol=3)
gridExtra::grid.arrange(pca_plot_TRUTH, tsne_plot_TRUTH, umap_plot_TRUTH, laplace_plot_TRUTH,
                        pca_plot_k, tsne_plot_k, umap_plot_k, laplace_plot_k, 
                        pca_plot_h, tsne_plot_h, umap_plot_h, laplace_plot_h, nrow=3)
dev.off()

#                  #
####################
####################
#Quality Evaluation#
library(coRanking)

Q.pca<-coranking(sm4, pca_data)
Q.tsne<-coranking(sm4, tsne_data)
Q.umap<-coranking(sm4, umap_data)
Q.laplace<-coranking(sm4, laplace_data)

pdf(paste0('OtherDR',yrs,'_QC.pdf'),height=5,width=20)
par(mfrow=c(1,4),pty='s')
imageplot(Q.pca, main = "Coranking matrix of PCA")
imageplot(Q.tsne, main = "Coranking matrix of t-SNE")
imageplot(Q.umap, main = "Coranking matrix of UMAP")
imageplot(Q.laplace, main = "Coranking matrix of Laplacian")
dev.off()

pdf(paste0('GCC',yrs,'_QC.pdf'),height=5,width=15)
par(mfrow=c(1,3),pty='s')

for(lambda in c(0,0.5,1)){
  coords<-read.table(paste0('us.congress.',yrs,'.votes_CircularCoordinates_',lambda,'_0.txt'))
  for(i in year_leading_vec[1:retainedD]){
    destfile<-paste0('us.congress.',yrs,'.votes_CircularCoordinates_',lambda,'_',i,'.txt')
    if(!file.exists(destfile)){
      break
    }
    coords<-cbind(coords,read.table(paste0('us.congress.',yrs,'.votes_CircularCoordinates_',lambda,'_',i,'.txt')))
  }
  Q.GCC<-coranking(sm4, coords)
  imageplot(Q.GCC, main = paste0('Coranking matrix of GCC \n (mod 23 , ',(1-lambda),'*L1 + ',lambda,'*L2)') )
}
dev.off()
#                  #
####################
