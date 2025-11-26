######################################################################################
## Create prediction models to improve the number of observations on each omic
## Apply PLS -path modelling to the multi-omics dataset
library(splines)
library(factoextra)
library(NbClust)
library('clValid')

## Prediction models to match -omics timepoints

x.rna <-  c(1,2,3,4,5,6,8,10,12,14,16,17,21,23,26,29)
x.net <- c(1,2,3,4,7,11,15,18,21,27,29)
load('Data/New_tmm_rnaMatrix.RData')
load('Data/New_tmm_clusters.RData')

rna.matrix.scaled <- cbind(rna.matrix.scaled[,c(2:16)], rna.matrix.scaled[,c(1)])

rna.predicted <- apply(rna.matrix.scaled, 1, function(y){predict.glm(glm(y ~ ns(x.rna,15)), data.frame(x.rna = c(1:30)), 
                                                                     interval = "prediction")})
rna.predicted <- scale(rna.predicted, center = T, scale = T)

plot(predict.glm(glm( rna.matrix.scaled[1,]~ ns(x.rna,15)), data.frame(x.rna = c(1:30)), interval = "prediction"), type = 'l')
points(predict.glm(glm( rna.matrix.scaled[1,]~ ns(x.rna,15)), data.frame(x.rna = c(1:30)), interval = "prediction"), col = 'red', pch = 19, lwd = 1)
points(x.rna, rna.matrix.scaled[1,], pch = 19, lwd = 2)

midist <- get_dist(t(rna.predicted[,names(rna.CLUSTERS.scaled$cluster)]), stand = FALSE, method = "euclidean")

pdf(file = 'dist_rna_heatmap.pdf', width = 6, height = 6)
fviz_dist(midist, show_labels = F,
          gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
dev.off()

fviz_nbclust(t(rna.predicted[,names(rna.CLUSTERS.scaled$cluster)]), kmeans)
fviz_nbclust(t(rna.predicted[,names(rna.CLUSTERS.scaled$cluster)]), kmeans, 'wss')


indices <- c( "kl", "ch", "hartigan",  "cindex", "db", "silhouette", "duda", "pseudot2", "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain", "dunn", "hubert", "sdindex", "dindex", "sdbw")

nc.summary <- c()
for (index in indices){
  print(index)
  res.nbclust <- NbClust(data = t(rna.predicted[,names(rna.CLUSTERS.scaled$cluster)]), diss = midist, distance = NULL, 
                         min.nc = 2, max.nc = 10, 
                         method = "kmeans", index =index)
  nc.summary <- c(nc.summary,res.nbclust$Best.nc[1])
}

table(nc.summary)
dev.off()

pdf('clustering/rna.pdf', width = 5, height = 5)
plot(table(nc.summary), frame.plot=F, ylab='Number of clustering algorithms', xlab='Number of clusters', main='RNA')
dev.off()
plot(silhouette(rna.CLUSTERS.scaled$cluster, midist), col=c('blue3','red3','green3'), border=NA, main = "K-MEDIAS")

###################################################################
## Check the outlier elements and classify in other groups

missclasified <- names(rna.CLUSTERS.scaled$cluster)[silhouette(rna.CLUSTERS.scaled$cluster, midist)[,3]<0]
rna.CLUSTERS.scaled$cluster[missclasified]

plot(t(rna.predicted)[missclasified[3],])

###################################################################


x.hms <- c(1,2,3,4,5,7,9,11,13,16,17,20,22,24,27,29)
load('Data/new_full_hms_matrixAndCluster.RData')
H3k9ac <- cbind(H3k9ac[,c(2:16)], H3k9ac[,c(1)])
H3k9ac <- t(scale(t(H3k9ac), scale = T, center = T))
k9.predicted <- apply(H3k9ac, 1, function(y){predict.glm(glm(y ~ ns(x.hms,15)), data.frame(x.hms = c(1:30)), 
                                                         interval = "prediction")})
k9.predicted <- scale(k9.predicted, center = T, scale = T)

plot(predict.glm(glm( H3k9ac[1,]~ ns(x.hms,15)), data.frame(x.hms = c(1:30)), interval = "prediction"), type = 'l')
points(predict.glm(glm( H3k9ac[1,]~ ns(x.hms,15)), data.frame(x.hms = c(1:30)), interval = "prediction"), col = 'red', pch = 19, lwd = 1)
points(x.hms, H3k9ac[1,], pch = 19, lwd = 2)

midist <- get_dist(t(k9.predicted[,names(new.k9.three.clusters$cluster)]), stand = FALSE, method = "euclidean")

pdf(file = 'dist_k9_heatmap.pdf', width = 6, height = 6)
fviz_dist(midist, show_labels = F,
          gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
dev.off()

fviz_nbclust(t(k9.predicted[,names(new.k9.three.clusters$cluster)]), kmeans)
fviz_nbclust(t(k9.predicted[,names(new.k9.three.clusters$cluster)]), kmeans, 'wss')


indices <- c( "kl", "ch", "hartigan",  "cindex", "db", "silhouette", "duda", "pseudot2", "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain", "dunn", "hubert", "sdindex", "dindex", "sdbw")

nc.summary <- c()
for (index in indices){
  print(index)
  res.nbclust <- NbClust(data = t(k9.predicted[,names(new.k9.three.clusters$cluster)]), diss = midist, distance = NULL, 
                         min.nc = 2, max.nc = 10, 
                         method = "kmeans", index =index)
  nc.summary <- c(nc.summary,res.nbclust$Best.nc[1])
}

table(nc.summary)
pdf('clustering/k9.pdf', width = 5, height = 5)
plot(table(nc.summary), frame.plot=F, ylab='Number of clustering algorithms', xlab='Number of clusters', main='H3K9ac')
dev.off()


H3k18ac <- cbind(H3k18ac[,c(2:16)], H3k18ac[,c(1)])
H3k18ac <- H3k18ac[as.logical(rowSums(H3k18ac != 0)), ]
H3k18ac <- t(scale(t(H3k18ac), scale = T, center = T))
k18.predicted <- apply(H3k18ac, 1, function(y){predict.glm(glm(y ~ ns(x.hms,15)), data.frame(x.hms = c(1:30)), 
                                                           interval = "prediction")})
k18.predicted <- scale(k18.predicted, center = T, scale = T)

plot(predict.glm(glm( H3k18ac[1,]~ ns(x.hms,15)), data.frame(x.hms = c(1:30)), interval = "prediction"), type = 'l')
points(predict.glm(glm( H3k18ac[1,]~ ns(x.hms,15)), data.frame(x.hms = c(1:30)), interval = "prediction"), col = 'red', pch = 19, lwd = 1)
points(x.hms, H3k18ac[1,], pch = 19)

midist <- get_dist(t(k18.predicted[,names(new.k18.three.clusters$cluster)]), stand = FALSE, method = "euclidean")

pdf(file = 'dist_k18_heatmap.pdf', width = 6, height = 6)
fviz_dist(midist, show_labels = F,
          gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
dev.off()

fviz_nbclust(t(k18.predicted[,names(new.k18.three.clusters$cluster)]), kmeans)
fviz_nbclust(t(k18.predicted[,names(new.k18.three.clusters$cluster)]), kmeans, 'wss')


indices <- c( "kl", "ch", "hartigan",  "cindex", "db", "silhouette", "duda", "pseudot2", "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain", "dunn", "hubert", "sdindex", "dindex", "sdbw")

nc.summary <- c()
for (index in indices){
  print(index)
  res.nbclust <- NbClust(data = t(k18.predicted[,names(new.k18.three.clusters$cluster)]), diss = midist, distance = NULL, 
                         min.nc = 2, max.nc = 10, 
                         method = "kmeans", index =index)
  nc.summary <- c(nc.summary,res.nbclust$Best.nc[1])
}

table(nc.summary)
pdf('clustering/k18.pdf', width = 5, height = 5)
plot(table(nc.summary), frame.plot=F, ylab='Number of clustering algorithms', xlab='Number of clusters', main='H3K18ac')
dev.off()




#########################################################################
##

x.metab <- c(1,2,3,4,5,6,8,9,11,13,15,16,17,18,20,22,23,25,28,29,30)
load('Data/metab_unique_matrix.RData')
load('Data/metabs_clusters.RData')

metab.matrix <- t(scale(median.metab.matrix, center = T, scale = T))
metab.matrix <- cbind(metab.matrix[,c(7:21)], metab.matrix[,c(1:6)])
metab.predicted <- apply(metab.matrix, 1, function(y){predict.glm(glm(y ~ ns(x.metab,20)), data.frame(x.metab = c(1:30)), 
                                                                  interval = "prediction")})
metab.predicted <- scale(metab.predicted, center = T, scale = T)

plot(predict.glm(glm( metab.matrix[1,]~ ns(x.metab,20)), data.frame(x.metab = c(1:30)), interval = "prediction"), type = 'l')
points(predict.glm(glm( metab.matrix[1,]~ ns(x.metab,20)), data.frame(x.metab = c(1:30)), interval = "prediction"), col = 'red', pch = 19, lwd = 1)
points(x.metab, metab.matrix[1,], pch = 19)

midist <- get_dist(t(metab.predicted[,names(three.metab.cluster$cluster)]), stand = FALSE, method = "euclidean")

pdf(file = 'dist_metab_heatmap.pdf', width = 6, height = 6)
fviz_dist(midist, show_labels = F,
          gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
dev.off()

fviz_nbclust(t(metab.predicted[,names(three.metab.cluster$cluster)]), kmeans)
fviz_nbclust(t(metab.predicted[,names(three.metab.cluster$cluster)]), kmeans, 'wss')


indices <- c( "kl", "ch", "hartigan",  "cindex", "db", "silhouette", "duda", "pseudot2", "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain", "dunn", "hubert", "sdindex", "dindex", "sdbw")

nc.summary <- c()
for (index in indices){
  print(index)
  res.nbclust <- NbClust(data = t(metab.predicted[,names(three.metab.cluster$cluster)]), diss = midist, distance = NULL, 
                         min.nc = 2, max.nc = 10, 
                         method = "kmeans", index =index)
  nc.summary <- c(nc.summary,res.nbclust$Best.nc[1])
}

table(nc.summary)
pdf('clustering/metabs.pdf', width = 5, height = 5)
plot(table(nc.summary), frame.plot=F, ylab='Number of clustering algorithms', xlab='Number of clusters', main='Metabolites')
dev.off()

metabs.clusters <- kmeans(t(metab.predicted[,names(three.metab.cluster$cluster)]), centers = 3)

plot(silhouette(metabs.clusters$cluster, midist), col=c('blue3','red3', 'green3'), border=NA, main = "K-MEDIAS")

par(mfrow=c(1,3))
plot(t(metab.predicted)['Lactate',], type='l', frame.plot=F)
plot(t(metab.predicted)['succinate',], type='l', frame.plot=F)
plot(t(metab.predicted)['Citrulline',], type='l', frame.plot=F)

plot(t(metab.predicted)['citrate',], type='l', frame.plot=F)
plot(t(metab.predicted)['TMP',], type='l', frame.plot=F)
plot(t(metab.predicted)['NAMN',], type='l', frame.plot=F)

plot(t(metab.predicted)['fumarate',], type='l', frame.plot=F)
plot(t(metab.predicted)['3-HBA',], type='l', frame.plot=F)
plot(t(metab.predicted)['acetyl CoA',], type='l', frame.plot=F)

###################################################################
## Check the outlier elements and classify in other groups

missclasified <- names(three.metab.cluster$cluster)[silhouette(three.metab.cluster$cluster, midist)[,3]<0]
three.metab.cluster$cluster[missclasified]

plot(t(metab.predicted)[missclasified[3],])

###################################################################



x.net <- c(1,2,3,4,7,11,15,18,21,27,29)
load('Data/Net_tmm_center.RData')
load('Data/net_clusters.RData')

net.data <- t(scale(t(net.center), center = T,scale = T))
net.data <- net.data[complete.cases(net.data),]
net.predicted <- apply(net.data, 1, function(y){predict.glm(glm(y ~ ns(x.net,10)), data.frame(x.net = c(1:30)), 
                                                            interval = "prediction")})
net.predicted <- scale(net.predicted, center = T, scale = T)

plot(predict.glm(glm( net.data[1,]~ ns(x.net,10)), data.frame(x.net = c(1:30)), interval = "prediction"), type = 'l')
points(predict.glm(glm( net.data[1,]~ ns(x.net,10)), data.frame(x.net = c(1:30)), interval = "prediction"), col = 'red', pch = 19, lwd = 1)
points(x.net, net.data[1,], pch = 19)

midist <- get_dist(t(net.predicted[,names(net.three.clusters$cluster)]), stand = FALSE, method = "euclidean")

pdf(file = 'dist_net_heatmap.pdf', width = 6, height = 6)
fviz_dist(midist, show_labels = F,
          gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
dev.off()

#midist <- get_dist(net.data, stand = FALSE, method = "euclidean")

pdf(file = 'dist_net_original_heatmap.pdf', width = 6, height = 6)
fviz_dist(midist, show_labels = F,
          gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
dev.off()

fviz_nbclust(t(net.predicted[,names(net.three.clusters$cluster)]), kmeans)
fviz_nbclust(t(net.predicted[,names(net.three.clusters$cluster)]), kmeans, 'wss')


indices <- c( "kl", "ch", "hartigan",  "cindex", "db", "silhouette", "duda", "pseudot2", "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain", "dunn", "hubert", "sdindex", "dindex", "sdbw")

nc.summary <- c()
for (index in indices){
  print(index)
  res.nbclust <- NbClust(data = t(net.predicted[,names(net.three.clusters$cluster)]), diss = midist, distance = NULL, 
                         min.nc = 2, max.nc = 10, 
                         method = "kmeans", index =index)
  nc.summary <- c(nc.summary,res.nbclust$Best.nc[1])
}

table(nc.summary)
pdf('clustering/net.pdf', width = 5, height = 5)
plot(table(nc.summary), frame.plot=F, ylab='Number of clustering algorithms', xlab='Number of clusters', main='NET')
dev.off()

##################################################################
load('Data/atac_masigpro.RData')
x.atac <- c(1,2,3,4,5,9,14,19,22,24,28,29)
load('Data/atac_batch_corrected.RData')
atac.data <- t(scale(t(batch.tmm.atac.matrix[,-c(11,14)][,c(c(6:12),c(1:5))]), center = T,scale = T))
atac.predicted <- apply(atac.data, 1, function(y){predict.glm(glm(y ~ ns(x.atac,11)), data.frame(x.atac = c(1:30)), 
                                                              interval = "prediction")})
atac.predicted <- scale(atac.predicted, center = T, scale = T)

plot(predict.glm(glm( atac.data[13,]~ ns(x.atac,11)), data.frame(x.atac = c(1:30)), interval = "prediction"), type = 'l')
points(predict.glm(glm( atac.data[13,]~ ns(x.atac,11)), data.frame(x.atac = c(1:30)), interval = "prediction"), col = 'red', pch = 19, lwd = 1)
points(x.atac, atac.data[13,], pch = 19)

midist <- get_dist(t(atac.predicted[,sigs$summary]), stand = FALSE, method = "euclidean")

pdf(file = 'dist_atac_heatmap.pdf', width = 6, height = 6)
fviz_dist(midist, show_labels = F,
          gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
dev.off()


indices <- c( "kl", "ch", "hartigan",  "cindex", "db", "silhouette", "duda", "pseudot2", 
              "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain", 
              "dunn", "hubert", "sdindex", "dindex", "sdbw")

nc.summary <- c()
for (index in indices){
  print(index)
  res.nbclust <- NbClust(data = t(atac.predicted[,sigs$summary]), diss = midist, distance = NULL, 
                         min.nc = 2, max.nc = 10, 
                         method = "kmeans", index =index)
  nc.summary <- c(nc.summary,res.nbclust$Best.nc[1])
}

table(nc.summary)
pdf('clustering/atac.pdf', width = 5, height = 5)
plot(table(nc.summary), frame.plot=F, ylab='Number of clustering algorithms', xlab='Number of clusters', main='ATAC')
dev.off()

atac.sig.genes <- sigs$summary

two.clusters <- kmeans(t(scale(t(sigs$sig.genes$sig.profiles), center = T, scale = T)), centers = 2)
par(mfrow=c(1,2))
plot(apply(t(atac.predicted)[names(two.clusters$cluster)[two.clusters$cluster==1],],2,mean), type='l', frame.plot=F, ylab='', main='ATAC Cluster1', lwd=1.5)
plot(apply(t(atac.predicted)[names(two.clusters$cluster)[two.clusters$cluster==2],],2,mean), type='l', frame.plot=F, ylab='', main='ATAC Cluster2', lwd=1.5)



save(rna.predicted, k9.predicted, k18.predicted, metab.predicted, net.predicted, atac.predicted, rna.CLUSTERS.scaled, net.three.clusters, 
     three.metab.cluster, new.k9.three.clusters, new.k18.three.clusters, atac.sig.genes, file = 'Data/all_plspm_datasets.RData')



############################################
## MORE

source('~/Desktop/Thesis/YMC/R/MORE/ComputeGLM_function.R')
source('~/Desktop/Thesis/YMC/R/MORE/MORE_GLM.R')
source('~/Desktop/Thesis/YMC/R/MORE/auxFunctions.R')

load('Data/all_plspm_datasets.RData')

## Associate omics

k9.genes <- intersect(names(new.k9.three.clusters$cluster), names(rna.CLUSTERS.scaled$cluster))
k9.association <- as.data.frame(cbind(k9.genes, paste0('k9_',k9.genes)))

k18.genes <- intersect(names(new.k18.three.clusters$cluster), names(rna.CLUSTERS.scaled$cluster))
k18.association <- as.data.frame(cbind(k18.genes, paste0('k18_',k18.genes)))

net.genes <- intersect(names(net.three.clusters$cluster), names(rna.CLUSTERS.scaled$cluster))
net.association <- as.data.frame(cbind(net.genes, paste0('net_',net.genes)))

atac.genes <- intersect(atac.sig.genes, names(rna.CLUSTERS.scaled$cluster))
atac.association <- as.data.frame(cbind(atac.genes, paste0('atac_',atac.genes)))

associations <- list('k18'=k18.association, 'k9'=k9.association, 'atac'=atac.association, 'net'=net.association)

mod.k18.predicted <- t(scale(k18.predicted, center = T, scale = T))[k18.genes,]; rownames(mod.k18.predicted) <- paste0('k18_', k18.genes)
mod.k9.predicted <- t(scale(k9.predicted, center = T, scale = T))[k9.genes,]; rownames(mod.k9.predicted) <- paste0('k9_', k9.genes)
mod.atac.predicted <- t(scale(atac.predicted, center = T, scale = T))[atac.genes,]; rownames(mod.atac.predicted) <- paste0('atac_', atac.genes)
mod.net.predicted <- t(scale(net.predicted, center = T, scale = T))[net.genes,]; rownames(mod.net.predicted) <- paste0('net_', net.genes)



omics <- list('k18'=mod.k18.predicted,'k9'=mod.k9.predicted,'atac'=mod.atac.predicted,'net'=mod.net.predicted)

GLMresults <- GetGLM(GeneExpression = t(rna.predicted)[names(rna.CLUSTERS.scaled$cluster),],
                     associations = associations,
                     data.omics = omics,
                     edesign = NULL, Res.df = 8,
                     epsilon = 0.00001, alfa = 0.05,
                     MT.adjust = 'fdr', family = gaussian(),
                     elasticnet = 2, center = F, scale = F,
                     stepwise = "two.ways.backward",
                     interactions.exp = T, interactions.reg = 1,
                     correlation = 0.95, action = 'mean',
                     cont.var = NULL, min.obs = 10, min)



k18.targets <- c(); k9.targets <- c()
atac.targets <- c(); net.targets <- c()

for (gene in names(GLMresults$ResultsPerGene)){
  sigRegs <- GLMresults$ResultsPerGene[[gene]]$significantRegulators
  if (length(grep('^k18', sigRegs))>0){
    k18.targets <- c(k18.targets,gene)
  }
  if (length(grep('^k9', sigRegs))>0){
    k9.targets <- c(k9.targets,gene)
  }
  if (length(grep('^atac', sigRegs))>0){
    atac.targets <- c(atac.targets,gene)
  }
  if (length(grep('^net', sigRegs))>0){
    net.targets <- c(net.targets,gene)  
  }
}

ox.nums<- c(sum(k18.targets%in%names(rna.CLUSTERS.scaled$cluster)[rna.CLUSTERS.scaled$cluster==2]),
            sum(k9.targets%in%names(rna.CLUSTERS.scaled$cluster)[rna.CLUSTERS.scaled$cluster==2]),
            sum(atac.targets%in%names(rna.CLUSTERS.scaled$cluster)[rna.CLUSTERS.scaled$cluster==2]),
            sum(net.targets%in%names(rna.CLUSTERS.scaled$cluster)[rna.CLUSTERS.scaled$cluster==2]))

rb.nums<- c(sum(k18.targets%in%names(rna.CLUSTERS.scaled$cluster)[rna.CLUSTERS.scaled$cluster==3]),
            sum(k9.targets%in%names(rna.CLUSTERS.scaled$cluster)[rna.CLUSTERS.scaled$cluster==3]),
            sum(atac.targets%in%names(rna.CLUSTERS.scaled$cluster)[rna.CLUSTERS.scaled$cluster==3]),
            sum(net.targets%in%names(rna.CLUSTERS.scaled$cluster)[rna.CLUSTERS.scaled$cluster==3]))

rc.nums<- c(sum(k18.targets%in%names(rna.CLUSTERS.scaled$cluster)[rna.CLUSTERS.scaled$cluster==1]),
            sum(k9.targets%in%names(rna.CLUSTERS.scaled$cluster)[rna.CLUSTERS.scaled$cluster==1]),
            sum(atac.targets%in%names(rna.CLUSTERS.scaled$cluster)[rna.CLUSTERS.scaled$cluster==1]),
            sum(net.targets%in%names(rna.CLUSTERS.scaled$cluster)[rna.CLUSTERS.scaled$cluster==1]))



cluster.results <- rbind(ox.nums,rb.nums,rc.nums)

rownames(cluster.results) <- c('OX','RB','RC')
colnames(cluster.results) <- c('K18ac','K9ac','ATAC','NET')


barplot(t(cluster.results), beside = T, col=c(rep('tomato3',4),rep('#50B547',4),rep('cornflowerblue',4)))

barplot(t(cluster.results/c(1163,662,1268)), beside = T, col=c(rep('tomato3',4),rep('#50B547',4),rep('cornflowerblue',4)))
barplot(cluster.results/c(1163,662,1268), beside = T, col=c(rep('tomato3',1),rep('#50B547',1),rep('cornflowerblue',1)))



###############################################################################


common <- intersect(atac.sig.genes, names(new.k9.three.clusters$cluster))
k9.non.diff <- k9.predicted[!rownames(k9.predicted)%in%paste('k9', names(new.k9.three.clusters$cluster), sep='_'),]
common.non.diff <- intersect(colnames(k9.non.diff), colnames(atac.predicted))


boxplot(abs(diag(cor(t(mod.k9.predicted[paste('k9', common, sep='_'),]),t(mod.atac.predicted[paste('atac', common, sep='_'),])))),
        abs(cor(k9.predicted[,common.non.diff],atac.predicted[,common.non.diff])))









