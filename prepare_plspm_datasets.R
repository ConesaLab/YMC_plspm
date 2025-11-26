#########################################################
## Prepare PLSPM objects


load('~/Desktop/Thesis/YMC/R/Data/all_plspm_datasets.RData')

atac.data <- t(scale(atac.predicted[,atac.sig.genes], center = T, scale = F))
atac.two.clusters <- kmeans(atac.data, centers = 2)

par(mfrow=c(1,2))
plot(apply(atac.data[names(which(atac.two.clusters$cluster==1)),], 2, mean), 
     type='l', lwd=2, frame.plot=F, xlab='Timepoints', ylab='Intensity')
plot(apply(atac.data[names(which(atac.two.clusters$cluster==2)),], 2, mean), 
     type='l', lwd=2, frame.plot=F, xlab='Timepoints', ylab='Intensity')

atac.two.clusters <- atac.two.clusters$cluster
atac.two.clusters[atac.two.clusters==1] <- "HOC"; atac.two.clusters[atac.two.clusters==2] <- "LOC"

ATAC.loc = t(atac.data[names(which(atac.two.clusters=="LOC")),])
ATAC.hoc = t(atac.data[names(which(atac.two.clusters=="HOC")),])

####
rna.data <- t(scale(rna.predicted[,names(rna.CLUSTERS.scaled$cluster)], center = T, scale = F))
rna.three.clusters <- kmeans(rna.data, centers = 3)
par(mfrow=c(1,3))
plot(apply(rna.data[names(which(rna.three.clusters$cluster==1)),], 2, mean), 
     type='l', lwd=2, frame.plot=F, xlab='Timepoints', ylab='Intensity')
plot(apply(rna.data[names(which(rna.three.clusters$cluster==2)),], 2, mean), 
     type='l', lwd=2, frame.plot=F, xlab='Timepoints', ylab='Intensity')
plot(apply(rna.data[names(which(rna.three.clusters$cluster==3)),], 2, mean), 
     type='l', lwd=2, frame.plot=F, xlab='Timepoints', ylab='Intensity')

rna.three.clusters <- rna.three.clusters$cluster
rna.three.clusters[rna.three.clusters==1] <- "OX"; rna.three.clusters[rna.three.clusters==2] <- "RB"
rna.three.clusters[rna.three.clusters==3] <- "RC"

RNA.ox = t(rna.data[names(which(rna.three.clusters=="OX")),])
RNA.rb = t(rna.data[names(which(rna.three.clusters=="RB")),])
RNA.rc = t(rna.data[names(which(rna.three.clusters=="RC")),])

####

k9.data <- t(scale(k9.predicted[,names(new.k9.three.clusters$cluster)], center = T, scale = F))
k9.two.clusters <- kmeans(k9.data, centers = 2)
par(mfrow=c(1,2))
plot(apply(k9.data[names(which(k9.two.clusters$cluster==1)),], 2, mean), 
     type='l', lwd=2, frame.plot=F, xlab='Timepoints', ylab='Intensity')
plot(apply(k9.data[names(which(k9.two.clusters$cluster==2)),], 2, mean), 
     type='l', lwd=2, frame.plot=F, xlab='Timepoints', ylab='Intensity')

k9.two.clusters <- k9.two.clusters$cluster
k9.two.clusters[k9.two.clusters==1] <- "HOC"; k9.two.clusters[k9.two.clusters==2] <- "LOC"

k9.loc = t(k9.data[names(which(k9.two.clusters=="LOC")),])
colnames(k9.loc) <- paste('k9', colnames(k9.loc), sep='_')
k9.hoc = t(k9.data[names(which(k9.two.clusters=="HOC")),])
colnames(k9.loc) <- paste('k9', colnames(k9.loc), sep='_')

####
k18.data <- t(scale(k18.predicted[,names(new.k18.three.clusters$cluster)], center = T, scale = F))
k18.two.clusters <- kmeans(k18.data, centers = 2)
par(mfrow=c(1,2))
plot(apply(k18.data[names(which(k18.two.clusters$cluster==1)),], 2, mean), 
     type='l', lwd=2, frame.plot=F, xlab='Timepoints', ylab='Intensity')
plot(apply(k18.data[names(which(k18.two.clusters$cluster==2)),], 2, mean), 
     type='l', lwd=2, frame.plot=F, xlab='Timepoints', ylab='Intensity')

k18.two.clusters <- k18.two.clusters$cluster
k18.two.clusters[k18.two.clusters==1] <- "HOC"; k18.two.clusters[k18.two.clusters==2] <- "LOC"

k18.loc = t(k18.data[names(which(k18.two.clusters=="LOC")),])
colnames(k18.loc) <- paste('k18', colnames(k18.loc), sep='_')
k18.hoc = t(k18.data[names(which(k18.two.clusters=="HOC")),])
colnames(k18.hoc) <- paste('k18', colnames(k18.hoc), sep='_')
####

net.data <- t(scale(net.predicted[,names(net.three.clusters$cluster)], center = T, scale = F))
net.two.clusters <- kmeans(net.data, centers = 2)
par(mfrow=c(1,2))
plot(apply(net.data[names(which(net.two.clusters$cluster==1)),], 2, mean), 
     type='l', lwd=2, frame.plot=F, xlab='Timepoints', ylab='Intensity')
plot(apply(net.data[names(which(net.two.clusters$cluster==2)),], 2, mean), 
     type='l', lwd=2, frame.plot=F, xlab='Timepoints', ylab='Intensity')

net.two.clusters <- net.two.clusters$cluster
net.two.clusters[net.two.clusters==1] <- "HOC"; net.two.clusters[net.two.clusters==2] <- "LOC"

net.loc = t(net.data[names(which(net.two.clusters=="LOC")),])
net.hoc = t(net.data[names(which(net.two.clusters=="HOC")),])

####
metab.data <- t(scale(metab.predicted[,names(three.metab.cluster$cluster)], center = T, scale = F))
metab.three.clusters <- kmeans(metab.data, centers = 3)
metab.two.clusters <- kmeans(metab.data, centers = 2)
par(mfrow=c(1,3))
plot(apply(metab.data[names(which(metab.three.clusters$cluster==1)),], 2, mean), 
     type='l', lwd=2, frame.plot=F, xlab='Timepoints', ylab='Intensity')
plot(apply(metab.data[names(which(metab.three.clusters$cluster==2)),], 2, mean), 
     type='l', lwd=2, frame.plot=F, xlab='Timepoints', ylab='Intensity')
plot(apply(metab.data[names(which(metab.three.clusters$cluster==3)),], 2, mean), 
     type='l', lwd=2, frame.plot=F, xlab='Timepoints', ylab='Intensity')

metab.three.clusters <- metab.three.clusters$cluster
metab.three.clusters[metab.three.clusters==1] <- "OX"; metab.three.clusters[metab.three.clusters==2] <- "RB"
metab.three.clusters[metab.three.clusters==3] <- "RC"


par(mfrow=c(1,2))
plot(apply(metab.data[names(which(metab.two.clusters$cluster==1)),], 2, mean), 
     type='l', lwd=2, frame.plot=F, xlab='Timepoints', ylab='Intensity')
plot(apply(metab.data[names(which(metab.two.clusters$cluster==2)),], 2, mean), 
     type='l', lwd=2, frame.plot=F, xlab='Timepoints', ylab='Intensity')

metab.two.clusters <- metab.two.clusters$cluster
metab.two.clusters[metab.two.clusters==1] <- "HOC"; metab.two.clusters[metab.two.clusters==2] <- "LOC"

metab.ox = t(metab.data[names(which(metab.three.clusters=="OX")),])
metab.ox1 = metab.ox[,-c(7,11,12)]
metab.ox2 = metab.ox[,c(7,11,12)]
metab.rb = t(metab.data[names(which(metab.three.clusters=="RB")),])
metab.rc = t(metab.data[names(which(metab.three.clusters=="RC")),])
metab.rc1 = metab.rc[,c(8,15,16)]
metab.rc2 = metab.rc[,-c(8,15,16)]
####


save(atac.predicted, k18.predicted, k9.predicted, metab.predicted, net.predicted, rna.predicted, 
     atac.two.clusters, rna.three.clusters, k9.two.clusters, k18.two.clusters, net.two.clusters,
     metab.data, metab.three.clusters, metab.two.clusters, file = '~/Desktop/Thesis/YMC/R/Data/plspm_datasets_clusters.RData')


save(ATAC.loc, ATAC.hoc, RNA.ox, RNA.rb, RNA.rc, k9.loc, k9.hoc, k18.loc, k18.hoc, 
     net.loc, net.hoc, metab.ox1, metab.ox2, metab.rb, metab.rc1, metab.rc2, file = 'Data/plspm_datasets_lvs.RData')

