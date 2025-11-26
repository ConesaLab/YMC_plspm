
load('../../Thesis/R_2/Data/plspm_datasets_clusters.RData')

k9.matrix.hoc <- t(k9.predicted)[names(k9.two.clusters)[k9.two.clusters=="HOC"],]
k9.matrix.loc <- t(k9.predicted)[names(k9.two.clusters)[k9.two.clusters=="LOC"],]
k18.matrix.hoc <- t(k18.predicted)[names(k18.two.clusters)[k18.two.clusters=="HOC"],]
k18.matrix.loc <- t(k18.predicted)[names(k18.two.clusters)[k18.two.clusters=="LOC"],]


common.genes.hoc <- intersect(rownames(k9.matrix.hoc), rownames(k18.matrix.hoc))
common.genes.loc <- intersect(rownames(k9.matrix.loc), rownames(k18.matrix.loc))
cors.list.hoc<-c(); cors.list.loc<-c()

for (gene in common.genes.hoc){
  cors.list.hoc<- c(cors.list.hoc, cor(k9.matrix.hoc[gene,], k18.matrix.hoc[gene,]))
}

for (gene in common.genes.loc){
  cors.list.loc<- c(cors.list.loc, cor(k9.matrix.loc[gene,], k18.matrix.loc[gene,]))
}

pdf('histone_correlation_hoc_density.pdf', width = 4, height = 4)
plot(density(cors.list.hoc), axes=F,main='HOC Correlation between H3K18ac and H3K9ac', frame.plot=F, xlim=c(0.2,1), ylim=c(0,8))
axis(side = 1, lwd = 2) # Thicker bottom axis
axis(side = 2, lwd = 2) # Thicker left axis
abline(v = median(cors.list.hoc), col='red', lwd=2, lty = 2)
dev.off()

pdf('histone_correlation_loc_density.pdf', width = 4, height = 4)
plot(density(cors.list.loc), axes=F, main='LOC Correlation between H3K18ac and H3K9ac', xlim=c(0.2,1), ylim=c(0,8))
abline(v = median(cors.list.loc), col='red', lwd=2, lty = 2)
axis(side = 1, lwd = 2) # Thicker bottom axis
axis(side = 2, lwd = 2) # Thicker left axis
dev.off()




