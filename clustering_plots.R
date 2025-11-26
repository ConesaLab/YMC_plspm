######################################################
## Plot the clusters - different clusters for each omic


load('~/Desktop/Thesis/YMC/R/Data/plspm_datasets_clusters.RData')


library(scales)
## RNA 


pdf('../Figures/Profiles/rnaseq_clusters.pdf', width = 5.5, height = 5)

#OX
rna.ox.median <- apply(t(rna.predicted)[names(which(rna.three.clusters=='OX')),], 2, median)
rna.ox.sd<- apply(t(rna.predicted)[names(which(rna.three.clusters=='OX')),], 2, sd)

plot(rna.ox.median, type = 'l', lwd=2, col='tomato3', frame.plot=F, ylim=c(-2,3), ylab='Expression value', xlab='Timepoint', main='RNA-seq')
points(rna.ox.median, lwd=1, pch=19, col='tomato3')
polygon(c(1:30, rev(1:30)), c(rna.ox.median-rna.ox.sd,rev(rna.ox.median+rna.ox.sd)), border = F, col=alpha('tomato3', 0.3))
points(rna.ox.median+rna.ox.sd, type = 'l', lwd=0.75, col='tomato3')
points(rna.ox.median-rna.ox.sd, type = 'l', lwd=0.75, col='tomato3')

#RB
rna.rb.median <- apply(t(rna.predicted)[names(which(rna.three.clusters=='RB')),], 2, median)
rna.rb.sd<- apply(t(rna.predicted)[names(which(rna.three.clusters=='RB')),], 2, sd)

lines(rna.rb.median, type = 'l', lwd=2, col='#50B547', ylim=c(-2,3))
points(rna.rb.median, lwd=1, pch=19, col='#50B547')
polygon(c(1:30, rev(1:30)), c(rna.rb.median-rna.rb.sd,rev(rna.rb.median+rna.rb.sd)), border = F, col=alpha('#50B547', 0.3))
points(rna.rb.median+rna.rb.sd, type = 'l', lwd=0.75, col='#50B547')
points(rna.rb.median-rna.rb.sd, type = 'l', lwd=0.75, col='#50B547')

#RC
rna.rc.median <- apply(t(rna.predicted)[names(which(rna.three.clusters=='RC')),], 2, median)
rna.rc.sd<- apply(t(rna.predicted)[names(which(rna.three.clusters=='RC')),], 2, sd)

lines(rna.rc.median, type = 'l', lwd=2, col='cornflowerblue', ylim=c(-2,3))
points(rna.rc.median, lwd=1, pch=19, col='cornflowerblue')
polygon(c(1:30, rev(1:30)), c(rna.rc.median-rna.rc.sd,rev(rna.rc.median+rna.rc.sd)), border = F, col=alpha('cornflowerblue', 0.3))
points(rna.rc.median+rna.rc.sd, type = 'l', lwd=0.75, col='cornflowerblue')
points(rna.rc.median-rna.rc.sd, type = 'l', lwd=0.75, col='cornflowerblue')

dev.off()


##################################
## H3K9ac


pdf('../Figures/Profiles/k9ac_clusters.pdf', width = 5.5, height = 5)

#HOC
k9ac.hoc.median <- apply(t(k9.predicted)[names(which(k9.two.clusters=='HOC')),], 2, median)
k9ac.hoc.sd<- apply(t(k9.predicted)[names(which(k9.two.clusters=='HOC')),], 2, sd)

plot(k9ac.hoc.median, type = 'l', lwd=2, col='darkorange', frame.plot=F, ylim=c(-2,3), ylab='Expression value', xlab='Timepoint', main='H3K9ac')
points(k9ac.hoc.median, lwd=1, pch=19, col='darkorange')
polygon(c(1:30, rev(1:30)), c(k9ac.hoc.median-k9ac.hoc.sd,rev(k9ac.hoc.median+k9ac.hoc.sd)), border = F, col=alpha('darkorange', 0.3))
points(k9ac.hoc.median+k9ac.hoc.sd, type = 'l', lwd=0.75, col='darkorange')
points(k9ac.hoc.median-k9ac.hoc.sd, type = 'l', lwd=0.75, col='darkorange')

#LOC
k9ac.loc.median <- apply(t(k9.predicted)[names(which(k9.two.clusters=='LOC')),], 2, median)
k9ac.loc.sd<- apply(t(k9.predicted)[names(which(k9.two.clusters=='LOC')),], 2, sd)

lines(k9ac.loc.median, type = 'l', lwd=2, col='yellow3', ylim=c(-2,3))
points(k9ac.loc.median, lwd=1, pch=19, col='yellow4')
polygon(c(1:30, rev(1:30)), c(k9ac.loc.median-k9ac.loc.sd,rev(k9ac.loc.median+k9ac.loc.sd)), border = F, col=alpha('yellow4', 0.3))
points(k9ac.loc.median+k9ac.loc.sd, type = 'l', lwd=0.75, col='yellow4')
points(k9ac.loc.median-k9ac.loc.sd, type = 'l', lwd=0.75, col='yellow4')

dev.off()



##################################
## H3K18ac


pdf('../Figures/Profiles/k18ac_clusters.pdf', width = 5.5, height = 5)

#HOC
k18ac.hoc.median <- apply(t(k18.predicted)[names(which(k18.two.clusters=='HOC')),], 2, median)
k18ac.hoc.sd<- apply(t(k18.predicted)[names(which(k18.two.clusters=='HOC')),], 2, sd)

plot(k18ac.hoc.median, type = 'l', lwd=2, col='darkorange', frame.plot=F, ylim=c(-2,3), ylab='Expression value', xlab='Timepoint', main='H3K18ac')
points(k18ac.hoc.median, lwd=1, pch=19, col='darkorange')
polygon(c(1:30, rev(1:30)), c(k18ac.hoc.median-k18ac.hoc.sd,rev(k18ac.hoc.median+k18ac.hoc.sd)), border = F, col=alpha('darkorange', 0.3))
points(k18ac.hoc.median+k18ac.hoc.sd, type = 'l', lwd=0.75, col='darkorange')
points(k18ac.hoc.median-k18ac.hoc.sd, type = 'l', lwd=0.75, col='darkorange')

#LOC
k18ac.loc.median <- apply(t(k18.predicted)[names(which(k18.two.clusters=='LOC')),], 2, median)
k18ac.loc.sd<- apply(t(k18.predicted)[names(which(k18.two.clusters=='LOC')),], 2, sd)

lines(k18ac.loc.median, type = 'l', lwd=2, col='yellow3', ylim=c(-2,3))
points(k18ac.loc.median, lwd=1, pch=19, col='yellow4')
polygon(c(1:30, rev(1:30)), c(k18ac.loc.median-k18ac.loc.sd,rev(k18ac.loc.median+k18ac.loc.sd)), border = F, col=alpha('yellow4', 0.3))
points(k18ac.loc.median+k18ac.loc.sd, type = 'l', lwd=0.75, col='yellow4')
points(k18ac.loc.median-k18ac.loc.sd, type = 'l', lwd=0.75, col='yellow4')

dev.off()


##################################
## NET-seq


pdf('../Figures/Profiles/net_clusters.pdf', width = 5.5, height = 5)

#HOC
net.hoc.median <- apply(t(net.predicted)[names(which(net.two.clusters=='HOC')),], 2, median)
net.hoc.sd<- apply(t(net.predicted)[names(which(net.two.clusters=='HOC')),], 2, sd)

plot(net.hoc.median, type = 'l', lwd=2, col='darkorange', frame.plot=F, ylim=c(-2,3), ylab='Expression value', xlab='Timepoint', main='NET-seq')
points(net.hoc.median, lwd=1, pch=19, col='darkorange')
polygon(c(1:30, rev(1:30)), c(net.hoc.median-net.hoc.sd,rev(net.hoc.median+net.hoc.sd)), border = F, col=alpha('darkorange', 0.3))
points(net.hoc.median+net.hoc.sd, type = 'l', lwd=0.75, col='darkorange')
points(net.hoc.median-net.hoc.sd, type = 'l', lwd=0.75, col='darkorange')

#LOC
net.loc.median <- apply(t(net.predicted)[names(which(net.two.clusters=='LOC')),], 2, median)
net.loc.sd<- apply(t(net.predicted)[names(which(net.two.clusters=='LOC')),], 2, sd)

lines(net.loc.median, type = 'l', lwd=2, col='yellow3', ylim=c(-2,3))
points(net.loc.median, lwd=1, pch=19, col='yellow4')
polygon(c(1:30, rev(1:30)), c(net.loc.median-net.loc.sd,rev(net.loc.median+net.loc.sd)), border = F, col=alpha('yellow4', 0.3))
points(net.loc.median+net.loc.sd, type = 'l', lwd=0.75, col='yellow4')
points(net.loc.median-net.loc.sd, type = 'l', lwd=0.75, col='yellow4')

dev.off()



##################################
## ATAC-seq


pdf('../Figures/Profiles/atac_clusters.pdf', width = 5.5, height = 5)

#HOC
atac.hoc.median <- apply(t(atac.predicted)[names(which(atac.two.clusters=='HOC')),], 2, median)
atac.hoc.sd<- apply(t(atac.predicted)[names(which(atac.two.clusters=='HOC')),], 2, sd)

plot(atac.hoc.median, type = 'l', lwd=2, col='darkorange', frame.plot=F, ylim=c(-2,3), ylab='Expression value', xlab='Timepoint', main='ATAC')
points(atac.hoc.median, lwd=1, pch=19, col='darkorange')
polygon(c(1:30, rev(1:30)), c(atac.hoc.median-atac.hoc.sd,rev(atac.hoc.median+atac.hoc.sd)), border = F, col=alpha('darkorange', 0.3))
points(atac.hoc.median+atac.hoc.sd, type = 'l', lwd=0.75, col='darkorange')
points(atac.hoc.median-atac.hoc.sd, type = 'l', lwd=0.75, col='darkorange')

#LOC
atac.loc.median <- apply(t(atac.predicted)[names(which(atac.two.clusters=='LOC')),], 2, median)
atac.loc.sd<- apply(t(atac.predicted)[names(which(atac.two.clusters=='LOC')),], 2, sd)

lines(atac.loc.median, type = 'l', lwd=2, col='yellow3', ylim=c(-2,3))
points(atac.loc.median, lwd=1, pch=19, col='yellow4')
polygon(c(1:30, rev(1:30)), c(atac.loc.median-atac.loc.sd,rev(atac.loc.median+atac.loc.sd)), border = F, col=alpha('yellow4', 0.3))
points(atac.loc.median+atac.loc.sd, type = 'l', lwd=0.75, col='yellow4')
points(atac.loc.median-atac.loc.sd, type = 'l', lwd=0.75, col='yellow4')

dev.off()






##################################
## Metabolites


## Two clusters

pdf('../Figures/Profiles/metab_two_clusters.pdf', width = 5.5, height = 5)

#HOC
metab.hoc.median <- apply(t(metab.predicted)[names(which(metab.two.clusters=='HOC')),], 2, median)
metab.hoc.sd<- apply(t(metab.predicted)[names(which(metab.two.clusters=='HOC')),], 2, sd)

plot(metab.hoc.median, type = 'l', lwd=2, col='darkorange', frame.plot=F, ylim=c(-2,3), ylab='Expression value', xlab='Timepoint', main='Metabolites')
points(metab.hoc.median, lwd=1, pch=19, col='darkorange')
polygon(c(1:30, rev(1:30)), c(metab.hoc.median-metab.hoc.sd,rev(metab.hoc.median+metab.hoc.sd)), border = F, col=alpha('darkorange', 0.3))
points(metab.hoc.median+metab.hoc.sd, type = 'l', lwd=0.75, col='darkorange')
points(metab.hoc.median-metab.hoc.sd, type = 'l', lwd=0.75, col='darkorange')

#LOC
metab.loc.median <- apply(t(metab.predicted)[names(which(metab.two.clusters=='LOC')),], 2, median)
metab.loc.sd<- apply(t(metab.predicted)[names(which(metab.two.clusters=='LOC')),], 2, sd)

lines(metab.loc.median, type = 'l', lwd=2, col='yellow3', ylim=c(-2,3))
points(metab.loc.median, lwd=1, pch=19, col='yellow4')
polygon(c(1:30, rev(1:30)), c(metab.loc.median-metab.loc.sd,rev(metab.loc.median+metab.loc.sd)), border = F, col=alpha('yellow4', 0.3))
points(metab.loc.median+metab.loc.sd, type = 'l', lwd=0.75, col='yellow4')
points(metab.loc.median-metab.loc.sd, type = 'l', lwd=0.75, col='yellow4')

dev.off()


##Three clusters

pdf('../Figures/Profiles/metab_three_clusters.pdf', width = 5.5, height = 5)

#OX
metab.ox.median <- apply(t(metab.predicted)[names(which(metab.three.clusters=='OX')),], 2, median)
metab.ox.sd<- apply(t(metab.predicted)[names(which(metab.three.clusters=='OX')),], 2, sd)

plot(metab.ox.median, type = 'l', lwd=2, col='tomato3', frame.plot=F, ylim=c(-2,3), ylab='Expression value', xlab='Timepoint', main='Metabolites')
points(metab.ox.median, lwd=1, pch=19, col='tomato3')
polygon(c(1:30, rev(1:30)), c(metab.ox.median-metab.ox.sd,rev(metab.ox.median+metab.ox.sd)), border = F, col=alpha('tomato3', 0.3))
points(metab.ox.median+metab.ox.sd, type = 'l', lwd=0.75, col='tomato3')
points(metab.ox.median-metab.ox.sd, type = 'l', lwd=0.75, col='tomato3')

#RB
metab.rb.median <- apply(t(metab.predicted)[names(which(metab.three.clusters=='RB')),], 2, median)
metab.rb.sd<- apply(t(metab.predicted)[names(which(metab.three.clusters=='RB')),], 2, sd)

lines(metab.rb.median, type = 'l', lwd=2, col='#50B547', ylim=c(-2,3))
points(metab.rb.median, lwd=1, pch=19, col='#50B547')
polygon(c(1:30, rev(1:30)), c(metab.rb.median-metab.rb.sd,rev(metab.rb.median+metab.rb.sd)), border = F, col=alpha('#50B547', 0.3))
points(metab.rb.median+metab.rb.sd, type = 'l', lwd=0.75, col='#50B547')
points(metab.rb.median-metab.rb.sd, type = 'l', lwd=0.75, col='#50B547')

#RC
metab.rc.median <- apply(t(metab.predicted)[names(which(metab.three.clusters=='RC')),], 2, median)
metab.rc.sd<- apply(t(metab.predicted)[names(which(metab.three.clusters=='RC')),], 2, sd)

lines(metab.rc.median, type = 'l', lwd=2, col='cornflowerblue', ylim=c(-2,3))
points(metab.rc.median, lwd=1, pch=19, col='cornflowerblue')
polygon(c(1:30, rev(1:30)), c(metab.rc.median-metab.rc.sd,rev(metab.rc.median+metab.rc.sd)), border = F, col=alpha('cornflowerblue', 0.3))
points(metab.rc.median+metab.rc.sd, type = 'l', lwd=0.75, col='cornflowerblue')
points(metab.rc.median-metab.rc.sd, type = 'l', lwd=0.75, col='cornflowerblue')

dev.off()