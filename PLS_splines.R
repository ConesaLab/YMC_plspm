#################################################################################
## All pairwise integration with spline models
##
## 1. RNA - metabolites
## 2. Histones - RNA -> Histones need the spline in DE gene-associated HMs*
#################################################################################
library(mixOmics)

load('../../Thesis/R/plspm_matrices.RData')
load('../../Thesis/R/New_tmm_clusters.RData')
load('HMs_rna_MOREassociation.RData')

source("calculate_R2_Q2.R")

library(scales)

#############################
## Metabolites - RNA

metabolites <- t(cbind(Ox.metab, Rb.metab, Rc.metab))
rna <- t(cbind(Ox.rna, Rb.rna, Rc.rna))

myresult <- pls(t(metabolites), t(rna), ncomp = 3, mode = "regression", scale = TRUE) # scale = TRUE means autoscaling

myresult$prop_expl_var
gene.col <- rna.CLUSTERS.scaled$cluster
gene.col[which(gene.col==2)] <- "tomato3"
gene.col[which(gene.col=='3')] <- "#50B547"
gene.col[which(gene.col=='1')] <- "cornflowerblue"

gene.col <- gene.col[rownames(rna)]


####################################################

library(wordcloud)

metab.col <- c(rep(alpha('darkred',0.7), dim(Ox.metab)[2]), rep(alpha('darkgreen',0.7), dim(Rb.metab)[2]), rep(alpha('darkblue',0.7), dim(Rc.metab)[2]))


jpeg(filename = 'metab_RNA_PLS_2.jpg', width = 5, height = 5, units = 'in', res = 500)

par(mar=c(5,5,5,4)+.1, font.axis = 2, font.lab = 2)

plot(scale(myresult$loadings$Y, center = F, scale = T), col = alpha(gene.col, 0.4), pch = 19, cex = 1.2, 
     frame.plot = F, cex.lab = 1, xlab = 'C1 (53.7% & 38.1%)', ylab = 'C2 (41% & 27.3%)', 
     cex.main = 1, cex.axis = 1, main = 'Metabolites - RNA-seq PLS', xlim = c(-2, 2), ylim = c(-2,2))
points(scale(myresult$loadings$X, center = F, scale = T), col = metab.col, pch = 17, cex=1.5)
text(scale(myresult$loadings$X, center = F, scale = T)[,1], scale(myresult$loadings$X, center = F, scale = T)[,2], 
         rownames(myresult$loadings$X), pos = 3, cex = 1.2)

axis(1, lwd = 4, labels = F)
axis(2, lwd = 4, labels = F)

dev.off()
myresult$explained_variance

source('auxiliary_functions.R')

r.square <- tune.comp(t(metabolites), t(rna), factor = c(1:30), fold = nrow(t(metabolites)),ncomp = 3, rep = 1, option = "R")
q.square <- tune.comp(t(metabolites), t(rna), factor = c(1:30), fold = nrow(t(metabolites)),ncomp = 3, rep = 1, option = "Q")




## Subset genes / metabolites by coordinates in the PLS components

comp.trans <- scale(myresult$loadings$Y, center = F, scale = T)

genes.in.ox.rb <- rownames(comp.trans)[comp.trans[,1]>1&comp.trans[,2]>-0.5&comp.trans[,2]<0.6]
genes.in.rb.rc <- rownames(comp.trans)[comp.trans[,2]<(-1)&comp.trans[,1]>-0.4&comp.trans[,1]<0.4]
genes.in.rc.ox <- rownames(comp.trans)[comp.trans[,1]<(-0.5)&comp.trans[,2]>0.4&comp.trans[,2]<1.2]

write.table(genes.in.ox.rb, file = '../Data/ox_rb_genes.txt', quote = F, col.names = F, row.names = F)
write.table(genes.in.rb.rc, file = '../Data/rb_rc_genes.txt', quote = F, col.names = F, row.names = F)
write.table(genes.in.rc.ox, file = '../Data/rc_ox_genes.txt', quote = F, col.names = F, row.names = F)


# Plot doughnut

plot(scale(myresult$loadings$Y, center = F, scale = T), col = alpha(gene.col, 0.4), pch = 19, cex = 1.2, 
     frame.plot = F, cex.lab = 1.5, xlab = 'C1 (53.7% & 38.1%)', ylab = 'C2 (41% & 27.3%)', 
     cex.main = 1.8, cex.axis = 1.3, main = 'Metabolites - Gene Expression PLS', xlim = c(-2, 2), ylim = c(-2,2))

points(scale(myresult$loadings$Y, center = F, scale = T)[genes.in.ox.rb,], col=alpha('black', 0.2), pch=19)
points(scale(myresult$loadings$Y, center = F, scale = T)[genes.in.rb.rc,], col=alpha('black', 0.2), pch=19)
points(scale(myresult$loadings$Y, center = F, scale = T)[genes.in.rc.ox,], col=alpha('black', 0.2), pch=19)

text(scale(myresult$loadings$X, center = F, scale = T), rownames(myresult$loadings$X), pos = 3, cex = 1.2)

# load genomic locations
saccer.gtf <- read.table('../../Thesis/Datasets/Saccharomyces_cerevisiae.R64-1-1.110.gtf', sep='\t')
saccer.gtf <- saccer.gtf[saccer.gtf$V3=='transcript',]
saccer.gtf$V9 <- unlist(lapply(saccer.gtf$V9, function(x){gsub(pattern = 'gene_id ', replacement = '', x = unlist(strsplit(x, split = ';'))[1])}))

ox.rb.chroms <- table(saccer.gtf[saccer.gtf$V9%in%genes.in.ox.rb,1])
rb.rc.chroms <- table(saccer.gtf[saccer.gtf$V9%in%genes.in.rb.rc,1])
rc.ox.chroms <- table(saccer.gtf[saccer.gtf$V9%in%genes.in.rc.ox,1])
table(saccer.gtf$V1)

ox.rb.chroms/table(saccer.gtf$V1)[names(ox.rb.chroms)]
rb.rc.chroms/table(saccer.gtf$V1)[names(rb.rc.chroms)]
rc.ox.chroms/table(saccer.gtf$V1)[names(rc.ox.chroms)]

for (chr in names(ox.rb.chroms)){
  saccer.gtf.chr <- saccer.gtf[saccer.gtf$V1==chr,]
  saccer.gtf.chr <- saccer.gtf.chr[order(saccer.gtf.chr$V4, decreasing = F),]
  
  gene.cluster <- c()
  
  for (gene in saccer.gtf.chr$V9){
    if (gene%in%genes.in.ox.rb){
      gene.cluster <- c(gene.cluster, 1)
    }
    else if (gene%in%genes.in.rb.rc){
      gene.cluster <- c(gene.cluster, 2)
    }
    else if (gene%in%genes.in.rc.ox){
      gene.cluster <- c(gene.cluster, 3)
    }
    else{
      gene.cluster <- c(gene.cluster, 0)
    }
  }
  
  names(gene.cluster) <- saccer.gtf.chr$V9
  
  write.table(gene.cluster, file = paste0('../Datasets/chr_order_datasets/chr_',chr,'_genes_betweenClusters.txt'), col.names = F, quote=F)

}


## Create heatmaps of the atac profile for the genes encoded

load('../../Thesis/R_2/Data/plspm_datasets_lvs.RData')

atac.genes <- c(colnames(ATAC.hoc), colnames(ATAC.loc))
k18.genes <- substr(c(colnames(k18.hoc), colnames(k18.loc)), start = 5, stop=20)
k9.genes <- c(colnames(k9.hoc), substr(colnames(k9.loc), start = 7, stop=20))

write.table(t(cbind(ATAC.hoc,ATAC.loc))[atac.genes%in%genes.in.ox.rb,], file = '../R/atac_in_oc-rb_transition.txt', quote=F, col.names = T, row.names = T, sep='\t')
write.table(t(cbind(ATAC.hoc,ATAC.loc))[atac.genes%in%genes.in.rb.rc,], file = '../R/atac_in_rb-rc_transition.txt', quote=F, col.names = T, row.names = T, sep='\t')
write.table(t(cbind(ATAC.hoc,ATAC.loc))[atac.genes%in%genes.in.rc.ox,], file = '../R/atac_in_rc-ox_transition.txt', quote=F, col.names = T, row.names = T, sep='\t')



## Write the atac, k18ac and k9ac profiles to  "see" the chromatin dynamics across the cycle

load('../../Thesis/R_2/Data/all_plspm_datasets.RData')

genes.order <- read.table('~/Downloads/order_genes.txt')$V1

write.table(t(rna.predicted[,genes.order]), file = '../R/fullrna_in_oc-rb_transition.txt', quote=F, col.names = T, row.names = T, sep='\t')
write.table(t(atac.predicted[,genes.order]), file = '../R/fullatac_in_oc-rb_transition.txt', quote=F, col.names = T, row.names = T, sep='\t')
write.table(t(k18.predicted[,genes.order[genes.order%in%colnames(k18.predicted)]]), file = '../R/fullk18_in_oc-rb_transition.txt', quote=F, col.names = T, row.names = T, sep='\t')
write.table(t(k9.predicted[,genes.order]), file = '../R/fullk9_in_oc-rb_transition.txt', quote=F, col.names = T, row.names = T, sep='\t')

common.genes <- intersect(intersect(intersect(colnames(rna.predicted), atac.genes), k18.genes), k9.genes)
write.table(t(rna.predicted[,common.genes][,genes.order]), file = '../R/fullrna__common_ordered.txt', quote=F, col.names = T, row.names = T, sep='\t')
write.table(t(atac.predicted[,common.genes][,genes.order]), file = '../R/fullatac_common_ordered.txt', quote=F, col.names = T, row.names = T, sep='\t')
write.table(t(k18.predicted[,common.genes][,genes.order]), file = '../R/fullk18_common_ordered.txt', quote=F, col.names = T, row.names = T, sep='\t')
write.table(t(k9.predicted[,common.genes][,genes.order]), file = '../R/fullk9_common_ordered.txt', quote=F, col.names = T, row.names = T, sep='\t')


## Test correlation 
library(scales)
gene.col <- rna.CLUSTERS.scaled$cluster
gene.col[which(gene.col==2)] <- "tomato3"
gene.col[which(gene.col=='3')] <- "#50B547"
gene.col[which(gene.col=='1')] <- "cornflowerblue"
gene.col <- gene.col[genes.order]

pdf('../Figures/heatmaps_chromatin/rna_omic_correlation_heatmaps.pdf', width = 9, height = 5)
par(mfrow=c(1,3))
plot(rna.predicted[,common.genes][,genes.order], atac.predicted[,common.genes][,genes.order], frame.plot=F, pch=19, col=alpha(gene.col, 0.3), xlab='RNA', ylab='ATAC')
cor.test(rna.predicted[,common.genes][,genes.order], atac.predicted[,common.genes][,genes.order])
plot(rna.predicted[,common.genes][,genes.order], k18.predicted[,common.genes][,genes.order], frame.plot=F, pch=19, col=alpha(gene.col, 0.3), xlab='RNA', ylab='K18')
cor.test(rna.predicted[,common.genes][,genes.order], k18.predicted[,common.genes][,genes.order])
plot(rna.predicted[,common.genes][,genes.order], k9.predicted[,common.genes][,genes.order], frame.plot=F, pch=19, col=alpha(gene.col, 0.3), xlab='RNA', ylab='K9')
cor.test(rna.predicted[,common.genes][,genes.order], k9.predicted[,common.genes][,genes.order])
dev.off()

#########
## which phases are the common
load('../../Thesis/R_2/Data/net_clusters.RData')

table(rna.CLUSTERS.scaled$cluster[genes.order])

#############################
## Metabolites - NET

net.two.clusters$cluster

metabolites <- t(cbind(Ox.metab, Rb.metab, Rc.metab))
net <- t(cbind(Ox.net, Rb.net, Rc.net))

net <- net[c(names(which(net.two.clusters$cluster==1)), names(which(net.two.clusters$cluster==2))),]

myresult <- pls(t(metabolites), t(net), ncomp = 3, mode = "regression", scale = TRUE) # scale = TRUE means autoscaling


net.col <- c(rep('purple', length(which(net.two.clusters$cluster==1))), rep('orange', length(which(net.two.clusters$cluster==2))))

metab.col <- c(rep(alpha('darkred',0.7), dim(Ox.metab)[2]), rep(alpha('darkgreen',0.7), dim(Rb.metab)[2]), rep(alpha('darkblue',0.7), dim(Rc.metab)[2]))


jpeg(filename = 'metab_NET_PLS_2.jpg', width = 5, height = 5, units = 'in', res = 500)

par(mar=c(5,5,5,4)+.1, font.axis = 2, font.lab = 2)

plot(scale(myresult$loadings$Y, center = F, scale = T), col = alpha(net.col, 0.4), pch = 19, cex = 1.2, 
     frame.plot = F, cex.lab = 1, xlab = 'C1 (54.8% & 38.4%)', ylab = 'C2 (34.5% & 27.1%)', 
     cex.main = 1, cex.axis = 1, main = 'Metabolites - NET-seq PLS', xlim = c(-2, 2), ylim = c(-2,2))
points(scale(myresult$loadings$X, center = F, scale = T), col = metab.col, pch = 17, cex = 1.5)
text(scale(myresult$loadings$X, center = F, scale = T), rownames(myresult$loadings$X), pos = 3, cex = 1.2)

axis(1, lwd = 4, labels = F)
axis(2, lwd = 4, labels = F)

dev.off()





metabolites <- t(cbind(Ox.metab, Rb.metab, Rc.metab))
net <- t(cbind(Ox.net, Rb.net, Rc.net))


net.col <- c(rep('tomato3', dim(Ox.net)[2]), rep('#50B547', dim(Rb.net)[2]), rep('cornflowerblue', dim(Rc.net)[2]))

metab.col <- c(rep('purple', dim(Ox.metab)[2]), rep('darkorange', dim(Rb.metab)[2]), rep('yellow2', dim(Rc.metab)[2]))



jpeg(filename = 'metab_NET_PLS_2.jpg', width = 10, height = 9, units = 'in', res = 500)

par(mar=c(5,5,5,4)+.1, font.axis = 2, font.lab = 2)

plot(scale(myresult$loadings$Y, center = F, scale = T), col = alpha(net.col, 0.4), pch = 19, cex = 1.2, 
     frame.plot = F, cex.lab = 1.5, xlab = 'C1 (54.8% & 38.4%)', ylab = 'C2 (34.5% & 27.1%)', 
     cex.main = 1.8, cex.axis = 1.3, main = 'Metabolites - Gene Expression PLS', xlim = c(-2, 2), ylim = c(-2,2))
points(scale(myresult$loadings$X, center = F, scale = T), col = metab.col)
points(scale(myresult$loadings$X, center = F, scale = T), col = metab.col, pch = 19, cex = 1.2)
text(scale(myresult$loadings$X, center = F, scale = T), rownames(myresult$loadings$X), pos = 3, cex = 1.2)

axis(1, lwd = 4, labels = F)
axis(2, lwd = 4, labels = F)

dev.off()

myresult$explained_variance

source('auxiliary_functions.R')

r.square <- tune.comp(t(metabolites), t(net), factor = c(1:30), fold = nrow(t(metabolites)),ncomp = 3, rep = 1, option = "R")
q.square <- tune.comp(t(metabolites), t(net), factor = c(1:30), fold = nrow(t(metabolites)),ncomp = 3, rep = 1, option = "Q")

#############################
## HMs - RNA

#############################
##Get significant associations


library('mixOmics')

source('MORE/auxFunctions.R')
source('MORE/ComputeGLM_function.R')
source('MORE/MORE_GLM.R')

gene.colours <- rna.CLUSTERS.scaled$cluster
gene.colours[which(gene.colours==2)] <- 'tomato3'
gene.colours[which(gene.colours==1)] <- 'cornflowerblue'
gene.colours[which(gene.colours==3)] <- '#50B547'

rna <- t(cbind(Ox.rna, Rb.rna, Rc.rna))

H3k18ac <- t(cbind(Ox.k18, Rb.k18, Rc.k18))
k18.rna <- rownames(rna)[which(rownames(rna)%in%rownames(H3k18ac))]
k18.association <- as.data.frame(cbind(k18.rna, paste('H3k18ac', k18.rna, sep = '_')))
rownames(H3k18ac) <- paste('H3k18ac', rownames(H3k18ac), sep = '_')
H3k9ac <- t(cbind(Ox.k9, Rb.k9, Rc.k9))
k9.rna <- rownames(rna)[which(rownames(rna)%in%rownames(H3k9ac))]
k9.association <- as.data.frame(cbind(k9.rna, paste('H3k9ac', k9.rna, sep = '_')))
rownames(H3k9ac) <- paste('H3k9ac', rownames(H3k9ac), sep='_')


associations <- list('H3K9ac'=k9.association, 'H3K18ac'=k18.association)
omics <- list('H3K9ac'=H3k9ac, 'H3K18ac'=H3k18ac)

GLMresults <- GetGLM(GeneExpression = rna,
                     associations = associations,
                     data.omics = omics,
                     edesign = NULL, Res.df = 8,
                     epsilon = 0.00001, alfa = 0.05,
                     MT.adjust = 'fdr', family = gaussian(),
                     elasticnet = 2, center = F, scale = F,
                     stepwise = "two.ways.backward",
                     interactions.exp = T, interactions.reg = 1,
                     correlation = 0.95, action = 'mean',
                     cont.var = NULL, min.obs = 10)

k9Genes <- c()
k18Genes <- c()

for (gene in GLMresults$ResultsPerGene){
  if ('significantRegulators'%in%names(gene)){
    if (length(grep('H3k9ac', gene$significantRegulators)) == 1){
      k9Genes <- c(k9Genes, gene$allRegulators$gene)
    }
    if (length(grep('H3k18ac', gene$significantRegulators)) == 1){
      k18Genes <- c(k18Genes, gene$allRegulators$gene)
    }
  }
}

save(k9Genes, k18Genes, file = 'hist_rna_MORE.RData')


#############################
load('hist_rna_MORE.RData')
## H3K18ac
H3k18ac <- t(cbind(Ox.k18, Rb.k18, Rc.k18))
rna <- t(cbind(Ox.rna, Rb.rna, Rc.rna))
myresult <- pls(t(H3k18ac), t(rna), ncomp = 3, mode = "regression", scale = TRUE) # scale = TRUE means autoscaling

gene.col <- rna.CLUSTERS.scaled$cluster
gene.col[which(gene.col==2)] <- "tomato3"
gene.col[which(gene.col=='3')] <- "#50B547"
gene.col[which(gene.col=='1')] <- "cornflowerblue"

gene.col <- gene.col[rownames(rna)]

jpeg(filename = 'rna_k18_PLS.jpg', width = 14, height = 8, units = 'in', res = 500)

par(mfrow = c(1,2))
plot(scale(myresult$loadings$Y, center = F, scale = T), xlim = c(-2, 2), ylim = c(-2,2), col = gene.col, main = 'Genes', frame.plot = F)
points(scale(myresult$loadings$Y, center = F, scale = T), col = gene.col, pch = 19, cex = 1.2)
points(scale(myresult$loadings$Y, center = F, scale = T)[k18Genes[which(k18Genes%in%(rownames(H3k18ac)))],], col = 'gray36', cex = 1.2)

axis(1, lwd = 4, labels = F)
axis(2, lwd = 4, labels = F)

##
k18.col <- c(rep('tomato3', dim(Ox.k18)[2]), rep('#50B547', dim(Rb.k18)[2]), rep('cornflowerblue', dim(Rc.k18)[2]))

plot(scale(myresult$loadings$X, center = F, scale = T), xlim = c(-2, 2), ylim = c(-2,2), col = k18.col, main = 'H3K18ac', frame.plot = F)
points(scale(myresult$loadings$X, center = F, scale = T), col = k18.col, pch = 19, cex = 1.2)
points(scale(myresult$loadings$X, center = F, scale = T)[k18Genes[which(k18Genes%in%(rownames(H3k18ac)))],], col = 'gray36', cex = 1.2)

axis(1, lwd = 4, labels = F)
axis(2, lwd = 4, labels = F)

mtext("H3K18ac", side = 3, line = -2, outer = TRUE, cex = 2)

dev.off()
myresult$explained_variance

source('auxiliary_functions.R')

r.square <- tune.comp(t(H3k18ac), t(rna), factor = c(1:30), fold = nrow(t(H3k18ac)),ncomp = 3, rep = 1, option = "R")
q.square <- tune.comp(t(H3k18ac), t(rna), factor = c(1:30), fold = nrow(t(H3k18ac)),ncomp = 3, rep = 1, option = "Q")

## H3K9ac
H3k9ac <- t(cbind(Ox.k9, Rb.k9, Rc.k9))
rna <- t(cbind(Ox.rna, Rb.rna, Rc.rna))
myresult <- pls(t(H3k9ac), t(rna), ncomp = 3, mode = "regression", scale = TRUE) # scale = TRUE means autoscaling

gene.col <- rna.CLUSTERS.scaled$cluster
gene.col[which(gene.col==2)] <- "tomato3"
gene.col[which(gene.col=='3')] <- "#50B547"
gene.col[which(gene.col=='1')] <- "cornflowerblue"

gene.col <- gene.col[rownames(rna)]

jpeg(filename = 'rna_k9_PLS.jpg', width = 14, height = 8, units = 'in', res = 500)
par(mfrow = c(1,2))
plot(scale(myresult$loadings$Y, center = F, scale = T), xlim = c(-2, 2), ylim = c(-2,2), col = gene.col, main = 'Genes', frame.plot = F)
points(scale(myresult$loadings$Y, center = F, scale = T), col = gene.col, pch = 19, cex = 1.2)
points(scale(myresult$loadings$Y, center = F, scale = T)[k9Genes,], col = 'gray36', cex = 1.2)

##
k9.col <- c(rep('tomato3', dim(Ox.k9)[2]), rep('#50B547', dim(Rb.k9)[2]), rep('cornflowerblue', dim(Rc.k9)[2]))

plot(scale(myresult$loadings$X, center = F, scale = T), xlim = c(-2, 2), ylim = c(-2,2), col = k9.col, main = 'H3K9ac', frame.plot = F)
points(scale(myresult$loadings$X, center = F, scale = T), col = k9.col, pch = 19, cex = 1.2)
points(scale(myresult$loadings$X, center = F, scale = T)[k9Genes,], col = 'gray36', cex = 1.2)
mtext("H3K9ac", side = 3, line = -2, outer = TRUE, cex = 2)

dev.off()

source('auxiliary_functions.R')

r.square <- tune.comp(t(H3k9ac), t(rna), factor = c(1:30), fold = nrow(t(H3k9ac)),ncomp = 3, rep = 1, option = "R")
q.square <- tune.comp(t(H3k9ac), t(rna), factor = c(1:30), fold = nrow(t(H3k9ac)),ncomp = 3, rep = 1, option = "Q")

r.square
q.square
