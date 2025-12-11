###########################################################
## Create PLSPM models

## Model 1: Metab - K9/K18 - NET - GE - Metab 
## Model 2: Metab - K9/K18 - ATAC / NET - GE - Metab



####################################################
## Model 1

load('Data/plspm_datasets_lvs.RData')

source('plspm_automate.R')

Rb.genes <- c(0,0,0,0,0,0,0,0,0,0,0)
Rb.metabs <- c(1,0,0,0,0,0,0,0,0,0,0)
LOC.hms <- c(0,1,0,0,0,0,0,0,0,0,0)
LOC.NET <- c(0,0,1,0,0,0,0,0,0,0,0)
Rc.genes <- c(0,0,0,1,0,0,0,0,0,0,0)
Rc.metabs1 <- c(0,0,0,0,1,0,0,0,0,0,0)
Rc.metabs2 <- c(0,0,0,0,1,0,0,0,0,0,0)
HOC.hms <- c(0,0,0,0,0,1,1,0,0,0,0)
HOC.NET <- c(0,0,0,0,0,0,0,1,0,0,0)
Ox.genes <- c(0,0,0,0,0,0,0,0,1,0,0)
Ox.metabs1 <- c(0,0,0,0,0,0,0,0,0,1,0)


ymc.model = rbind(Rb.genes, Rb.metabs, LOC.hms, LOC.NET, Rc.genes, 
      Rc.metabs1, Rc.metabs2, HOC.hms, HOC.NET, Ox.genes, Ox.metabs1)

colnames(ymc.model) = rownames(ymc.model)

innerplot(ymc.model)

indicators_matrix <- cbind(RNA.rb, metab.rb, k9.loc, k18.loc, net.loc, RNA.rc, metab.rc1, metab.rc2,
                           k9.hoc, k18.hoc, net.hoc, RNA.ox, metab.ox1)

indicators_blocks <- indicators_block(blocks_list = list(RNA.rb, metab.rb, cbind(k9.loc, k18.loc), net.loc, 
                                                         RNA.rc, metab.rc1, metab.rc2, cbind(k9.hoc, k18.hoc),
                                                         net.hoc, RNA.ox, metab.ox1))

indicators_list = list(RNA.rb, metab.rb, cbind(k9.loc, k18.loc), net.loc, 
                       RNA.rc, metab.rc1, metab.rc2, cbind(k9.hoc, k18.hoc),
                       net.hoc, RNA.ox, metab.ox1)

latent_modes = rep("A", 11)

lv.fitted.plspm <- plspm.fit.lvs(indicators_matrix = indicators_matrix, model = ymc.model, blocks = indicators_blocks, 
                                 modes = latent_modes, indicators_list = indicators_list)

plot(lv.fitted.plspm$plspm)

rel.fitted.plspm <- plspm.rel.fit(plspm.obj = lv.fitted.plspm$plspm, indicators_matrix = list.cbind(lv.fitted.plspm$indicators), 
                                  model = ymc.model, blocks = indicators_block(lv.fitted.plspm$indicators), modes = latent_modes, 
                                  indicators_list = lv.fitted.plspm$indicators)

plot(rel.fitted.plspm)

rel.fitted.plspm$gof

rel.fitted.plspm$path_coefs

plot.omic(RNA.rb, "Rb.genes", rel.fitted.plspm)
plot.omic(metab.rb, "Rb.metabs", rel.fitted.plspm)
plot.omic(cbind(k9.loc, k18.loc), "LOC.hms", rel.fitted.plspm)

plot(apply(cbind(k9.loc, k18.loc)[,!colnames(cbind(k9.loc, k18.loc))%in%as.character(subset(rel.fitted.plspm$outer_model, block == "LOC.hms")$name)],1,mean), type='l', lwd=1.5)

plot(apply(t(RNA.rb)[as.character(subset(rel.fitted.plspm$outer_model, block == "Rb.genes")$name),], 2, mean))
subset(rel.fitted.plspm$outer_model, block == "Rb.metabs")
subset(rel.fitted.plspm$outer_model, block == "Rc.metabs1")
subset(rel.fitted.plspm$outer_model, block == "Rc.metabs2")

####################################################
## Model2

load('Data/plspm_datasets_lvs.RData')

source('plspm_automate.R')

Rb.genes <- c(0,0,0,0,0,0,0,0,0,0,0,0,0)
Rb.metabs <- c(1,0,0,0,0,0,0,0,0,0,0,0,0)
LOC.hms <- c(0,1,0,0,0,0,0,0,0,0,0,0,0)
LOC.atac <- c(0,0,1,0,0,0,0,0,0,0,0,0,0)
LOC.net <- c(0,0,0,1,0,0,0,0,0,0,0,0,0)
Rc.genes <- c(0,0,0,0,1,0,0,0,0,0,0,0,0)
Rc.metabs1 <- c(0,0,0,0,0,1,0,0,0,0,0,0,0)
Rc.metabs2 <- c(0,0,0,0,0,1,0,0,0,0,0,0,0)
HOC.hms <- c(0,0,0,0,0,0,1,1,0,0,0,0,0)
HOC.atac <- c(0,0,0,0,0,0,0,0,1,0,0,0,0)
HOC.net <- c(0,0,0,0,0,0,0,0,0,1,0,0,0)
Ox.genes <- c(0,0,0,0,0,0,0,0,0,0,1,0,0)
Ox.metabs1 <- c(0,0,0,0,0,0,0,0,0,0,0,1,0)




ymc.model = rbind(Rb.genes, Rb.metabs, LOC.hms, LOC.atac, LOC.net, Rc.genes, 
                  Rc.metabs1, Rc.metabs2, HOC.hms, HOC.atac, HOC.net, Ox.genes, Ox.metabs1)

colnames(ymc.model) = rownames(ymc.model)

innerplot(ymc.model)

indicators_matrix <- cbind(RNA.rb, metab.rb, k9.loc, k18.loc, ATAC.loc, net.loc, RNA.rc, metab.rc1, metab.rc2,
                           k9.hoc, k18.hoc, ATAC.hoc, net.hoc, RNA.ox, metab.ox1)

indicators_blocks <- indicators_block(blocks_list = list(RNA.rb, metab.rb, cbind(k9.loc, k18.loc), ATAC.loc, net.loc,
                                                         RNA.rc, metab.rc1, metab.rc2, cbind(k9.hoc, k18.hoc),
                                                         ATAC.hoc, net.hoc, RNA.ox, metab.ox1))

indicators_list = list(RNA.rb, metab.rb, cbind(k9.loc, k18.loc), ATAC.loc, net.loc,
                       RNA.rc, metab.rc1, metab.rc2, cbind(k9.hoc, k18.hoc),
                       ATAC.hoc, net.hoc, RNA.ox, metab.ox1)

latent_modes = rep("A", 13)


lv.fitted.plspm <- plspm.fit.lvs(indicators_matrix = indicators_matrix, model = ymc.model, blocks = indicators_blocks, 
                                 modes = latent_modes, indicators_list = indicators_list)

plot(lv.fitted.plspm$plspm)

rel.fitted.plspm <- plspm.rel.fit(plspm.obj = lv.fitted.plspm$plspm, indicators_matrix = list.cbind(lv.fitted.plspm$indicators), 
                                  model = ymc.model, blocks = indicators_block(lv.fitted.plspm$indicators), modes = latent_modes, 
                                  indicators_list = lv.fitted.plspm$indicators)

plot(rel.fitted.plspm)

rel.fitted.plspm$gof

rel.fitted.plspm$path_coefs

subset(rel.fitted.plspm$outer_model, block == "Ox.metabs1")
subset(rel.fitted.plspm$outer_model, block == "Rb.metabs")
subset(rel.fitted.plspm$outer_model, block == "Rc.metabs1")
subset(rel.fitted.plspm$outer_model, block == "Rc.metabs2")


dim(subset(rel.fitted.plspm$outer_model, block == "HOC.atac"))


save(rel.fitted.plspm, file = 'plspm_model2_results_240507.RData')

####################################################
####################################################
####################################################

#########################################################################


used.loc.hms <- as.character(subset(rel.fitted.plspm$outer_model, block == "LOC.hms")$name)
not.used.loc.hms <- as.character(colnames(cbind(k9.loc, k18.loc))[!colnames(cbind(k9.loc, k18.loc))%in%used.loc.hms])



Rb.genes <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0)
Rb.metabs <- c(1,0,0,0,0,0,0,0,0,0,0,0,0,0)
LOC.hms1 <- c(0,1,0,0,0,0,0,0,0,0,0,0,0,0)
LOC.hms2 <- c(0,1,0,0,0,0,0,0,0,0,0,0,0,0)
LOC.atac <- c(0,0,1,1,0,0,0,0,0,0,0,0,0,0)
LOC.net <- c(0,0,0,0,1,0,0,0,0,0,0,0,0,0)
Rc.genes <- c(0,0,0,0,0,1,0,0,0,0,0,0,0,0)
Rc.metabs1 <- c(0,0,0,0,0,0,1,0,0,0,0,0,0,0)
Rc.metabs2 <- c(0,0,0,0,0,0,1,0,0,0,0,0,0,0)
HOC.hms <- c(0,0,0,0,0,0,0,1,1,0,0,0,0,0)
HOC.atac <- c(0,0,0,0,0,0,0,0,0,1,0,0,0,0)
HOC.net <- c(0,0,0,0,0,0,0,0,0,0,1,0,0,0)
Ox.genes <- c(0,0,0,0,0,0,0,0,0,0,0,1,0,0)
Ox.metabs1 <- c(0,0,0,0,0,0,0,0,0,0,0,0,1,0)




ymc.model = rbind(Rb.genes, Rb.metabs, LOC.hms1, LOC.hms2, LOC.atac, LOC.net, Rc.genes, 
                  Rc.metabs1, Rc.metabs2, HOC.hms, HOC.atac, HOC.net, Ox.genes, Ox.metabs1)

colnames(ymc.model) = rownames(ymc.model)

innerplot(ymc.model)

indicators_matrix <- cbind(RNA.rb, metab.rb, cbind(k9.loc, k18.loc)[,used.loc.hms], cbind(k9.loc, k18.loc)[,not.used.loc.hms], 
                           ATAC.loc, net.loc, RNA.rc, metab.rc1, metab.rc2,
                           k9.hoc, k18.hoc, ATAC.hoc, net.hoc, RNA.ox, metab.ox1)

indicators_blocks <- indicators_block(blocks_list = list(RNA.rb, metab.rb, cbind(k9.loc, k18.loc)[,used.loc.hms], cbind(k9.loc, k18.loc)[,not.used.loc.hms], 
                                                         ATAC.loc, net.loc, RNA.rc, metab.rc1, metab.rc2, cbind(k9.hoc, k18.hoc),
                                                         ATAC.hoc, net.hoc, RNA.ox, metab.ox1))

indicators_list = list(RNA.rb, metab.rb, cbind(k9.loc, k18.loc)[,used.loc.hms], cbind(k9.loc, k18.loc)[,not.used.loc.hms], 
                       ATAC.loc, net.loc, RNA.rc, metab.rc1, metab.rc2, cbind(k9.hoc, k18.hoc),
                       ATAC.hoc, net.hoc, RNA.ox, metab.ox1)

latent_modes = rep("A", 14)


lv.fitted.plspm <- plspm.fit.lvs(indicators_matrix = indicators_matrix, model = ymc.model, blocks = indicators_blocks, 
                                 modes = latent_modes, indicators_list = indicators_list)

plot(lv.fitted.plspm$plspm)

rel.fitted.plspm <- plspm.rel.fit(plspm.obj = lv.fitted.plspm$plspm, indicators_matrix = list.cbind(lv.fitted.plspm$indicators), 
                                  model = ymc.model, blocks = indicators_block(lv.fitted.plspm$indicators), modes = latent_modes, 
                                  indicators_list = lv.fitted.plspm$indicators)

plot(rel.fitted.plspm)

traitors.matrix <- detect.traitors(rel.fitted.plspm)

rel.fitted.plspm$gof

rel.fitted.plspm$path_coefs

save(rel.fitted.plspm, file = 'plspm_model2_results_230604.RData')

for (lv in c('Rb.genes', 'Rb.metabs', 'LOC.hms1', 'LOC.hms2', 'LOC.atac', 'LOC.net', 'Rc.genes', 
             'Rc.metabs1', 'Rc.metabs2', 'HOC.hms', 'HOC.atac', 'HOC.net', 'Ox.genes', 'Ox.metabs1')){
  outer.model.names <- as.character(subset(rel.fitted.plspm$outer_model, block == lv)$name)
  if(length(grep(pattern = 'hms', lv))>0){
    hm.outer.model.names <- gsub(pattern = 'k9_','',outer.model.names)
    hm.outer.model.names <- gsub(pattern = 'k18_','',hm.outer.model.names)
    write.table(hm.outer.model.names, file = paste0('lv_geneNames/',paste(lv, 'geneNames.txt', sep='_')), col.names = F, row.names = F, quote = F, sep = '\t')
  }
  else if(length(grep(pattern = 'metab', lv))>0){
    write.table(outer.model.names, file = paste0('lv_geneNames/',paste(lv, 'metabNames.txt', sep='_')), col.names = F, row.names = F, quote = F, sep = '\t')
  }
  else{
    write.table(outer.model.names, file = paste0('lv_geneNames/',paste(lv, 'geneNames.txt', sep='_')), col.names = F, row.names = F, quote = F, sep = '\t')
  }
  
}


loc.atac <- as.character(subset(rel.fitted.plspm$outer_model, block == 'LOC.atac')$name)
hoc.atac <- as.character(subset(rel.fitted.plspm$outer_model, block == 'HOC.atac')$name)

rb.genes.names <- as.character(subset(rel.fitted.plspm$outer_model, block == 'Rb.genes')$name)
rc.genes.names <- as.character(subset(rel.fitted.plspm$outer_model, block == 'Rc.genes')$name)
ox.genes.names <- as.character(subset(rel.fitted.plspm$outer_model, block == 'Ox.genes')$name)


barplot(c(length(intersect(ox.genes.names, loc.atac)), length(intersect(rb.genes.names, loc.atac)), 
        length(intersect(rc.genes.names, loc.atac)), length(setdiff(setdiff(setdiff(loc.atac, rb.genes.names), ox.genes.names), rc.genes.names))), 
        col=c('red3','green3','blue3', 'gray55'), main='LOC atac', names=c('Ox','Rb','Rc', 'Undetermined'))

barplot(c(length(intersect(ox.genes.names, hoc.atac)), length(intersect(rb.genes.names, hoc.atac)), 
          length(intersect(rc.genes.names, hoc.atac)), length(setdiff(setdiff(setdiff(loc.atac, rb.genes.names), ox.genes.names), rc.genes.names))), 
col=c('red3','green3','blue3', 'gray55'), main='HOC atac', names=c('Ox','Rb','Rc', 'Undetermined'))

#########################################################################


barplot(c(length(intersect(colnames(RNA.ox), colnames(ATAC.loc))), length(intersect(colnames(RNA.rb), colnames(ATAC.loc))), 
          length(intersect(colnames(RNA.rc), colnames(ATAC.loc))), length(setdiff(setdiff(setdiff(colnames(ATAC.loc), colnames(RNA.rb)), colnames(RNA.ox)), colnames(RNA.rc)))), 
        col=c('red3','green3','blue3', 'gray55'), main='LOC atac', names=c('Ox','Rb','Rc', 'Undetermined'))

barplot(c(length(intersect(colnames(RNA.ox), colnames(ATAC.hoc))), length(intersect(colnames(RNA.rb), colnames(ATAC.hoc))), 
          length(intersect(colnames(RNA.rc), colnames(ATAC.hoc))), length(setdiff(setdiff(setdiff(colnames(ATAC.hoc), colnames(RNA.rb)), colnames(RNA.ox)), colnames(RNA.rc)))), 
        col=c('red3','green3','blue3', 'gray55'), main='HOC atac', names=c('Ox','Rb','Rc', 'Undetermined'))



#########################################################################

pdf('traitor_elements_classification.pdf', width = 6, height = 6)
par(mfrow=c(3,3))
barplot(table(colnames(traitors.matrix[,c(3:16)])[apply(traitors.matrix[traitors.matrix$block=='Rb.genes',c(3:16)], 1, function(x){which.max(abs(x))})]), main = 'Rb.genes', las=2)
barplot(table(colnames(traitors.matrix[,c(3:16)])[apply(traitors.matrix[traitors.matrix$block=='LOC.hms1',c(3:16)], 1, function(x){which.max(abs(x))})]), main = 'LOC.hms1', las=2)
barplot(table(colnames(traitors.matrix[,c(3:16)])[apply(traitors.matrix[traitors.matrix$block=='LOC.hms2',c(3:16)], 1, function(x){which.max(abs(x))})]), main = 'LOC.hms2', las=2)
barplot(table(colnames(traitors.matrix[,c(3:16)])[apply(traitors.matrix[traitors.matrix$block=='LOC.atac',c(3:16)], 1, function(x){which.max(abs(x))})]), main = 'LOC.atac', las=2)
barplot(table(colnames(traitors.matrix[,c(3:16)])[apply(traitors.matrix[traitors.matrix$block=='LOC.net',c(3:16)], 1, function(x){which.max(abs(x))})]), main = 'LOC.net', las=2)
barplot(table(colnames(traitors.matrix[,c(3:16)])[apply(traitors.matrix[traitors.matrix$block=='Rc.genes',c(3:16)], 1, function(x){which.max(abs(x))})]), main = 'Rc.genes', las=2)
barplot(table(colnames(traitors.matrix[,c(3:16)])[apply(traitors.matrix[traitors.matrix$block=='Rc.metabs1',c(3:16)], 1, function(x){which.max(abs(x))})]), main = 'Rc.metabs1', las=2)
barplot(table(colnames(traitors.matrix[,c(3:16)])[apply(traitors.matrix[traitors.matrix$block=='Rc.metabs2',c(3:16)], 1, function(x){which.max(abs(x))})]), main = 'Rc.metabs2', las=2)
barplot(table(colnames(traitors.matrix[,c(3:16)])[apply(traitors.matrix[traitors.matrix$block=='HOC.hms',c(3:16)], 1, function(x){which.max(abs(x))})]), main = 'HOC.hms', las=2)
barplot(table(colnames(traitors.matrix[,c(3:16)])[apply(traitors.matrix[traitors.matrix$block=='HOC.atac',c(3:16)], 1, function(x){which.max(abs(x))})]), main = 'HOC.atac', las=2)
barplot(table(colnames(traitors.matrix[,c(3:16)])[apply(traitors.matrix[traitors.matrix$block=='HOC.net',c(3:16)], 1, function(x){which.max(abs(x))})]), main = 'HOC.net', las=2)
barplot(table(colnames(traitors.matrix[,c(3:16)])[apply(traitors.matrix[traitors.matrix$block=='Ox.genes',c(3:16)], 1, function(x){which.max(abs(x))})]), main = 'Ox.genes', las=2)
dev.off()

write.table(traitors.matrix, file = 'traitors_matrix.txt', col.names = T, row.names = T, quote = F, sep = '\t')

####################################################
####################################################
####################################################

## Model3

load('Data/plspm_datasets_lvs.RData')

source('plspm_automate.R')

Rb.genes <- c(0,0,0,0,0,0,0,0,0,0,0,0)
Rb.metabs <- c(1,0,0,0,0,0,0,0,0,0,0,0)
LOC.hms <- c(0,1,0,0,0,0,0,0,0,0,0,0)
LOC.net <- c(0,0,1,0,0,0,0,0,0,0,0,0)
Rc.genes <- c(0,0,0,1,0,0,0,0,0,0,0,0)
Rc.metabs1 <- c(0,0,0,0,1,0,0,0,0,0,0,0)
Rc.metabs2 <- c(0,0,0,0,1,0,0,0,0,0,0,0)
HOC.hms <- c(0,0,0,0,0,1,1,0,0,0,0,0)
HOC.atac <- c(0,0,0,0,0,1,1,0,0,0,0,0)
HOC.net <- c(0,0,0,0,0,0,0,1,1,0,0,0)
Ox.genes <- c(0,0,0,0,0,0,0,0,0,1,0,0)
Ox.metabs1 <- c(0,0,0,0,0,0,0,0,0,0,1,0)




ymc.model = rbind(Rb.genes, Rb.metabs, LOC.hms, LOC.net, Rc.genes, 
                  Rc.metabs1, Rc.metabs2, HOC.hms, HOC.atac, HOC.net, Ox.genes, Ox.metabs1)

colnames(ymc.model) = rownames(ymc.model)

innerplot(ymc.model)

indicators_matrix <- cbind(RNA.rb, metab.rb, k9.loc, k18.loc, net.loc, RNA.rc, metab.rc1, metab.rc2,
                           k9.hoc, k18.hoc, ATAC.hoc, net.hoc, RNA.ox, metab.ox1)

indicators_blocks <- indicators_block(blocks_list = list(RNA.rb, metab.rb, cbind(k9.loc, k18.loc), net.loc,
                                                         RNA.rc, metab.rc1, metab.rc2, cbind(k9.hoc, k18.hoc),
                                                         ATAC.hoc, net.hoc, RNA.ox, metab.ox1))

indicators_list = list(RNA.rb, metab.rb, cbind(k9.loc, k18.loc), net.loc,
                       RNA.rc, metab.rc1, metab.rc2, cbind(k9.hoc, k18.hoc),
                       ATAC.hoc, net.hoc, RNA.ox, metab.ox1)

latent_modes = rep("A", 12)


lv.fitted.plspm <- plspm.fit.lvs(indicators_matrix = indicators_matrix, model = ymc.model, blocks = indicators_blocks, 
                                 modes = latent_modes, indicators_list = indicators_list)

plot(lv.fitted.plspm$plspm)

rel.fitted.plspm <- plspm.rel.fit(plspm.obj = lv.fitted.plspm$plspm, indicators_matrix = list.cbind(lv.fitted.plspm$indicators), 
                                  model = ymc.model, blocks = indicators_block(lv.fitted.plspm$indicators), modes = latent_modes, 
                                  indicators_list = lv.fitted.plspm$indicators)

plot(rel.fitted.plspm)

rel.fitted.plspm$gof

rel.fitted.plspm$path_coefs

subset(rel.fitted.plspm$outer_model, block == "Ox.metabs1")
subset(rel.fitted.plspm$outer_model, block == "Rb.metabs")
subset(rel.fitted.plspm$outer_model, block == "Rc.metabs1")
subset(rel.fitted.plspm$outer_model, block == "Rc.metabs2")


dim(subset(rel.fitted.plspm$outer_model, block == "HOC.atac"))






plot.omic <- function(matrix, omic, plspm){
  plot(apply(t(matrix)[as.character(subset(plspm$outer_model, block == omic)$name),], 2, mean), type='l', lwd=1.5, frame.plot=F, ylab='')
}
