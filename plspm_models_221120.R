###########################################################
## Create PLSPM models

####################################################
## Model

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
