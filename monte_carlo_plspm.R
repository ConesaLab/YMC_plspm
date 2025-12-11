####################################################
## automate iteration of plspm models for 

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


load('plspm_model2_results_240507.RData')








###############################################################

load('Data/plspm_datasets_lvs.RData')

source('random_model.R')
source('plspm_automate.R')

#Map displaying which element of the model corresponds to which matrix of indicators
indicators.match.list <- list('Rc.metabs1' = c('metab.rc1'), 'Rc.metabs2' = c('metab.rc2'), 'HOC.hms' = c('k9.hoc', 'k18.hoc'),
                              'HOC.net' = c('net.hoc'), 'HOC.atac' = c('ATAC.hoc'), 'Ox.genes' = c('RNA.ox'),
                              'Ox.metabs1' = c('metab.ox1'), 'Rb.genes' = c('RNA.rb'),
                              'Rb.metabs' = c('metab.rb'), 'LOC.hms' = c('k9.loc', 'k18.loc'),
                              'LOC.net' = c('net.loc'), 'LOC.atac' = c('ATAC.loc'),'Rc.genes' = c('RNA.rc'))

# List of indicators with all original matrices separated instead of bound, this applies to HMs, which are separated

indicators.list.sep = list(RNA.rb, metab.rb, k9.loc, k18.loc, ATAC.loc, net.loc, RNA.rc, 
                           metab.rc1, metab.rc2, k9.hoc, k18.hoc, ATAC.hoc, net.hoc, RNA.ox, metab.ox1)

names(indicators.list.sep) <- c('RNA.rb', 'metab.rb', 'k9.loc', 'k18.loc', 'ATAC.loc', 'net.loc', 'RNA.rc', 
                                'metab.rc1', 'metab.rc2', 'k9.hoc', 'k18.hoc', 'ATAC.hoc', 'net.hoc', 'RNA.ox', 'metab.ox1')

#Randomize the elements of the model
model.elements <- c('Rc.metabs1', 'Rc.metabs2', 'HOC.hms', 'HOC.net', 'HOC.atac', 'Ox.genes', 
                    'Ox.metabs1', 'Rb.genes', 'Rb.metabs', 'LOC.hms', 'LOC.net', 'LOC.atac',  'Rc.genes')




total.gof <- c()
rels.list <- list()
good.models <- list()
bad.models <- list()
latent_modes = rep("A", 13)
number.rels <- c()

bb = 0; zz = 0

#Iterate n number of times to generate n number of random models and retrieve rels and gof stats
for (y in 1:1500){
  print(paste('Random iteration number:', y))
  uuu <- randomize.model(model.elements = model.elements, indicators.match.list = indicators.match.list, indicators.list.sep = indicators.list.sep)
  flag <- TRUE
  
  number.rels <- c(number.rels, sum(rowSums(uuu$random.model)))
  
  tryCatch({
    lv.fitted.plspm <- plspm.fit.lvs(indicators_matrix = list.cbind(uuu$indicators.list), model = uuu$random.model, 
                                     blocks = indicators_block(uuu$indicators.list), 
                                     modes = latent_modes, indicators_list = uuu$indicators.list)
    }, error = function(e){
                                       if (e$message == 'One block lacks indicators'){
                                         flag<<-FALSE
                                         print('One block lacks indicators')
                                       }
                                       else if (test_null_weights(weights, specs) == TRUE){
                                         flag<<-FALSE
                                         print('Maximum number of iterations')
                                       }
                                       else{
                                         stop(e)
                                       }
                                     }
  )
  
  if (!flag){
    print('Error: this iteration could not be processed - 1')
    zz = zz + 1
    bad.models[[zz]] = uuu
    next
  }
  
  tryCatch({
    rel.fitted.plspm <- plspm.rel.fit(plspm.obj = lv.fitted.plspm$plspm, indicators_matrix = list.cbind(lv.fitted.plspm$indicators), 
                                      model = uuu$random.model, blocks = indicators_block(lv.fitted.plspm$indicators), modes = latent_modes, 
                                      indicators_list = lv.fitted.plspm$indicators)}, error = function(e) {
                                        if (e$message == "missing value where TRUE/FALSE needed"){
                                          flag<<-FALSE
                                          print('One LV has no relationships')
                                        }
                                      }
  )
  if (!flag){
    print('Error: this iteration could not be processed - 2') 
    zz = zz + 1
    bad.models[[zz]] = uuu
    next
  } 
  
  for (n in 1:length(rel.fitted.plspm$inner_model)){
    omic <- names(rel.fitted.plspm$inner_model)[n]
    rels <- rel.fitted.plspm$inner_model[[n]]
    if (!omic%in%names(rels.list)){
      rels.list[[omic]] <- list()
    }
    for (rel.elem in rownames(rel.fitted.plspm$inner_model[[omic]])[-1]){
      if (rel.elem%in%names(rels.list[[omic]])){
        rels.list[[omic]][[rel.elem]] <- c(rels.list[[omic]][[rel.elem]], 
                                           rel.fitted.plspm$inner_model[[omic]][rel.elem,4])
      }
      else{
        rels.list[[omic]][[rel.elem]] <- c(rel.fitted.plspm$inner_model[[omic]][rel.elem,4])
      }
    }
  }
  total.gof <- c(total.gof, rel.fitted.plspm$gof)
  if (rel.fitted.plspm$gof > 0.7353491){
    bb = bb + 1
    good.models[[bb]] = rel.fitted.plspm
  }
}

save(good.models, bad.models, total.gof, rels.list, number.rels, file = 'divideLVs_thousand_random_iter11.RData')


hist(number.rels)
length(bad.models)
length(good.models)



