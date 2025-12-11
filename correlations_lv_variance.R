##################################################
## Calculate correlations with the mean of a matrix
## for each latent variable


load(file = 'plspm_model2_results_240507.RData')
load('Data/plspm_datasets_lvs.RData')


calc.cor <- function(x){
  x.mean <- apply(t(x), 2, mean)
  return(cor(x, x.mean)[,1])
}


Rb.genes.bef <- calc.cor(RNA.rb)
Rb.genes.af <- calc.cor(rel.fitted.plspm$data[,rel.fitted.plspm$outer_model$name[rel.fitted.plspm$outer_model$block=='Rb.genes']])

Rb.metab.bef <- calc.cor(metab.rb)
Rb.metab.af <- calc.cor(rel.fitted.plspm$data[,rel.fitted.plspm$outer_model$name[rel.fitted.plspm$outer_model$block=='Rb.metabs']])

loc.k18.bef <- calc.cor(k18.loc)
loc.k9.bef <- calc.cor(k9.loc)
loc.hms.af <- calc.cor(rel.fitted.plspm$data[,rel.fitted.plspm$outer_model$name[rel.fitted.plspm$outer_model$block=='LOC.hms']])

ATAC.loc.bef <- calc.cor(ATAC.loc)
ATAC.loc.af <- calc.cor(rel.fitted.plspm$data[,rel.fitted.plspm$outer_model$name[rel.fitted.plspm$outer_model$block=='LOC.atac']])

NET.loc.bef <- calc.cor(net.loc)
NET.loc.af <- calc.cor(rel.fitted.plspm$data[,rel.fitted.plspm$outer_model$name[rel.fitted.plspm$outer_model$block=='LOC.net']])

Rc.genes.bef <- calc.cor(RNA.rc)
Rc.genes.af <- calc.cor(rel.fitted.plspm$data[,rel.fitted.plspm$outer_model$name[rel.fitted.plspm$outer_model$block=='Rc.genes']])

Rc.metab2.bef <- calc.cor(metab.rc2)
Rc.metab2.af <- calc.cor(rel.fitted.plspm$data[,rel.fitted.plspm$outer_model$name[rel.fitted.plspm$outer_model$block=='Rc.metabs2']])

Rc.metab1.bef <- calc.cor(metab.rc1)
Rc.metab1.af <- calc.cor(rel.fitted.plspm$data[,rel.fitted.plspm$outer_model$name[rel.fitted.plspm$outer_model$block=='Rc.metabs1']])

hoc.k18.bef <- calc.cor(k18.hoc)
hoc.k9.bef <- calc.cor(k9.hoc)
hoc.hms.af <- calc.cor(rel.fitted.plspm$data[,rel.fitted.plspm$outer_model$name[rel.fitted.plspm$outer_model$block=='HOC.hms']])

ATAC.hoc.bef <- calc.cor(ATAC.hoc)
ATAC.hoc.af <- calc.cor(rel.fitted.plspm$data[,rel.fitted.plspm$outer_model$name[rel.fitted.plspm$outer_model$block=='HOC.atac']])

NET.hoc.bef <- calc.cor(net.hoc)
NET.hoc.af <- calc.cor(rel.fitted.plspm$data[,rel.fitted.plspm$outer_model$name[rel.fitted.plspm$outer_model$block=='HOC.net']])

Ox.genes.bef <- calc.cor(RNA.ox)
Ox.genes.af <- calc.cor(rel.fitted.plspm$data[,rel.fitted.plspm$outer_model$name[rel.fitted.plspm$outer_model$block=='Ox.genes']])

Ox.metab.bef <- calc.cor(metab.ox1)
Ox.metab.af <- calc.cor(rel.fitted.plspm$data[,rel.fitted.plspm$outer_model$name[rel.fitted.plspm$outer_model$block=='Ox.metabs1']])


boxplot(abs(Rb.genes.bef), abs(Rb.genes.af), abs(loc.k18.bef), abs(loc.k9.bef), abs(loc.hms.af), abs(ATAC.loc.bef), abs(ATAC.loc.af), abs(NET.loc.bef), abs(NET.loc.af),
        abs(Rc.genes.bef), abs(Rc.genes.af), abs(Rc.metab2.bef), abs(Rc.metab2.af), abs(Rc.metab1.bef), abs(Rc.metab1.af), abs(hoc.k18.bef), abs(hoc.k9.bef), abs(hoc.hms.af),
        abs(ATAC.hoc.bef), abs(ATAC.hoc.af), abs(NET.hoc.bef), abs(NET.hoc.af), abs(Ox.genes.bef), abs(Ox.genes.af), abs(Ox.metab.bef), abs(Ox.metab.af), frame.plot=F, 
        col=c(rep()))

library(scales)

pdf('figures_heatmaps/correlations_lvs_befaf.pdf', width = 6, height = 5)
boxplot(abs(Rb.genes.bef), abs(Rb.metab.bef), c(abs(loc.k18.bef), abs(loc.k9.bef)), abs(ATAC.loc.bef), abs(NET.loc.bef), abs(Rc.genes.bef), abs(Rc.metab2.bef),
        abs(Rc.metab1.bef), c(abs(hoc.k18.bef), abs(hoc.k9.bef)), abs(ATAC.hoc.bef), abs(NET.hoc.bef), abs(Ox.genes.bef), abs(Ox.metab.bef), abs(Rb.genes.af), 
        abs(Rb.metab.af), abs(loc.hms.af), abs(ATAC.loc.af), abs(NET.loc.af), abs(Rc.genes.af), abs(Rc.metab2.af),
        abs(Rc.metab1.af), abs(hoc.hms.af), abs(ATAC.hoc.af), abs(NET.hoc.af), abs(Ox.genes.af), abs(Ox.metab.af), frame.plot=F, 
        col=c(alpha("#73BF44", 1), alpha("#5FA038", 1), alpha("#FFCD2E", 1), alpha("#FCB216", 1), 
              alpha("#E7A423", 1), alpha("#6D91CB", 1), alpha("#5C7AAA", 1), alpha("#49628B", 1), 
              alpha("#B46AAB", 1), alpha("#975890", 1), alpha("#80497A", 1), alpha("#F2634B", 1), alpha("#CE4F38", 1), 
              "#73BF44", "#5FA038", "#FFCD2E", "#FCB216", "#E7A423", "#6D91CB", "#5C7AAA", "#49628B", 
              "#B46AAB", "#975890", "#80497A", "#F2634B", "#CE4F38"))
dev.off()
calc.cor <- function(x){
  x.mean <- apply(t(x), 2, mean)
  return(cor(x, x.mean)[,1])
}
