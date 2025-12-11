################################################################################
## Plot the metabolite LVs to show the shifts between early and late


load(file = 'plspm_model2_results_240507.RData')
load('Data/plspm_datasets_lvs.RData')

pdf('../YMC_paper/metab_profiles_plspm_phased.pdf', height=4, width=5)
par(mfrow=c(2,3))

plot(1, type="n", xlab="", ylab="", xlim=c(0, 30), ylim=c(-2.5, 3), frame.plot=F)
apply(metab.ox1, 2, function(x){lines(c(1:30),x,col='tomato3', lwd=1.5)})

plot(1, type="n", xlab="", ylab="", xlim=c(0, 30), ylim=c(-2.5, 3), frame.plot=F)
apply(metab.rc1, 2, function(x){lines(c(1:30),x,col='blue3', lwd=1.5)})

plot(1, type="n", xlab="", ylab="", xlim=c(0, 30), ylim=c(-2.5, 3), frame.plot=F)
apply(metab.rb, 2, function(x){lines(c(1:30),x,col='green2', lwd=1.5)})

plot(1, type="n", xlab="", ylab="", xlim=c(0, 30), ylim=c(-2.5, 3), frame.plot=F)
apply(metab.ox2, 2, function(x){lines(c(1:30),x,col='tomato3', lwd=1.5)})

plot(1, type="n", xlab="", ylab="", xlim=c(0, 30), ylim=c(-2.5, 3), frame.plot=F)
apply(metab.rc2, 2, function(x){lines(c(1:30),x,col='blue3', lwd=1.5)})

dev.off()
