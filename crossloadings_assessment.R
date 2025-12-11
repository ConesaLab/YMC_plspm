'''
Study the cross-loadings
'''

load(file = 'plspm_model2_results_240507.RData')

cross.loads.matrix <- matrix(0,nrow=13, ncol=13)
rownames(cross.loads.matrix) <- unique(rel.fitted.plspm$crossloadings$block)
colnames(cross.loads.matrix) <- unique(rel.fitted.plspm$crossloadings$block)

for (lv in unique(rel.fitted.plspm$crossloadings$block)){
  cross.loads <- rel.fitted.plspm$crossloadings[rel.fitted.plspm$crossloadings$block==lv,c(3:15)]
  rownames(cross.loads) <- rel.fitted.plspm$crossloadings[rel.fitted.plspm$crossloadings$block==lv,1]
  cross.loads.table <- table(apply(cross.loads, 1, which.max))
  cross.loads.matrix[as.numeric(names(cross.loads.table)), lv] <- cross.loads.table
}

par(mar=c(8,3,3,3))

barplot(t(t(cross.loads.matrix)/colSums(cross.loads.matrix)), las=2, col = inferno(13), 
        names = c('Rb Genes','Rb Metabs','LOC HMs','LOC ATAC','LOC NET','RC Genes','Late Rc metabs','Early Rc Metabs',
                  'HOC HMs','HOC ATAC','HOC NET','OX Genes','Ox Metabs'))


legend(1,1,c('Rb Genes','Rb Metabs','LOC HMs','LOC ATAC','LOC NET','RC Genes','Late Rc metabs','Early Rc Metabs',
         'HOC HMs','HOC ATAC','HOC NET','OX Genes','Ox Metabs'),col = inferno(13), fill=inferno(13), cex=3)


non.crossloading <- rownames(cross.loads[apply(cross.loads, 1, which.max) == 13,])

write.table(non.crossloading, '~/Desktop/example.txt', quote = F, col.names = F, row.names = F, sep = '\t')
