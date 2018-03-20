#!/usr/bin/env R

cols <- c('iBAQ.E4','iBAQ.E5','iBAQ.E8','iBAQ.E3','iBAQ.E6','iBAQ.E2','iBAQ.E1','iBAQ.E7','iBAQ.E15','iBAQ.E10','iBAQ.E9','iBAQ.E14','iBAQ.E11','iBAQ.E16','iBAQ.E12','iBAQ.E13')
f <- factor(c(rep('Before',8),rep('After',8)),
levels=c('Before','After'))
design <- model.matrix(~0+f)
colnames(design) <- c('Before','After')
contrast.matrix <- makeContrasts(
"After-Before",levels=design)
