#!/usr/bin/env R

cols <- c('iBAQ.E15','iBAQ.E14','iBAQ.E16','iBAQ.E13','iBAQ.E4','iBAQ.E3','iBAQ.E2','iBAQ.E1','iBAQ.E5','iBAQ.E8','iBAQ.E6','iBAQ.E7','iBAQ.E10','iBAQ.E9','iBAQ.E11','iBAQ.E12')
f <- factor(c(rep('After_B',4),rep('Before_A',4),rep('Before_B',4),rep('After_A',4)),
levels=c('After_B','Before_A','Before_B','After_A'))
design <- model.matrix(~0+f)
colnames(design) <- c('After_B','Before_A','Before_B','After_A')
contrast.matrix <- makeContrasts(
"After_A-Before_A","After_B-Before_B",levels=design)
