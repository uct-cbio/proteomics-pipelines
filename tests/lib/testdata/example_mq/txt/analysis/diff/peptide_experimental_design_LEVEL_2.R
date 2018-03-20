#!/usr/bin/env R

cols <- c('Intensity.E15','Intensity.E14','Intensity.E16','Intensity.E13','Intensity.E4','Intensity.E3','Intensity.E2','Intensity.E1','Intensity.E5','Intensity.E8','Intensity.E6','Intensity.E7','Intensity.E10','Intensity.E9','Intensity.E11','Intensity.E12')
f <- factor(c(rep('After_B',4),rep('Before_A',4),rep('Before_B',4),rep('After_A',4)),
levels=c('After_B','Before_A','Before_B','After_A'))
design <- model.matrix(~0+f)
colnames(design) <- c('After_B','Before_A','Before_B','After_A')
contrast.matrix <- makeContrasts(
"After_A-Before_A","After_B-Before_B",levels=design)
