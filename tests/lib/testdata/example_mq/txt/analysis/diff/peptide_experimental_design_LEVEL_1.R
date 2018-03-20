#!/usr/bin/env R

cols <- c('Intensity.E4','Intensity.E5','Intensity.E8','Intensity.E3','Intensity.E6','Intensity.E2','Intensity.E1','Intensity.E7','Intensity.E15','Intensity.E10','Intensity.E9','Intensity.E14','Intensity.E11','Intensity.E16','Intensity.E12','Intensity.E13')
f <- factor(c(rep('Before',8),rep('After',8)),
levels=c('Before','After'))
design <- model.matrix(~0+f)
colnames(design) <- c('Before','After')
contrast.matrix <- makeContrasts(
"After-Before",levels=design)
