#!/usr/bin/env python3

import sys
import importlib.machinery
from collections import defaultdict
import yaml

config = yaml.load(open(sys.argv[1]).read())
output = sys.argv[2]

vals=[]

d = '#!/usr/bin/env R\n\n'
vals.append(d)

d = "cols <- c("
vals.append(d)

samples = config['samples']
experiment = defaultdict(list)
strain_group = {}
for sample in samples:
    group = samples[sample]['STRAIN_GROUP']
    strain_group[group] = samples[sample]['STRAIN']
    experiment[group].append(sample)

groups = []
reps = []
for group in experiment:
    groups.append(group)
    for sample in experiment[group]:
        rep = "'iBAQ.{}'".format(sample)
        reps.append(rep)
reps = ','.join(reps)
d = reps
vals.append(d)

vals.append(')\n')

d = 'f <- factor(c('
vals.append(d)

reps = []
for group in groups:
    rep = "rep('{}',{})".format(group, len(experiment[group]))
    reps.append(rep)
d = ','.join(reps)
vals.append(d)

d = "),\n"
vals.append(d)

levels = ','.join(["'" +group + "'" for group in groups])
d ="levels=c({}))\n".format(levels)
vals.append(d)

d='design <- model.matrix(~0+f)\n'
vals.append(d)

d = 'colnames(design) <- c({})\n'.format(levels)
vals.append(d)

d = 'contrast.matrix <- makeContrasts(\n'
vals.append(d)

inter_comps = []
intra_comps = [] 
comps = []

for comparison in config['comparisons']:
    strains = []
    strains.append(strain_group[comparison[0]])
    strains.append(strain_group[comparison[1]])
    strains = list(set(strains))
    comp = '"{}-{}"'.format(comparison[0], comparison[1])
    comps.append(comp)
    if len(strains) == 1:
        intra_comps.append(comp)
    elif len(strains) == 2:
        inter_comps.append(comp)

d = ','.join(comps)
vals.append(d)
d = ',levels=design)\n'
vals.append(d)

intra_comps  = ','.join(intra_comps)
d = 'intra_comps <- c({})\n'.format(intra_comps)
vals.append(d)

inter_comps  = ','.join(inter_comps)
d = 'inter_comps <- c({})\n'.format(inter_comps)
vals.append(d)

template = ''.join(vals)

w = open(output +'/experimental_design.R','w')
w.write(template)
w.close()

