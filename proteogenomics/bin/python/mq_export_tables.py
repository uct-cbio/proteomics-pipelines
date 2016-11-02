#!/usr/bin/env python

import numpy as np
import pandas as pd
import importlib.machinery
import sys
import os
import shutil



loader = importlib.machinery.SourceFileLoader('config', sys.argv[1])
config = loader.load_module()
output = sys.argv[2]

combined = pd.read_csv(output + '/combined.csv')

spescols = [i for i in combined.columns if i.startswith('_')]

try:
    os.mkdir(output +'/tables')
except:
    shutil.rmtree(output +'/tables')
    os.mkdir(output + '/tables')

def filter_annotation_type(val, type):
    try:
        if type in val:
            return True
        else:
            return False
    except:
        return False

types = ['Putative upstream TSS (upstream non-TSS peptide)', 'Downstream TSS identified', 'Annotated TSS validated', 'No putative upstream TSS (upstream non-TSS peptide)', 'Annotated start site non-enzymative cleavage site', 'Upstream non-enzymative cleavage site', 'Upstream TSS identified', 'Downstream non-enzymative cleavage site']

annotations = pd.DataFrame()

for strain in config.strains:
    cols= [ 'Identifier',
            '_all.peptides.strain.{}'.format(strain),
            '_specific.peptides.strain.{}'.format(strain),
            '_exclusive.peptides.strain.{}'.format(strain),
            '_exclusive.peptides.mapped.reference.blast.strain.{}'.format(strain),
            '_protein.group.reference.match',
            '_identified.polymorphism.mapped.reference.feature.overlap.strain.{}'.format(strain),
            '_orfs.mapped.frameshift.validated.strain.{}'.format(strain),
            '_orfs.mapped.frameshift.evidence.strain.{}'.format(strain),
            '_unmapped.peptides.strain.{}'.format(strain),
            '_translated.orfs.strain.{}'.format(strain),
            '_variant.orfs.strain.{}'.format(strain),
            '_annotation.type.strain.{}'.format(strain),
            '_reference.entries.mapped',
            '_translated.orfs.mapped.reference.best.blast.evalue.strain.{}'.format(strain),
            '_translated.orfs.mapped.reference.best.blast.match.{}'.format(strain) ]
    strain_dir = output +'/' + strain
    strain_df = combined[cols]
    for type in types:
        type_df = strain_df[strain_df['_annotation.type.strain.{}'.format(strain)].apply(lambda x: filter_annotation_type(x, type)) == True]
        type = '_'.join(type.split())
        type = ''.join(''.join(type.split('(')).split(')'))
        type = type.lower()
        type_df.to_csv(strain_dir + '/' + type)
        annotations.loc[type, strain ]  = len(type_df)


annotations.to_csv(output + '/tables/annotations.csv')
combined = combined[spescols]

fs = combined[combined['_frameshift'].notnull()]
fs.to_csv(output + '/tables/frameshifts.csv')

var = combined[combined['_exclusive.peptide.strains'].notnull()]
var.to_csv(output +'/tables/variants.csv')


