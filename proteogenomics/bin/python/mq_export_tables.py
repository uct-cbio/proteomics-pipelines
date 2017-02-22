#!/usr/bin/env python3

import numpy as np
import pandas as pd
import importlib.machinery
import sys
import os
import shutil
import yaml


loader = importlib.machinery.SourceFileLoader('config', sys.argv[1])
config = yaml.load(open(sys.argv[1]).read())
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

for strain in config['strains']:
    cols= [ 'Identifier',
            "All peptides strain {}".format(strain),
            "Specific peptides strain {}".format(strain),
            "Exclusive peptides strain {}".format(strain),
            "Exclusive peptides reference BLAST strain {}".format(strain),
            "Exclusive peptides polymorphism reference feature overlap strain {}".format(strain),
            "Novel specific peptides strain {}".format(strain),
            "Annotated specific peptides strain {}".format(strain),
            "Non-genomic peptides strain {}".format(strain),
            "Translated orfs strain {}".format(strain),
            'Variant orfs strain {}'.format(strain),
            'Annotation type strain {}'.format(strain),
            'Frameshift validated strain {}'.format(strain),
             'Frameshift evidence strain {}'.format(strain),
             "Best orf-reference blast evalue strain {}".format(strain),
             "Best orf-reference blast match {}".format(strain),
            "Reference BLAST strain {}".format(strain)]
    strain_dir = output +'/strains/' + strain
    strain_cols = [col for col in cols if col in combined.columns]
    strain_df = combined[strain_cols]
    for type in types:
        type_df = strain_df[strain_df['Annotation type strain {}'.format(strain)].apply(lambda x: filter_annotation_type(x, type)) == True]
        type = '_'.join(type.split())
        type = ''.join(''.join(type.split('(')).split(')'))
        type = type.lower()
        type_df.to_csv(strain_dir + '/' + type + '.csv')
        annotations.loc[type, strain ]  = len(type_df)

annotations.to_csv(output + '/tables/annotations.csv')
combined = combined[spescols]

try:
    fs = combined[combined['Frameshift'].notnull()]
    fs.to_csv(output + '/tables/frameshifts.csv')
except:
    pass

try:
    var = combined[combined['Exclusive peptide strains'].notnull()]
    var.to_csv(output +'/tables/variants.csv')
except:
    pass
