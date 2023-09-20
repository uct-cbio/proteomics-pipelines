#!/usr/bin/env python3

import numpy as np
import pandas as pd
import importlib.machinery
import sys
import os
import shutil
import yaml


#loader = importlib.machinery.SourceFileLoader('config', sys.argv[1])

config = yaml.load(open(sys.argv[1]).read(),Loader=yaml.Loader)

output = sys.argv[2]

for reference in config['reference']:
    combined = pd.read_csv(output + '/{}_combined.csv'.format(reference), sep='\t')
    outdir = output + '/tables/' + reference
    
    
    
    try:
        os.makedirs(outdir)
    except:
        shutil.rmtree(outdir)
        os.makedirs(outdir)
    
    spescols = [i for i in combined.columns if i.startswith('_')]




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
        mult = combined[combined['Specific peptides strain {}'.format(strain)].apply(lambda x : len(str(x).split('\n')) ) > 1]
        novel = len(mult[mult['Reference BLAST strain {}'.format(strain)].isnull()])
        annotations.loc['novel', strain ]  = novel
        
        for type in types:
            type_df = strain_df[strain_df['Annotation type strain {}'.format(strain)].apply(lambda x: filter_annotation_type(x, type)) == True]
            type = '_'.join(type.split())
            type = ''.join(''.join(type.split('(')).split(')'))
            type = type.lower()
            type_df.to_csv(outdir + '/' + type + '.csv')
            annotations.loc[type, strain ]  = len(type_df)

    annotations.to_csv(outdir + '/annotations.csv')
    combined = combined[spescols]

    try:
        fs = combined[combined['Frameshift'].notnull()]
        fs.to_csv(outdir + '/frameshifts.csv')
    except:
        pass

    try:
        var = combined[combined['Exclusive peptide strains'].notnull()]
        var.to_csv(outdir +'/variants.csv')
    except:
        pass
