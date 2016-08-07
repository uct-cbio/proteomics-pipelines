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
try:
    os.mkdir(output +'/tables')
except:
    shutil.rmtree(output +'/tables')
    os.mkdir(output +'/tables')

# Frameshifts
frams = combined[combined['_frameshift'].notnull()]
fs_anno = frams[(((frams['_annotated.frameshift.peptides'].notnull())|(frams['_combined.specific.novel.peptides']=='')) & (frams['_reference.proteins.mapped.count']==1))]
fs_novel = frams[~(((frams['_annotated.frameshift.peptides'].notnull())|(frams['_combined.specific.novel.peptides']=='')) & (frams['_reference.proteins.mapped.count']==1))]
fs_anno.to_csv(output +'/tables/annotated_ICDSs.csv')
fs_novel.to_csv(output +'/tables/novel_ICDSs.csv')

# Mutants
mutants = combined[combined['_exclusive.peptide.strains'].notnull()]
mutants.to_csv(output +'/tables/mutants.csv')

