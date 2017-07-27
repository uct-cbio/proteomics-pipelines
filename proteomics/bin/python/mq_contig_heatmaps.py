#!/usr/bin/env python

import pandas as pd
import sys
import multiprocessing
import numpy as np
import pandas as pd
import sys
import importlib.machinery
import Bio; from Bio import SeqIO
import sequtils
import shutil
import algo
import os
import subprocess
from collections import defaultdict
import json
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Seq import translate
import sys
from collections import Counter
from io import StringIO
import uniprot
import pickle
from io import StringIO
from collections import defaultdict

import plotly.plotly as py
import plotly.graph_objs as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot

loader = importlib.machinery.SourceFileLoader('config', sys.argv[1])
config = loader.load_module()
output = sys.argv[2]

pg = pd.read_csv(output + '/diff/msnbase/cleaned.csv')
intensity_cols = [i for i in pg.columns if 'iBAQ.' in i]

def data_export(df):

    global contig_dict
    global strain
    global intensity_columns

    intensity_vals = []
    for i in intensity_cols:
        intensity_vals.append(df[i])
    

    refmapped =  df['X_reference.entries.mapped']

    try:
        orfs = df['X_translated.orfs.strain.{}'.format(strain)].split('\n')
        orfs = [i for i in orfs if i.startswith('>')]
    except:
        orfs = []
    
    for orf in orfs: 
        contig = orf.split('_')[1]
        orf_start = int(orf.split('|')[-1].split(')')[1].split(':')[0])
        orf_end = int(orf.split('|')[-1].split(')')[1].split(':')[1].split()[0])
        orf_strand = orf.split('|')[-1].split('(')[1].split(')')[0]
        
        if not contig in contig_dict:
            contig_dict[contig] = {}
        if not orf_strand in contig_dict[contig]:
            contig_dict[contig][orf_strand] = defaultdict()
        
        orflabel= '{} {} {}'.format(orf.split('|')[1], orf.split('|')[2].split()[0], refmapped)
        
        contig_dict[contig][orf_strand][orf_start] = (orflabel, intensity_vals)

    
    
for strain in config.strains:
    d = pg.copy()

    contig_dict = {}

    d.apply(data_export, axis=1)
    
    if not os.path.exists( output +'/{}/contig_heatmaps'.format(strain)):
        os.mkdir(output + '/{}/contig_heatmaps'.format(strain))

    for contig in contig_dict:
        
        sample_dict = {}

        for strand in contig_dict[contig]:
            vals = contig_dict[contig][strand]
            intensities = []
            ids = []
            positions = list(vals.keys())

            positions.sort()
            positions = positions[::-1]
            for position in positions:
                tup = vals[position]
                
                ids.append(tup[0])
                #coords = tup[0].split()[0].split('|')[-1].split(')')[1].split(':')
                #coords = (int(coords[0]), int(coords[1]))

                intensities.append(tup[1])
            
            data = [go.Heatmap( z=intensities, x = intensity_cols, y = ids)]
            x = plot(data,filename='{}/{}/contig_heatmaps/{}_({})_strand.html'.format(output,strain, contig, strand), auto_open=False)


