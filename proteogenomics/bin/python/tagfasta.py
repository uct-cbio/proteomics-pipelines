#!/usr/bin/env python

import sys
import pandas as pd
import tagmatch
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pickle
import sqlite3

infile = sys.argv[1]
outfile = sys.argv[2]

def calcmw(df):
    charge = df['Identification Charge'] 
    mz = df['Measured m/z'] 
    mw = tagmatch.mz2mw(mz, charge)
    return mw

def create_fasta(df):
    sequence = df['Longest AminoAcid sequence']
    filename = df['File Name']
    scan = df['_scan']
    #mw = df['Theoretic m/z']
    mw = df['_calc.MW']
    ngap = df['N-Gap']
    cgap = df['C-Gap']
    seq = Seq(sequence)
    id = '{};scan={};mw={};ngap={};cgap={}'.format(filename, scan, mw, ngap, cgap)
    #scanid = '{};scan={}'.format(filename, scan)
    rec = SeqRecord(id=id, seq=seq)
    return rec.format('fasta')

tags = pd.read_csv(infile, sep='\t', nrows=100)
columns = tags.columns.tolist()
assert(columns[-1] == 'Identification Charge')
tags = tags.reset_index()

assert len(tags['Identification Charge'].dropna()) == 0
del tags['Identification Charge']
tags.columns = columns

tags['_calc.MW'] = tags.apply(calcmw, axis=1)
tags['_scan'] = tags['Spectrum Title'].apply(lambda x : x.split('scan=')[1].split('"')[0])

tags['_record'] = tags.apply(create_fasta, axis=1)


print(tags.head())
recs = ''.join(tags['_record'])

w = open(outfile, 'w')
w.write(recs)
w.close()
