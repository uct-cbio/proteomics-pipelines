#!/usr/bin/env python

import collections
from collections import defaultdict
import pandas as pd
import sys
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import algo
import re
import tagmatch
import proteomics
import pickle
import numpy as np

db = SeqIO.parse(sys.argv[1],'fasta')
print('loaded db')

pepdict = pickle.load(open(sys.argv[2], 'rb'))
print('loaded peptide dict')

#Trie = algo.Trie(list(pepdict.keys()))

prectol = float(sys.argv[3])

specificity=sys.argv[4]

enzymes = [i for i in sys.argv[5].split(';') if i != '']
max_missed_cleavages = int(sys.argv[6])
gap_tol = float(sys.argv[7])
fixed_modifications = [i for i in sys.argv[8].split(';') if i != '']
variable_modifications = [i for i in sys.argv[9].split(';') if i != '']


Trie = pickle.load(open(sys.argv[10],'rb'))
print('loaded trie')

out = sys.argv[11] + '/' + sys.argv[1].split('/')[-1] +'.csv'
cv = []
scounts = []
length = []
peptides = []
scans = []
num_peptides = []
fasta = []
accession = []
samples = [] 

for rec in db:
    qrec = str(rec.seq)
    rec_peps = set()
    rec_pep_ids = set()
    TrieMatch=algo.TrieMatch(Trie, qrec)
    peptides_mapped=TrieMatch.trie_export()
    covered = []

    for pepstr in peptides_mapped:
        #trec = pepdict[pepstr]
        ids = pepdict[pepstr]
        #tseq = str(trec.seq)
        tseq = pepstr
        #ids = trec.description.split('scans=')[1].split('|')
        prec_ions = defaultdict(set)
        
        for id in ids:
            _ = id.split(';')
            scan = _[0] + ';' + _[1]
            mw = np.round(float(_[2].split('mw=')[1]),6)
            ngap = np.round(float(_[3].split('ngap=')[1]),6)
            cgap = np.round(float(_[4].split('cgap=')[1]),6)
            prec_ions[(ngap, mw, cgap)].add(scan)
        
        prec_array = list(prec_ions.keys())
        print(pepstr, prec_array, qrec, prectol, specificity)
        tm  = tagmatch.TagMatch(pepstr, prec_array, qrec, prec_tol = prectol, specificity=specificity, enzymes=enzymes, max_missed_cleavages=max_missed_cleavages, gap_tol=gap_tol , fixed_modifications=fixed_modifications, variable_modifications=variable_modifications)
        rec_peps.update(tm.validated_peptides)
        
        for _ in tm.validated_precursor:
            validated_scans = prec_ions[_]
            rec_pep_ids.update(validated_scans)

    ln = len(qrec)
    sample_spectra = defaultdict(list)
    for id in rec_pep_ids:
        sample = id.split(';')[0]
        spectrum = id
        sample_spectra[sample].append(spectrum)

    #cov = TrieMatch.trie_coverage()
    sc = len(rec_pep_ids)
    if len(rec_pep_ids) > 0: 
        samples.append(sample_spectra.copy())
        peptides.append('\n'.join(rec_peps))
        num_peptides.append(len(rec_peps))
        scans.append('\n'.join(rec_pep_ids))
        scounts.append(sc)
        #cv.append(cov)
        accession.append(rec.id)
        fasta.append(rec.format('fasta'))
        length.append(ln)
  
table = pd.DataFrame()

table['ID'] = pd.Series(accession)
table['Records'] = pd.Series(fasta)
table['SpectralCounts'] = pd.Series(scounts)
table['MSMS'] = pd.Series(scans)
table['PeptideSequences'] = pd.Series(peptides)
table['PeptideCounts'] = pd.Series(num_peptides)
#table['SequenceCoverage_%'] = pd.Series(cv)
table['SequenceLength'] = pd.Series(length)

table = table.reset_index()

del table['index']

for row in table.iterrows():
    indx = row[0]
    sample_data = samples[indx]
    for sample in sample_data:
        table.loc[indx, 'Sample ' + sample + ' (msms)'] = len(sample_data[sample])

table.to_csv(out)
