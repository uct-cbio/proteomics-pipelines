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

db = SeqIO.parse(sys.argv[1],'fasta')
print('loaded db')

peprecs = SeqIO.parse(sys.argv[2], 'fasta')

pepdict = {}
for peprec in peprecs:
    pepdict[str(peprec.seq)]=peprec

print('created peptide dict')

Trie = algo.Trie(list(pepdict.keys()))
print('created trie')

fragtol = float(sys.argv[3])
out = sys.argv[4] + '/' + sys.argv[1].split('/')[-1] +'.csv'
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
        trec = pepdict[pepstr]
        tseq = str(trec.seq)
        ids = trec.description.split('scans=')[1].split('|')
        prec_ions = defaultdict(set)
        
        for id in ids:
            _ = id.split(';')
            mz = float(_[-2])
            charge = int(_[-1][0])
            mw = tagmatch.mz2mw(mz, charge)
            prec_ions[mw].add(id)
        
        prec_array = list(prec_ions.keys())
        tm  = tagmatch.TagMatch(pepstr, prec_array, qrec, tol = fragtol)
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
