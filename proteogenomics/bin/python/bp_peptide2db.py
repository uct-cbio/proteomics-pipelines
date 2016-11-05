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

out = sys.argv[3] + '/' + sys.argv[1].split('/')[-1] +'.csv'
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
    rec_peps = []
    rec_pep_ids = set()
    peptides_mapped=Trie.trie_export(qrec)
    covered = []
    print(qrec, peptides_mapped)
    for pepstr in peptides_mapped:
        trec = pepdict[pepstr]
        tseq = str(trec.seq)
        rec_peps.append(tseq)
        ids = trec.description.split('scans=')[1].split('|')
        rec_pep_ids.update(ids)
    ln = len(qrec)
    sample_spectra = defaultdict(list)
    for id in rec_pep_ids:
        sample = id.split(';')[0]
        sample_spectra[sample].append(id)
    cov = Trie.trie_coverage(qrec)
    sc = len(rec_pep_ids)
    if cov > 0: 
        samples.append(sample_spectra.copy())
        peptides.append('\n'.join(rec_peps))
        num_peptides.append(len(rec_peps))
        scans.append('\n'.join(rec_pep_ids))
        scounts.append(sc)
        cv.append(cov)
        accession.append(rec.id)
        fasta.append(rec.format('fasta'))
        length.append(ln)
    
table = pd.DataFrame()

table['ID'] = pd.Series(accession)
table['Records'] = pd.Series(fasta)
table['SpectralCounts'] = pd.Series(scounts)
table['MSMS'] = pd.Series(scans)
table['SubjectSequences'] = pd.Series(peptides)
table['SubjectCounts'] = pd.Series(num_peptides)
table['SequenceCoverage_%'] = pd.Series(cv)
table['SequenceLength'] = pd.Series(length)

table = table.reset_index()
del table['index']

for row in table.iterrows():
    indx = row[0]
    sample_data = samples[indx]
    for sample in sample_data:
        table.loc[indx, 'Sample ' + sample + ' (msms)'] = len(sample_data[sample])

table.to_csv(out)
