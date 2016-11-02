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

db = list(SeqIO.parse(sys.argv[1],'fasta'))

peprecs = list(SeqIO.parse(sys.argv[2], 'fasta'))

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
    
    covered = []
    for trec in peprecs:
        tseq = str(trec.seq)
        if tseq in qrec:
            rec_peps.append(tseq)
            ids = trec.description.split('scans=')[1].split('|')
            rec_pep_ids.update(ids)

    ln = len(qrec)
    
    sample_spectra = defaultdict(list)
    for id in rec_pep_ids:
        sample = id.split(';')[0]
        sample_spectra[sample].append(id)
    
    for rp in rec_peps:
        starts = [m.start() for m in re.finditer('(?={})'.format(rp), qrec)]
        for s in starts:
            for _ in range(s, s + len(rp)):
                covered.append(_)
    
    cov = len(list(set(covered)))/float(ln) * 100
    sc = len(rec_pep_ids)
    #empai = proteomics.emPAI(sc, qrec,  min_len, max_len)

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
