#!/usr/bin/env python

import pandas as pd
import Bio
from collections import defaultdict
import sys
from Bio import SeqIO
import sequtils
from Bio.Align import AlignInfo
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
import tempfile
from Bio.Seq import Seq
import shutil
import multiprocessing

fasta=sys.argv[1]
out=sys.argv[2]
method=sys.argv[3]

methods=['multiple', 'consensus', 'concatenate']
assert method in methods

inpt = list(SeqIO.parse(fasta,'fasta'))

records = defaultdict(list)
for rec in inpt:
    desc = ';'.join(rec.description.split(';')[:-2])
    desc = '.'.join(desc.split())
    records[desc].append(rec)

newrecords = []

keys =list(records.keys())

def process(key):
    options = records[key]
    count = 1
    for option in options:
        option.id = 'record_'+ str(count)
        count += 1
    seqs = [str(i.seq) for i in options]
    
    
    if method == 'consensus':
        
        if len(seqs) > 1:
            
            temp=tempfile.mkdtemp()
            aln, tree= sequtils.clustalw(temp + '/temp.fasta', options)
            summary_align = AlignInfo.SummaryInfo(aln)
            consensus = summary_align.dumb_consensus()
            peptide = str(consensus)
            shutil.rmtree(temp)
        else:
            peptide =  seqs[0]
    
    elif method =='concatenated':
        peptide = '-'.join(list(set(seqs)))

    new_record = SeqRecord(id= key, seq=Seq(peptide))
    return new_record

if method !='multiple':
    pool = multiprocessing.Pool(56)
    newrecords = pool.map(process, keys)
else:
    newrecords=[]
    for scan in records:
        for option in records[scan]:
            new_record = SeqRecord(id= scan, seq=option.seq)
            newrecords.append(new_record)
SeqIO.write(newrecords, out, 'fasta')

